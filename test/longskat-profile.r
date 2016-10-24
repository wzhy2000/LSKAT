longskat_gene_plink_profile<-function( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, file.phe.long,  file.phe.cov, file.phe.time=NULL,
			gene.names, options=list() )
{
	cat( "[ LONGSKAT_GENE_PLINK ] Procedure.\n");
	cat( "Checking the optional items......\n");

	if (missing(options))
		options <- get_default_options()
	else
	{
		options0 <- get_default_options();
        options0[names(options)] <- options;
        options <- options0;
    }

	cat( "* Covariate Time Effect: ",  options$time.cov, "\n");
	cat( "* Parallel Computing: ", ifelse( options$n.cpu>1, "Yes,", "No,"), options$n.cpu,  "CPU(s)\n");
	cat( "* Debug Output: ", ifelse( options$verbose, "Yes", "No"),"\n");
	cat( "* SNP Impute: ", options$snp.impute,"\n");
	cat( "* C/C++ Module Used Output: ", ifelse( options$run.cpp, "Yes", "No"), "\n");
	cat( "* Beta Weights for Common SNPs: ",  options$w.common[1], options$w.common[2], "\n");
	cat( "* Beta Weights for Rare SNPs: ",  options$w.rare[1], options$w.rare[2], "\n");

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	chk.phe <- check_pheno_file( file.phe.long, file.phe.time, chk.plink$family )
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, chk.plink$family )
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	cat( "Starting to load all data files......\n");

	PF <- read_gen_phe_cov ( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time, file.phe.cov );
	PF.pars <- list(file.plink.bed = file.plink.bed,
			file.plink.bim = file.plink.bim,
			file.plink.fam = file.plink.fam,
			file.phe.long  = file.phe.long,
			file.phe.time  = file.phe.time,
			file.phe.cov   = file.phe.cov,
			file.gene.set  = file.gene.set,
			cov.count    = NCOL(PF$phe.cov)-1,
			time.cov     = options$time.cov,
			w.common       = options$w.common,
			w.rare         = options$w.rare );

	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	r.model <- longskat_est_model(PF$phe.long, PF$phe.cov, PF$phe.time, time.cov=options$time.cov, g.maxiter=options$g.maxiter, verbose=options$verbose);

	if( class(r.model) != "LSKAT.null.model" )
		stop("! Failed to estimate the parameters of Covariance Compoment.");

	cat("* SIGMA_A =", h0$sig.a, "\n");
	cat("* SIGMA_B =", h0$sig.b, "\n");
	cat("* SIGMA_E =", h0$sig.e, "\n");
	cat("* RHO =", h0$rho, "\n");
	cat("* MU =", h0$u, "\n");
	cat("* Beta(Cov) =", h0$par_cov, "\n");
	cat("* Beta(Time) =", h0$par_t, "\n");
	cat("* L(min) =", h0$val, "\n");

	cat("* GENE.NAMES =", length(gene.names), "-", gene.names, "\n");
	cat("Starting the LSKAT estimate for each gene......\n");

	tm.start <- proc.time();

	gen.list <- read_gen_dataset(file.gene.set, file.plink.bim);

	rs.lst <- list();
	rs.lst.idx <- 0;

	for(k in 1:length(gene.names) )
	{
		gen <- get_gen_family ( gen.list, gene.names[k] );
		if (options$verbose) cat(" Finding", length(gen$snps), "SNPs...\n");
		gen.mat <- get_snp_mat( PF$snp.mat, gen, snp.impute );

		if (is.null(gen.mat)) next;

		if (length(which(gen.mat$maf>0.5))>0) browser();

		if (length(gen.mat$maf)==0) next;

		ls <- longskat_gene_test(r.model, gen.mat$snp,
							    weights.common= options$weights.common, weights.rare  = options$weights.rare, run.cpp=options$run.cpp, verbose = options$verbose );

		if(is.null(ls))
		{
			cat("! Failed to calculate Gene[", k, "]=", as.character(gen$name), "\n");
			next;
		}

		rs.lst.idx <- rs.lst.idx + 1;

		rs.name  <- c(as.character(gen$name));
		rs.lskat <- c(0, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv );

		if(options$verbose) cat(" [", k, "]", as.character(gen$name), ls$snp.total, ls$snp.rare, ls$qv, ls$pv, "\n");

		if(NROW(gen.mat$info)>1)
		{
			for(i in 1:NROW(gen.mat$info) )
			{
				snp.mat <- gen.mat$snp[-i,,drop=F]
				ls <- longskat_gene_test(r.model,
						    time.cov = options$time.cov,
						    weights.rare  = options$weights.rare,
						    weights.common= options$weights.common,
						    run.cpp=options$run.cpp,
						    verbose = options$verbose);

				if(is.null(ls))
				{
					cat("! Failed to calculate SNP[", i, "]=", gen.mat$info[i,1], "\n");
					next;
				}

				rs.name  <- c(rs.name, paste("-", gen.mat$info[i,1], sep=""));
				rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
			}

			if(NROW(gen.mat$info)>2)
			{
				pv <- c();
				for(i in 1:NROW(gen.mat$info) )
				{
					snp.mat <- gen.mat$snp[i,,drop=F]
					ls <- longskat_gene_test(Y.delt, X, snp.mat, par_null,
							    time.cov = options$time.cov,
							    weights.rare  = options$weights.rare,
							    weights.common= options$weights.common,
							    run.cpp=options$run.cpp,
							    verbose = options$verbose);

					rs.name  <- c(rs.name, paste("*", gen.mat$info[i,1], sep=""));
					rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					pv <- c(pv, ls$pv);
				}

				r.pv <- sort.int(pv, decreasing=T, index.return=T);
				r.sel <- c();
				for(i in r.pv$ix )
				{
					r.sel <- c(r.sel, i);
					snp.mat <- gen.mat$snp[r.sel,,drop=F]
					ls <- longskat_gene_test(Y.delt, X, snp.mat, par_null,
							    time.cov = options$time.cov,
							    weights.rare  = options$weights.rare,
							    weights.common= options$weights.common,
							    run.cpp=options$run.cpp,
							    verbose = options$verbose);

					rs.name  <- c(rs.name, paste("*", gen.mat$info[r.sel,1], sep="", collapse="&"));
					rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					pv <- c(pv, ls$pv);
				}

			}

			if(0)
			{
				if(NROW(gen.mat$info)>3)
				{
					p <- combn(1:NROW(gen.mat$info), 2);
					for(i in 1:NCOL(p) )
					{
						snp.mat <- gen.mat$snp[ p[,i],, drop=F ]
						ls <- longskat_gene_test(Y.delt, X, snp.mat, par_null,
								    time.cov = options$time.cov,
								    weights.rare  = options$weights.rare,
								    weights.common= options$weights.common,
								    run.cpp=options$run.cpp,
								    verbose = options$verbose);

						rs.name  <- c(rs.name, paste("*", gen.mat$info[p[1,i],1], "&", gen.mat$info[p[2,i],1], sep=""));
						rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					}

					for(i in 1:NCOL(p) )
					{
						snp.mat <- gen.mat$snp[ -p[,i],, drop=F ]
						ls <- longskat_gene_test(Y.delt, X, snp.mat, par_null,
								    time.cov = options$time.cov,
								    weights.rare  = options$weights.rare,
								    weights.common= options$weights.common,
								    run.cpp=options$run.cpp,
								    verbose = options$verbose);

						rs.name  <- c(rs.name, paste("-", gen.mat$info[p[1,i],1], "&", gen.mat$info[p[2,i],1], sep=""));
						rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					}
				}
			}

			rownames(rs.lskat)<-rs.name;
			colnames(rs.lskat)<-c("idx", "chr", "min.pos", "snp", "rare", "Q", "pv");
		}

		rs.lst[[rs.lst.idx]] <- rs.lskat;
	}

	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");

	return(rs.lst);
}
