#public:
longskat_gene_test <- function( r.model, snp.mat, weights.common=c(0.5,0.5), weights.rare=c(1,25),
					rare.cutoff=NULL, test.type="Joint", snp.impute="mean", verbose=F )
{
	run.cpp <- FALSE;

	if(!(tolower(snp.impute) %in% c("random", "mean")))
		stop("The option 'snp.impute' has 2 values: random, mean.");

	if(!(tolower(test.type) %in% c("joint", "common.only", "rare.only")))
		stop("The option 'test.type' has 3 values: Joint, Common.Only, Rare.Only");

	Y.delt   <- r.model$phe.delt;
	Y.time   <- r.model$phe.time;
	par_null <- c( r.model$par$sig.a,  r.model$par$sig.b, r.model$par$sig.e, r.model$par$rho );

	idx.snp <- match(rownames(Y.delt), rownames(snp.mat))
	if(length(which(is.na(idx.snp)))>0)
		warning("Unconsistent Individual set between SNP matrix and NULL Model.");

	snp.mat <- snp.mat[idx.snp,,drop=F];
	snp.NMISS <- unlist(apply(snp.mat, 1, function(snp){length(which(is.na(snp)))}));
	
	if (any(is.na(colMeans(snp.mat)/2)))
		snp.mat <- snp_impute( snp.mat, impute=snp.impute);

	X <- r.model$phe.cov;
	if(r.model$par$intercept)
		X <- cbind(1, r.model$phe.cov);

	if(is.null(rare.cutoff))
		rare.cutoff <- 1/sqrt(2*length(idx.snp));

	## SKAT Mode
 	Z.scale.skat <- SKAT_Scale_Genotypes( X,
 					snp.mat,
 					weights.common = weights.common,
 					weights.rare   = weights.rare,
 					rare.cutoff    = rare.cutoff,
 					test.type      = test.type,
 					r.corr.common  = 0,
 					r.corr.rare    = 0 )

	## Burden Mode
 	Z.scale.burden <- SKAT_Scale_Genotypes( X,
 					snp.mat,
 					weights.common = weights.common,
 					weights.rare   = weights.rare,
 					rare.cutoff    = rare.cutoff,
 					test.type      = test.type,
 					r.corr.common  = 1,
 					r.corr.rare    = 1 )

	Q <- try( est_gen_Q( par_null,
					X,
 					Y.delt,
					Y.time,
					Z.scale.skat,
					Z.scale.burden,
					r.model$par$time.cov,
					run.cpp = run.cpp ) );

	if (class(Q)=="try-error")
		return(NULL);

	p.lskat <- get_Qv_pvalue(Q$v.lskat, Q$w.lskat);
	p.burden <- get_Qu_pvalue(Q$v.burden, Q$w.burden);

	if(verbose) cat( " MAF=", (Z.scale.skat$maf), "\n");

	r.lskat <- list( mle      = r.model,
					snp.NMISS = snp.NMISS,
					snp.MAF   = colMeans(snp.mat)/2,
					snp.total = NCOL(snp.mat),
					snp.rare  = Z.scale.skat$rare,
					q.lskat   = Q$v.lskat,
					p.lskat   = p.lskat$p.value,
					q.lburden = Q$v.burden,
					p.lburden = p.burden$p.value);
	class(r.lskat) <- "LSKAT.gen.ret";

	return(r.lskat)
}

est_gen_Q.R<-function( par_null, X, Y.delt, Y.time, maf, Z, time.cov )
{
	Get_SKAT_Residuals.Get_X1 = function(X1){

		qr1<-qr(X1)
		q1<-ncol(X1)
		if(qr1$rank < q1){

			X1.svd<-svd(X1)
			X1 = X1.svd$u
		}

		return(X1)
	}

	X <- Get_SKAT_Residuals.Get_X1(X);

	sig.a2  <- par_null[1]^2;
	sig.b2  <- par_null[2]^2;
	sig.e2  <- par_null[3]^2;
	par_rho <- par_null[4];

	n <- NROW(Y.delt);
	m <- NCOL(Y.delt);
	k <- NCOL(Z);

	# SKAT   qv,pv
	Q.v<-0;
	# burden qu,pu
	Q.u<-0;
	Q.w<-c();

	AR.1 <- matrix(rep(1:m,m),ncol=m, byrow=T)
	AR.1 <- par_rho^abs(t(AR.1)-AR.1)

	V    <- array(1, dim=c(m,m)) * sig.a2 +  AR.1 * sig.b2 + diag(1, m) * sig.e2
	V_1  <- solve(V);

	V.i <- lapply(1:n, function(i) { idx<- !is.na(Y.delt[i,]); return(solve(V[idx, idx])); });
	Y.i <- lapply(1:n, function(i) { idx<- !is.na(Y.delt[i,]); return(Y.delt[i,idx]); });
	M.i <- unlist(lapply(1:n, function(i) { return(length(Y.i[[i]])); }));

	for(i in 1:k)
	{
		Q.i <- lapply(1:n, function(j){ Z[j,i]*t(rep(1, M.i[j] )) %*% V.i[[j]] %*% (Y.i[[j]]);});
		Q.v <- Q.v + (sum(unlist(Q.i)))^2;
		Q.u <- Q.u + (sum(unlist(Q.i)));
	}

	Q.v <- Q.v/2;
	Q.u <- Q.u^2/2;

	n.x <- NCOL(X) + time.cov ;
	W0 <- array(0, dim=c(k,k));
	W1 <- array(0, dim=c(k,n.x));
	W2 <- array(0, dim=c(n.x,n.x));
	W3 <- array(0, dim=c(n.x,k));

	for (i in 1:n)
	{
		m.i <- M.i[i] ;
		kr.Z<- kronecker( Z[i,,drop=F], array(1, dim=c( m.i, 1 )) )
		kr.X<- kronecker( X[i,,drop=F], array(1, dim=c( m.i,1)) )

		if( time.cov>0 )
		{
			yt.i <- Y.time[i, !is.na(Y.delt[i,]) ];

			for(t in 1:time.cov)
				kr.X<- cbind( kr.X, yt.i**t);
		}

		W0 <- W0 + t(kr.Z) %*% V.i[[i]] %*% kr.Z;
		W1 <- W1 + t(kr.Z) %*% V.i[[i]] %*% kr.X;
		W2 <- W2 + t(kr.X) %*% V.i[[i]] %*% kr.X;
		W3 <- W3 + t(kr.X) %*% V.i[[i]] %*% kr.Z;
	}

	Q.w <- (W0 - W1 %*% solve(W2) %*% W3)/2;

#cat("v=", Q.v, "u=", Q.u, "\n");
#show( Q.w );

	return(list(v=Q.v, w=Q.w, u=Q.u));
}

est_gen_Q<-function( par_null, X, Y.delt, Y.time, Z.scale.skat, Z.scale.burden, time.cov, run.cpp=T)
{
	t0 <- proc.time()
	r.lskat<- list(v=0, w=0);

	if(!is.null(Z.scale.skat$new))
		if(!run.cpp)
			r.lskat<- est_gen_Q.R(
						par_null,
						X,
						Y.delt,
						Y.time,
						Z.scale.skat$maf,
						Z.scale.skat$new,
						time.cov)
		else
			lskat<- .Call( "est_gen_Q_C",
						as.vector( par_null*1.0),
						as.matrix( X*1.0 ),
						as.matrix( Y.delt*1.0 ),
						as.matrix( Y.time*1.0 ),
						as.vector( Z.scale.skat$maf*1.0),
						as.matrix( Z.scale.skat$new*1.0),
						time.cov);


	r.burden<- list(v=0, w=0);
	if(!is.null(Z.scale.skat$new))
		if(!is.null(Z.scale.burden$new))
			r.burden<- est_gen_Q.R(
						par_null,
						X,
						Y.delt,
						Y.time,
						Z.scale.burden$maf,
						Z.scale.burden$new,
						time.cov)
		else
			r.burden<- est_gen_Q.R(
						as.vector( par_null*1.0),
						as.matrix( X*1.0 ),
						as.matrix( Y.delt*1.0 ),
						as.matrix( Y.time*1.0 ),
						as.vector( Z.scale.burden$maf*1.0),
						time.cov);

	r <- list(  v.lskat = r.lskat$v,
				w.lskat  = r.lskat$w,
				v.burden = r.burden$u,
				w.burden = r.burden$w);

	t1 <- proc.time();

	# r1<- est_gen_Q.R(Y.delt, maf, Z, X, par_null)
	# if( any(c(r1$v) != c(r0$v) ) ) cat("V1<>V0");
	# if( any(c(r1$w) != c(r0$w) ) ) cat("W1<>W0");
	# cat("EST_GEN_Q Runtime:", (t1-t0)[3], "\n");

	return(r);
}

#public:
print.LSKAT.gen.ret<-function(r.lskat, useS4=FALSE)
{
	cat("  LSKAT/Gene Summary\n");
	cat("[1] MLE Results:\n");
	cat("* SIGMA_A =", 		   r.lskat$mle$par$sig.a, "\n");
	cat("* SIGMA_B =",         r.lskat$mle$par$sig.b, "\n");
	cat("* SIGMA_E =",         r.lskat$mle$par$sig.e, "\n");
	cat("* RHO =",             r.lskat$mle$par$rho, "\n");
	cat("* MU =",              r.lskat$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        r.lskat$mle$par$cov.effect, "\n");
	cat("* Beta(Time)=",       r.lskat$mle$par$time.effect, "\n");
	cat("* L(min) =",          r.lskat$mle$likelihood, "\n");

	cat("[2] LSKAT:\n");
	cat("* SNP total = ",      r.lskat$snp.total, "\n");
	cat("* Rare SNPs = ",      r.lskat$snp.rare, "\n");
	cat("* Q @ LSKAT = ",      r.lskat$q.lskat, "\n");
	cat("* p-value   = ",      r.lskat$p.lskat, "\n");
	cat("* Q @ L-burden= ",      r.lskat$q.lburden, "\n");
	cat("* p-value   = ",      r.lskat$p.lburden, "\n");
}

#private
longskat_gene_task<-function( r.model, PF.gen, PF.phe, gene, weights.common, weights.rare, run.cpp, snp.impute, rare.cutoff, test.type="Joint", verbose=FALSE)
{
	rs.name <-c();
	rs.lst <- vector("list", length(which(!is.na(gene))) );
	rs.lst.idx <- 0;

	for(i in gene )
	{
		if (is.na(i)) next;

		gen <- try( get_gen_mat( PF.gen, i, snp.impute ) );
		if( is.null(gen) || class(gen)=="try-error" || length(gen$maf)==0 )
		{
			if (verbose) cat("! No SNPS for Gene[", i, "]=", i, "\n");
			next;
		}
		else
			if (verbose) cat("  Finding", NROW(gen$snp), "SNPs...\n");

		ls <- try( longskat_gene_test( r.model,
						gen$snp,
						weights.common = weights.common,
						weights.rare   = weights.rare,
						#run.cpp        = run.cpp,
						rare.cutoff    = rare.cutoff,
						test.type      = test.type,
						verbose        = verbose) );

		rs.lst.idx <- rs.lst.idx + 1;
		rs.name <- c(rs.name, as.character(gen$name));

		if(is.null(ls) || class(ls) == "try-error" )
		{
			rs.lst[[rs.lst.idx]] <- c( i,
						min(gen$info[,2]),
						min(gen$info[,3]),
						NROW(gen$mat),
						NA, NA, NA, NA, NA );
			cat("! Failed to calculate Gene[", i, "]=", as.character(gen$name), "\n");
		}
		else
		{
			rs.lst[[rs.lst.idx]] <- c(i,
						min(gen$info[,2]),
						min(gen$info[,3]),
						ls$snp.total,
						ls$snp.rare,
						ls$q.lskat,
						ls$p.lskat,
						ls$q.lburden,
						ls$p.lburden );

			if(verbose) cat(" [", i, "]",
						as.character(gen$name),
						ls$snp.total,
						ls$snp.rare,
						ls$q.lskat,
						ls$p.lskat,
						ls$q.lburden,
						ls$p.lburden, "\n");
		}
	}

	rs <- do.call( "rbind", rs.lst );
	ret <- data.frame( index = rs[,1],
					gene.name= rs.name,
					chr      = rs[,2],
					min.pos  = rs[,3],
					snp.total= rs[,4],
					snp.rare = rs[,5],
					q.lskat  = rs[,6],
					p.lskat  = rs[,7],
					q.lburden= rs[,8],
					p.lburden= rs[,9] );

	return(ret);
}

#public:
longskat_gene_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam,
				file.phe.long,
				file.phe.cov,
				file.phe.time=NULL,
				file.gene.set,
				gene.set = NULL,
				options=list(),
				verbose=FALSE)
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

	show_options( options );

	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	chk.phe <- check_pheno_file( file.phe.long, file.phe.time, file.plink.fam )
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, file.plink.fam )
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	cat( "Starting to load all data files......\n");

	PF.gen <- read_gen_plink ( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, options$plink.path )

	PF.phe <- read_phe_cov ( file.phe.long, file.phe.time, 	file.phe.cov, PF.gen );

	PF.par <- list( file.plink.bed = file.plink.bed,
					file.plink.bim = file.plink.bim,
					file.plink.fam = file.plink.fam,
					file.phe.long  = file.phe.long,
					file.phe.time  = file.phe.time,
					file.phe.cov   = file.phe.cov,
					file.gene.set  = file.gene.set,
					cov.count      = NCOL(PF.phe$phe.cov),
					time.cov       = options$time.cov,
					weights.common = options$weights.common,
					weights.rare   = options$weights.rare,
					rare.cutoff    = options$rare.cutoff);


	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	r.model <- longskat_est_model(PF.phe$phe.long,
					PF.phe$phe.cov,
					PF.phe$phe.time,
					time.cov   = options$time.cov,
					g.maxiter  = options$g.maxiter,
					intercept  = options$intercept,
					verbose    = options$verbose,
					par.init   = options$par.init,
					method     = options$est.method);

	if( class(r.model) != "LSKAT.null.model" )
		stop("! Failed to estimate the parameters of Covariance Compoment.");

	cat("* SIGMA_A =",   r.model$par$sig.a, "\n");
	cat("* SIGMA_B =",   r.model$par$sig.b, "\n");
	cat("* SIGMA_E =",   r.model$par$sig.e, "\n");
	cat("* RHO =",       r.model$par$rho, "\n");
	cat("* MU =",        r.model$par$mu, "\n");
	cat("* Beta(Cov) =", r.model$par$cov.effect, "\n");
	cat("* Beta(Time) =",r.model$par$time.effect, "\n");
	cat("* L(min) =",    r.model$likelihood, "\n");
	cat("* Invalid IDs =",    length(r.model$na.rm.id), "\n");

	if(is.null(PF.gen$gen.list))
		stop("! Failed to load gene definition file.");

	if (is.null(gene.set))
		gene.set<- c(1:PF.gen$gen.list$len)

	if(is.character(gene.set))
	{
		gene.idx <- match(gene.set, PF.gen$gen.list$names);
		if( sum(!is.na(gene.idx))==0)
			stop("All genes in the parameter 'gene.set' are available in current gene defintion file.\n");
			
		gene.set <- gene.idx[!is.na(gene.idx)];
	}

	if( length(which( gene.set > PF.gen$gen.list$len))>0 )
	{
		warning("The gene range is out of data set.");
		gene.set <- gene.set[- which( gene.set > PF.gen$gen.list$len ) ];
	}

	cat("* GENE.RANGE =", min(gene.set),"-", max(gene.set), "[", length(gene.set), "]\n");

	cpu.fun <- function( sect )
	{
		g.range0 <- gene.set[((sect-1)*n.perjob+1):(sect*n.perjob)];
		if ( length(gene.set) < sect*n.perjob )
			g.range0 <- gene.set[((sect-1)*n.perjob+1):length(gene.set)];

		PF.gen.new <- clone_plink_refer(PF.gen)

		res.cluster <- longskat_gene_task( r.model,
						PF.gen.new,
						PF.phe,
						g.range0,
						options$weights.common,
						options$weights.rare,
						options$run.cpp,
						options$snp.impute,
						options$rare.cutoff,
						options$test.type,
						options$verbose );

		return(res.cluster);
	}

	tm.start <- proc.time();
	n.perjob <- 100;
	n.subjob <- ceiling( length(gene.set)/100 );
	res.cluster <- list();

	if( options$n.cpu>1 )
	{
		cat("Starting parallel computing, snowfall/snow......\n");
		sfInit(parallel=TRUE, cpus=options$n.cpu, type="SOCK")

		sfExport("n.perjob", "gene.set", "r.model", "PF.gen", "PF.phe", "options", "PF.par" );
		res.cluster <- sfClusterApplyLB( 1:n.subjob, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n");

	}
	else
	{
		cat("Starting the LSKAT estimate for each gene......\n");
		for(k in 1:n.subjob)
			res.cluster[[k]] <- cpu.fun(k);
	}

	lskat.ret <- do.call("rbind", res.cluster);

	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");

	ret = list( result=lskat.ret, par=PF.par, mle=r.model);
	class(ret) <- "LSKAT.gen.plink";

	return(ret);
}

#public:
summary.LSKAT.gen.plink<-function(r.lskat)
{
	df.gene <- cbind( r.lskat$result[, c(2,3,4,5,6,8 ),drop=F], bonferroni=p.adjust( r.lskat$result[, "p.lskat" ], method="bonferroni"));
	df.gene <- df.gene[ order(df.gene[,6]), ];
	colnames(df.gene) <- c("Gene", "Chr.", "Pos.", "Total", "Rare", "p-value", "Bonferroni");

	return(df.gene[c(1:ifelse(NROW(df.gene)>20, 20, NROW(df.gene))),]);
}

#public:
print.LSKAT.gen.plink<-function(r.lskat, useS4=FALSE)
{
	cat("  LSKAT/Gene Summary\n");
	cat("[1] Paramaters:\n");
	cat("* PHE.LOG.FILE = ",   r.lskat$par$file.phe.long, "\n");
	cat("* PHE.LOG.TIME = ",   r.lskat$par$file.phe.time, "\n");
	cat("* PHE.COV.FILE = ",   r.lskat$par$file.phe.cov, "\n");
	cat("* Covariate Count = ",r.lskat$par$cov.count, "\n");
	cat("* Time Effect = ",    r.lskat$par$time.cov, "\n");
	cat("* Weight Rare = ",    r.lskat$par$weights.rare, "\n");
	cat("* Weight Common = ",  r.lskat$par$weights.common, "\n");

	cat("[2] MLE Results:\n");
	cat("* SIGMA_A =", 		   r.lskat$mle$par$sig.a, "\n");
	cat("* SIGMA_B =",         r.lskat$mle$par$sig.b, "\n");
	cat("* SIGMA_E =",         r.lskat$mle$par$sig.e, "\n");
	cat("* RHO =",             r.lskat$mle$par$rho, "\n");
	cat("* MU =",              r.lskat$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        r.lskat$mle$par$cov.effect, "\n");
	cat("* Beta(Time)=",       r.lskat$mle$par$time.effect, "\n");
	cat("* L(min) =",          r.lskat$mle$likelihood, "\n");

	#r.lskat$result
	# [1] index
	# [2] gene.name
	# [3] chr.
	# [4] min.pos
	# [5] SNP count
	# [6] Rare Counts
	# [7] Q of LSKAT
	# [8] p-value of LSKAT
	# [9] Q of L-Burden
	# [10] p-value of L-Burden
	df.gene <- cbind( r.lskat$result[, c(2,3,4,5,6,8),drop=F], bonferroni=p.adjust( r.lskat$result[, 8], method="bonferroni"));
	df.gene <- df.gene[ order(df.gene[,6]), ];
	colnames(df.gene) <- c("Gene", "Chr.", "Pos.", "Total", "Rare", "p-value", "Bonferroni");

	n.sig <- c();
	n.sig.i <- c();
	for( i in 20:4 )
	{
		n <- length( which( df.gene$bonferroni <10^(-i)) );
		if(n>0)
		{
			n.sig   <- c(n.sig, n);
			n.sig.i <- c(n.sig.i, i);
		}
	}

	cat("[3] Longitudinal SKAT Results:\n");
	cat("* Gene Count = ", NROW(r.lskat$result), "\n");
	cat("* Total SNP = ", sum(r.lskat$result[,5]), "\n");
	cat("* Rare SNP = ",  sum(r.lskat$result[,6]), "\n");

	if (length(n.sig.i)>0)
	{
		cat("* Significant Genes (Bonferroni Correct):\n")
		for(i in 1:length(n.sig.i))
			cat("* <= 10^(-", n.sig.i[i], ") Levels:", n.sig[i], "\n")
	}

	cat("* Top 20 Genes (Bonferroni Correct):\n");
	show(df.gene[c(1:ifelse(NROW(df.gene)>20, 20, NROW(df.gene))),]);
}

#public:
plot.LSKAT.gen.plink<-function( r.lskat, pdf.file=NA, title="", y.max=NA, bonferroni=F )
{
	if(is.na(pdf.file))
		pdf.file<-paste(r.lskat$par$file.phe.long, ".pdf", sep="");

	#CHR, POS, NAME, PV
	df.gene <- r.lskat$result[,c(3,4,2,8)];
	df.gene <- df.gene[ order(df.gene[,1], df.gene[,2]), ];

	if(bonferroni) df.gene[,4] <- df.gene[,4] * NROW(df.gene);

	pdf( pdf.file, width=6, height=3.5)

	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)

	draw_manhattan( df.gene[,c(1,4)],
					map.title = title,
					ifelse(bonferroni, 0.0001, 0.05)/NROW(df.gene),
					0.7,
					y.max = y.max );

	dev.off();

	cat("* Manhattan plot is save into ", pdf.file, "\n");
}


get_default_options<-function()
{
 	options <- list( rare.cutoff   = NULL,
					time.cov    = 0,
					g.maxiter     = 20,
					weights.common= c(0.5,0.5),
					weights.rare  = c(1,25),
					run.cpp       = F,
					verbose       = F,
					n.cpu         = 1,
					snp.impute    = "mean",
					intercept     = F,
					plink.path    = NULL,
					test.type     = "Joint",
					est.method    = "REML");

	return(options);
}

show_options<-function(options)
{
	cat( "* Covariate Time Effect: ",  options$time.cov, "\n");
	cat( "* Parallel Computing: ", ifelse( options$n.cpu>1, "Yes,", "No,"), options$n.cpu,  "CPU(s)\n");
	cat( "* Debug Output: ", ifelse( options$verbose, "Yes", "No"),"\n");
	cat( "* SNP Impute: ",  options$snp.impute, "\n");
	#cat( "* C/C++ Module Used Output: ", ifelse( options$run.cpp, "Yes", "No"), "\n");
	cat( "* Beta Weights for Common SNPs: ",  options$weights.common[1], options$weights.common[2], "\n");
	cat( "* Beta Weights for Rare SNPs: ",  options$weights.rare[1], options$weights.rare[2], "\n");
	cat( "* Common-Rare cutoff:", options$rare.cutoff, "\n");
	cat( "* Test type:", options$test.type, "\n");
	cat( "* Intercept:", options$intercept, "\n");
	cat( "* PLINK Path:", options$plink.path, "\n");
	cat( "* NULL Estimate Method:", options$est.method, "\n");
}

