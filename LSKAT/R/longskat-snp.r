#private:
est_snp_Q.R<-function(Y.delt, maf, Z, Y.t, X, par_null, time.effect )
{
	sig_a2  <- par_null[1]^2
	sig_b2  <- par_null[2]^2
	sig_e2  <- par_null[3]^2
	par_rho <- par_null[4]

	x.col <- NCOL(X);
	n <- dim(Y.delt)[1];
	m <- dim(Y.delt)[2];
	k <- 1;

	AR.1 <- array(0,dim=c(m,m));
	for(i in 1:m)
	for(j in 1:m)
		AR.1[i,j] <- par_rho^abs(i-j);

	V.j <- array(1, dim=c(m,m)) * sig_a2 +  AR.1 * sig_b2 + diag(1, m) * sig_e2
	V.j_1 <- solve(V.j);

	V.j_x <- list();
	Y.j_x <- list();
	M.j_x <- c();

	for(i in 1:n)
	{
		y.i <- Y.delt[i,];
		y.j_i <- V.j[!is.na(y.i), !is.na(y.i)];
		V.j_x[[i]] <- solve(y.j_i);
		Y.j_x[[i]] <- y.i[!is.na(y.i)];
		M.j_x[i] <- length(Y.j_x[[i]]);
	}

	Q.i <- 0;
	for(i in 1:n)
	#	Q.i <- Q.i + (Z[i,1]*t(rep(1,m))%*%V.j_x[[i]]%*%t(y.delt[i,,drop=F]))
		Q.i <- Q.i + Z[i,1]*t(rep(1, M.j_x[i] ))%*%V.j_x[[i]]%*% ( Y.j_x[[i]] )
		
		
	Q.v1 <- (Q.i)^2/2
	Q.v <- Q.v1[1,1]

	W0 <- array(0, dim=c(k, k));
	W1 <- array(0, dim=c(k, x.col));
	W2 <- array(0, dim=c(x.col, x.col));
	W3 <- array(0, dim=c(x.col, k));
	for(i in 1:n)
	{
		kr <- array(1, dim=c(1, M.j_x[i])) %*% V.j_x[[i]] %*% array(1, dim=c(M.j_x[i],1))
		kr.x <- X[i,];
		
		#need to discuss?
		#if( time.effect>0 )	
		#{
		#	t.i <- Y.t[i,,drop=F ];
		#	if ( length( which(is.na(t.i)) )> 0) 
		#		t.i <- t.i[ -which(is.na(t.i)) ];
		#	for(t in 1:time.effect)
		#		kr.X<- cbind(kr.X, t.i^t);
		#}
			
		W0 <- W0 + Z[i,1]^2 * kr;
		W1 <- W1 + Z[i,1] * kr %*% t( kr.x );
		W2 <- W2 + kr[1,1] *  kr.x %*%t ( kr.x );	
		W3 <- W3 + Z[i,1] * kr.x %*% kr;
	}
		
	Q.w1 <- (W0 - W1 %*% solve(W2) %*% W3)/2;
	Q.w <- Q.w1[1,1];
		
	#cat("Q.v==Q.v1", Q.v==Q.v1, "Q.w==Q.w1", Q.w==Q.w1, Q.v, Q.v1, Q.w, Q.w1, "\n");	
	
	r   <- sqrt(Q.v)/Q.w;
	sig_r2 <- 0;
	for(i in 1:n)
	{
		KY <- t(Y.j_x[[i]] - Z[i,1]*r)
		sig_ri <- KY%*%V.j_x[[i]]%*%t(KY);
		sig_r2 <- sig_r2 + sig_ri[1,1];
	}
	
	sig_r2 <- sig_r2/( n - 1- (x.col-1) -1 )
	p.v <-  pchisq(Q.v/(Q.w*sig_r2), df=1, lower.tail=F);
	
	return(list(v=Q.v, w=Q.w, r=r, sig_r2=sig_r2, chi2=Q.v/(Q.w*sig_r2), p.v=p.v ));
}

est_snp_Q<-function(Y.delt, maf, Z, Y.t, X, par_null, time.effect, run.cpp=T)
{
	r0 <- NA;
	t0 <- Sys.time()

	#if(!run.cpp)
		r0<- est_snp_Q.R(Y.delt, maf, Z, Y.t, X, par_null, time.effect)
	#else	
	#	r0<-.Call("est_snp_Q_C", Y.delt, Z,  X, maf, par_null );

	t1<-Sys.time()
	
	return(r0);		
}

get_weights<-function(maf, n, beta.common=c(0.5, 0.5), beta.rare=c(1, 25),rare.cutoff=NULL )
{
	if(is.null(rare.cutoff)) rare.cutoff <- 1/sqrt(2*n);
	
	wj <- c();	
	for( i in 1:length(maf))
	{
		if ( maf[i] > rare.cutoff )
			wj<- c( wj, dbeta( maf[i], beta.common[1], beta.common[2] ) )
		else
			wj<- c( wj, dbeta( maf[i], beta.rare[1], beta.rare[2] ) )
	}

	return(wj);
}

SKAT_Scale_Genotypes_snp= function( Z, weights.common=c(0.5,0.5), weights.rare=c(1,25),rare.cutoff=NULL )
{
	n<-dim(Z)[1];
	Z.maf <- colMeans(Z)/2;
	if(is.null(rare.cutoff)) rare.cutoff <- 1/sqrt(2*n);
	
	wr <- get_weights( Z.maf, n, 
					beta.common = weights.common, 
					beta.rare   = weights.rare, 
					rare.cutoff = rare.cutoff);
	Z <- Z * wr;
	
	return( list( new=Z, maf=Z.maf, rare = ifelse(Z.maf<rare.cutoff,1,0) ) )
}

#public:
longskat_snp_run<-function(r.model, snp, weights.common=c(0.5,0.5), weights.rare=c(1,25), run.cpp=F, rare.cutoff=NULL, debug=debug)
{
	get_snp_info<-function(snp)
	{
		s.miss <- which( is.na(snp) );
		s0 <- which( snp == 0 );
		s1 <- which( snp == 1 );
		s2 <- which( snp == 2 );

		if (length(s.miss)>0)
			snp <- snp[-s.miss];

		if ( mean(snp) > 1 ) snp <- 2 - snp;
		snp.imp <- as.matrix( snp, dim=c(length(snp),1) ); 
		snp.maf <- sum(snp.imp)/(length(snp)*2);

		return(list( snp = snp.imp, maf = snp.maf, nmiss = length(s.miss), miss = s.miss ) );
	}

    snp.info <- snp;
	if (class(snp)!="list") 
        snp.info <- get_snp_info(snp);
        
	if (length(snp.info$miss)>0)
	{
		if(debug) cat("! Missing SNP:", length(snp.info$miss), "\n");
		Y.delt <- Y.delt[ -snp.info$miss, ,drop=F];
		X <- X[-snp.info$miss, ,drop=F];
	}

	Y.delt <- r.model$y.delt; 
	Y.time <- r.model$y.time; 
	par_null <- c( r.model$par$sig_a, r.model$par$sig_b, r.model$par$sig_e, r.model$par$rho );
	X <- matrix( cbind(1, r.model$y.cov));

	Z.scale <- SKAT_Scale_Genotypes_snp( snp.info$snp, 
					weights.common = weights.common, 
					weights.rare   = weights.rare, 
					rare.cutoff    = rare.cutoff )

	Q <- est_snp_Q( Y.delt, 
					Z.scale$maf, 
					Z.scale$new, 
					Y.time, 
					X, 
					par_null, 
					time.effect, 
					run.cpp = run.cpp );

	P <- get_Qv_pvalue(Q$v, Q$w);

	r.lskat<- list( snp.total = length(Z.scale$maf), 
					snp.rare  = Z.scale$rare, 
					qv        = Q$v, 
					pv        = P$p.value, 
					maf       = snp.info$maf, 
					nmiss     = snp.info$nmiss,
	            	mle       = list(par=r.model$par, likelihood=r.model$likelihood) ); 
	class(r.lskat) <- "LSKAT.snp.ret";

	return( r.lskat )
}

#public:
print.LSKAT.snp.ret<-function(r.lskat, useS4 = FALSE)
{
	cat("  LSKAT/SNP Summary\n");	
	cat("[1] MLE Results:\n");	
	cat("* SIGMA_A =",         r.lskat$mle$par$sig_a, "\n");
	cat("* SIGMA_B =",         r.lskat$mle$par$sig_b, "\n");
	cat("* SIGMA_E =",         r.lskat$mle$par$sig_e, "\n");
	cat("* RHO =",             r.lskat$mle$par$rho, "\n");
	cat("* MU =",              r.lskat$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        r.lskat$mle$par$par_cov, "\n");
	cat("* Beta(Time)=",       r.lskat$mle$par$par_t, "\n");
	cat("* L(min) =",          r.lskat$mle$likelihood, "\n");

	cat("[2] Paramaters:\n");	
	cat("* Rare = ",           r.lskat$snp.rare, "\n");	
	cat("* MAF = ",            r.lskat$maf, "\n");	
	cat("* N.miss = ",         r.lskat$nmiss, "\n");	
	cat("* Q value = ",        r.lskat$qv, "\n");	
	cat("* p-value = ",		   r.lskat$pv, "\n");	
}

#private:
longskat_snp_task<-function(r.model, snp.range, file.gene.set, PF, weights.common, weights.rare, run.cpp=F, snp.impute="mean", rare.cutoff=NULL, debug=debug)
{
	gen.tb <- read.table(file.gene.set, sep=" ", header=F);

	rs.name <-c();
	gene.name <-c();
	rs.lst <- vector("list", length(which(!is.na(snp.range))) );
	rs.lst.idx <- 0;
	
	for(i in snp.range )
	{
		if (is.na(i)) next;
		
		snp.vec <- PF$snp.mat$genotypes[, idx, drop=F ];
		snp.info <- get_snp_plink_info (i, snp.vec, gen.tb );

		if (any(is.na(snp.info)) ) next;
		if (length(which(snp.info$maf>0.5))>0) browser();
		if (length(snp.info$maf)==0) next;

		ls <- longskat_snp_run( r.model, 
						snp.info, 
						weights.common, 
						weights.rare, 
						run.cpp, 
						rare.cutoff = rare.cutoff, 
						debug       = debug);
		
		rs.lst.idx <- rs.lst.idx + 1;
		rs.lst[[rs.lst.idx]] <- c(i, 
						snp.info$chr, 
						snp.info$loc, 
						snp.info$maf, 
						length(snp.info$miss), 
						ls$snp.total, 
						ls$snp.rare, 
						ls$qv, 
						ls$pv );
						
		rs.name     <- c(rs.name, as.character(snp.info$name) );
		gene.name   <- c(gene.name, as.character(snp.info$gene) );

		if (debug) cat(" [", i, "]", 
						snp.info$chr, 
						snp.info$loc, 
						as.character(snp.info$name), 
						as.character(snp.info$gene), 
						snp.info$maf, 
						length(snp.info$miss), 
						ls$snp.total, 
						ls$snp.rare, 
						ls$qv, 
						ls$pv, 
						"\n");
	}	

	rs <- do.call( "rbind", rs.lst );
	ret <- data.frame( id = rs[,1], 
					chr   = rs[,2], 
					loc   = rs[,3], 
					name  = rs.name, 
					gene  = gene.name, 
					maf   = rs[,4], 
					miss  = rs[,5], 
					total = rs[,6], 
					rare  = rs[,7], 
					Q     = rs[,8], 
					pv    = rs[,9]);

	return(ret);
}

#public
longskat_snp_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, 
						file.gene.set, 
						file.phe.long, 
						file.phe.cov, 
						file.phe.time=NULL, 
						snp.range=NULL,  
						options = list() )
{	
	cat( "[ LONGSKAT_SNP_PLINK ] Procedure\n");
	cat( "Checking the optional items......\n");

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
		options0[names(options)] <- options;
		options <- options0;
	}
	
	show_options(options);
	
	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	chk.phe <- check_pheno_file( file.phe.long, file.phe.time, file.plink.fam ) 
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, file.plink.fam, options$y.cov.count )
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	cat( "Starting to load all data files......\n");

	PF <- read_gen_phe_cov( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time, file.phe.cov ); 
	if(missing(options$y.cov.count) || is.na(options$y.cov.count) ) options$y.cov.count <- NCOL(PF$phe.cov);

	PF.par <- list(file.plink.bed = file.plink.bed, 
					file.plink.bim = file.plink.bim, 
					file.plink.fam = file.plink.fam, 
					file.phe.long  = file.phe.long, 
					file.phe.cov   = file.phe.cov,
					file.gene.set  = file.gene.set,
					y.cov.count    = options$y.cov.count, 
					y.cov.time     = options$y.cov.time,
					weights.common = options$weights.common, 
					weights.rare   = options$weights.rare,
					rare.cutoff    = options$rare.cutoff);

	if( is.na(options$y.cov.count) )
		options$y.cov.count <- NCOL(PF$phe.cov)-1;

	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	r.model <- longskat_est_model(PF$phe.long, 
					PF$phe.cov, 
					PF$phe.time, 
					y.cov.time = options$y.cov.time, 
					g.maxiter  = options$g.maxiter, 
					debug      = options$debug);

	if( class(r.model) != "LSKAT.null.model" ) 
		stop("! Failed to estimate the parameters of Covariance Compoment.");

	cat("* SIGMA_A =",    r.model$par$sig_a, "\n");
	cat("* SIGMA_B =",    r.model$par$sig_b, "\n");
	cat("* SIGMA_E =",    r.model$par$sig_e, "\n");
	cat("* RHO =",        r.model$par$rho, "\n");
	cat("* MU =",         r.model$par$mu, "\n");
	cat("* Beta(Cov) =",  r.model$par$par_cov, "\n");
	cat("* Beta(Time) =", r.model$par$par_t, "\n");
	cat("* L(min) =",     r.model$likelihood, "\n");

	snp.len <- dim(PF$snp.mat$genotypes)[2];
	if (is.null(snp.range)) snp.range<- c(1:snp.len)
	if( length(which( snp.range > snp.len))>0 )
	{
		warning("The snp range is out of data set.");
		snp.range <- snp.range[- which( snp.range > snp.len ) ];
	}

	cat("* SNP.RANGE =", min(snp.range),"-", max(snp.range), "[", length(snp.range), "]\n");
	
	cpu.fun<-function( sect )
	{
		snp.range0 <- snp.range[((sect-1)*n.percpu+1):(sect*n.percpu)];
		
		PF <- read_gen_phe_cov ( PF.par$file.plink.bed, PF.par$file.plink.bim, PF.par$file.plink.fam, 
						PF.par$file.phe.long, 
						PF.par$file.phe.cov );
		
		ret.cluster <- longskat_snp_task( r.model, 
						snp.range0, 
						PF.par$file.gene.set, 
						PF, 
						weights.common = options$weights.common, 
						weights.rare   = options$weights.rare, 
						run.cpp        = options$run.cpp, 
						snp.impute     = options$snp.impute, 
						rare.cutoff    = options$rare.cutoff, 
						debug          = options$debug );

		return(ret.cluster);
	}

	lskat.ret<-c();
	tm.start <- proc.time();
	if( options$n.cpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel=TRUE, cpus=options$n.cpu, type="SOCK")

		n.percpu <- ceiling( length(snp.range)/options$n.cpu );
		sfExport("n.percpu", "snp.range", "Y.delt", "Y.t", "X", "par_null", "options", "PF.par" );
		
		ret.cluster <- sfClusterApplyLB( 1:options$n.cpu, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n"); 
		lskat.ret <- do.call("rbind", ret.cluster);
	}
	else
	{
		cat("Starting the LSKAT estimate for each gene......\n");
		lskat.ret <- longskat_snp_task( r.model, 
						snp.range, 
						file.gene.set, 
						PF, 
						weights.common = options$weights.common, 
						weights.rare   = options$weights.rare, 
						run.cpp        = options$run.cpp, 
						snp.impute     = options$snp.impute, 
						rare.cutoff    = options$rare.cutoff, 
						debug          = options$debug );
	}
	
	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");
	
	ret = list( snp=lskat.ret, par=PF.par, mle=r.model);
	class(ret) <- "LSKAT.snp.plink";

	return(ret);
}

#public:
print.LSKAT.snp.plink<-function(r.lskat, useS4 = FALSE)
{
	summary.LSKAT.snp.plink(r.lskat);
}

#public:
summary.LSKAT.snp.plink<-function(r.lskat)
{
	cat("  LSKAT/SNP Summary\n");	
	cat("[1] Paramaters:\n");	
	cat("* PHE.LOG.FILE = ",   r.lskat$par$file.phe.long, "\n");	
	cat("* PHE.LOG.TIME = ",   r.lskat$par$file.phe.time, "\n");	
	cat("* PHE.COV.FILE = ",   r.lskat$par$file.phe.cov, "\n");	
	cat("* Covariate Count = ",r.lskat$par$y.cov.count, "\n");	
	cat("* Time Effect = ",    r.lskat$par$y.cov.time, "\n");	
	cat("* Weight Rare = ",    r.lskat$par$weights.rare, "\n");	
	cat("* Weight Common = ",  r.lskat$par$weights.common, "\n");	

	cat("[2] MLE Results:\n");	
	cat("* SIGMA_A =",         r.lskat$mle$par$sig_a, "\n");
	cat("* SIGMA_B =",         r.lskat$mle$par$sig_b, "\n");
	cat("* SIGMA_E =",         r.lskat$mle$par$sig_e, "\n");
	cat("* RHO =",             r.lskat$mle$par$rho, "\n");
	cat("* MU =",              r.lskat$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        r.lskat$mle$par$par_cov, "\n");
	cat("* Beta(Time)=",       r.lskat$mle$par$par_t, "\n");
	cat("* L(min) =",          r.lskat$mle$likelihood, "\n");

	cat("[3] Longitudinal SKAT Results:\n");
	cat("* SNP Count = ", NROW(r.lskat$snp), "\n");	

	# 1st = Chr, 2nd= Loc, 3rd=SNP.name, 4th=Gene.name, 5th=MAF, 
	# 6th = NMISS, 7th=Total, 8th=Rare, 9th=Q, 10th= PV,
	
	df.snp <- cbind( r.lskat$snp[, c(1,2,3,4,5,10),drop=F], bonferroni=p.adjust( r.lskat$snp[,10], method="bonferroni"));
	df.snp <- df.snp[ order(df.snp[,6]),,drop=F ]; 
	colnames(df.snp) <- c("SNP", "Gene", "Chr.", "Pos.", "MAF", "p-value", "Bonferroni");
	
	n.sig <- c();
	n.sig.i <- c();
	for( i in 20:4 )
	{
		n <- length( which( df.snp$bonferroni < 10^(-i) ) );
		if(n>0)
		{
			n.sig   <- c(n.sig, n);
			n.sig.i <- c(n.sig.i, i);
		}	
	}

	cat("  Summary of LSKAT result:\n");	
	cat("* PHE.LOG.FILE = ", r.lskat$par$file.phe.long, "\n")	
	cat("* PHE.COV.FILE = ", r.lskat$par$file.phe.cov, "\n")	
	if (length(n.sig.i)>0)
	{
		cat("* Significant Genes (Bonferroni Correct):\n")	
		for(i in 1:length(n.sig.i))
			cat("* <= 10^(-", n.sig.i[i], ") Levels:", n.sig[i], "\n")	
	}
	
	cat("* Top 20 SNPs (Bonferroni Correct):\n");
	show(df.snp[c(1:20),]);

}

#public:
plot.LSKAT.snp.plink<-function( r.lskat, pdf.file=NA, title="",  y.max=NA )
{
	if(is.na(pdf.file)) pdf.file<-paste(r.lskat$par$file.phe.long, ".pdf", sep="");

	#CHR, POS, NAME, PV
	df.snp <- r.lskat$snp[,c(1,2,3,10)];
	df.snp <- df.snp[ order(df.snp[,1], df.snp[,2]), ]; 
	df.snp[,4] <- df.snp[,4]*NROW(df.snp);

	pdf( pdf.file, width=6, height=3.5)

	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	draw_manhattan( df.snp[,c(1,3,4)], map.title=title, 0.0001/NROW(df.snp), 0.7, y.max= y.max );

	dev.off();

	cat("* Manhattan plot is save into ", pdf.file, "\n");	
}
  