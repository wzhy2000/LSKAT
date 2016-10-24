#private:
est_snp_Q.R<-function(Y.delt, maf, Z, Y.t, X, par_null, time.effect )
{
	sig.a2  <- par_null[1]^2
	sig.b2  <- par_null[2]^2
	sig.e2  <- par_null[3]^2
	par_rho <- par_null[4]

	x.col <- NCOL(X);
	n <- dim(Y.delt)[1];
	m <- dim(Y.delt)[2];
	k <- 1;

	AR.1 <- array(0,dim=c(m,m));
	for(i in 1:m)
	for(j in 1:m)
		AR.1[i,j] <- par_rho^abs(i-j);

	V.j <- array(1, dim=c(m,m)) * sig.a2 +  AR.1 * sig.b2 + diag(1, m) * sig.e2
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
	sig.r2 <- 0;
	for(i in 1:n)
	{
		KY <- t(Y.j_x[[i]] - Z[i,1]*r)
		sig.ri <- KY%*%V.j_x[[i]]%*%t(KY);
		sig.r2 <- sig.r2 + sig.ri[1,1];
	}
	
	sig.r2 <- sig.r2/( n - 1- (x.col-1) -1 )
	p.v <-  pchisq(Q.v/(Q.w*sig.r2), df=1, lower.tail=F);
	
	return(list(v=Q.v, w=Q.w, r=r, sig.r2=sig.r2, chi2=Q.v/(Q.w*sig.r2), p.v=p.v ));
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
longskat_snp_test<-function(r.model, snp, 
				weights.common=c(0.5,0.5), 
				weights.rare=c(1,25), 
				snp.impute = "mean",
				rare.cutoff = NULL, 
				verbose = FALSE)
{
	run.cpp <- F;
	
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

    snp.obj <- snp;
	if (class(snp)!="list") 
        snp.obj <- get_snp_info(snp);
        
	Y.delt <- r.model$phe.delt; 
	Y.time <- r.model$phe.time; 
	X <- r.model$phe.cov; 

	if (!is.null(snp.obj$miss) && length(snp.obj$miss)>0)
	{
		if(verbose) cat("! Missing SNP:", length(snp.obj$miss), "\n");
		Y.delt <- Y.delt[ -snp.obj$miss, ,drop=F];
		X <- X[-snp.obj$miss, ,drop=F];
	}

	if(is.null(rare.cutoff))
		rare.cutoff <- 1/sqrt(2*NROW(X));

	par_null <- c( r.model$par$sig.a, r.model$par$sig.b, r.model$par$sig.e, r.model$par$rho );
	X <- matrix( cbind(1, r.model$phe.cov));

	Z.scale <- SKAT_Scale_Genotypes_snp( snp.obj$snp, 
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

	r.lskat<- list( qv        = Q$v, 
					pv        = P$p.value, 
					maf       = snp.obj$maf, 
					nmiss     = snp.obj$nmiss,
					rare      = ifelse(snp.obj$maf<=rare.cutoff, TRUE, FALSE),
	            	mle       = list(par=r.model$par, likelihood=r.model$likelihood) ); 
	class(r.lskat) <- "LSKAT.snp.ret";

	return( r.lskat )
}

#public:
print.LSKAT.snp.ret<-function(r.lskat, useS4 = FALSE)
{
	cat("  LSKAT/SNP Summary\n");	
	cat("[1] MLE Results:\n");	
	cat("* SIGMA_A =",         r.lskat$mle$par$sig.a, "\n");
	cat("* SIGMA_B =",         r.lskat$mle$par$sig.b, "\n");
	cat("* SIGMA_E =",         r.lskat$mle$par$sig.e, "\n");
	cat("* RHO =",             r.lskat$mle$par$rho, "\n");
	cat("* MU =",              r.lskat$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        r.lskat$mle$par$cov.effect, "\n");
	cat("* Beta(Time)=",       r.lskat$mle$par$time.effect, "\n");
	cat("* L(min) =",          r.lskat$mle$likelihood, "\n");

	cat("[2] Paramaters:\n");	
	cat("* Rare = ",           r.lskat$snp.rare, "\n");	
	cat("* MAF = ",            r.lskat$maf, "\n");	
	cat("* N.miss = ",         r.lskat$nmiss, "\n");	
	cat("* Q value = ",        r.lskat$qv, "\n");	
	cat("* p-value = ",		   r.lskat$pv, "\n");	
}

#private:
longskat_snp_task<-function(r.model, PF.gen, PF.phe, 
			snp.set, 
			weights.common, 
			weights.rare, 
			run.cpp = FALSE, 
			snp.impute = "mean", 
			rare.cutoff = NULL, 
			verbose = FALSE)
{
	rs.name <-c();
	gene.name <-c();
	rs.lst <- vector("list", length(which(!is.na(snp.set))) );
	rs.lst.idx <- 0;
	
	for(i in snp.set )
	{
		if (is.na(i)) next;
		
		snp.obj <- get_snp_mat (PF.gen, i, snp.impute);
		if(is.null(snp.obj)) 
			next;

		ls <- longskat_snp_test( r.model, 
						snp.obj, 
						weights.common, 
						weights.rare, 
						#run.cpp, 
						snp.impute  = snp.impute,
						rare.cutoff = rare.cutoff, 
						verbose     = verbose);
		
		rs.lst.idx <- rs.lst.idx + 1;
		rs.lst[[rs.lst.idx]] <- c(i, 
						snp.obj$chr, 
						snp.obj$loc, 
						snp.obj$maf, 
						ls$nmiss, 
						ls$rare, 
						ls$qv, 
						ls$pv );
						
		rs.name     <- c(rs.name, as.character(snp.obj$name) );
		gene.name   <- c(gene.name, as.character(snp.obj$gene) );

		if (verbose) cat(" [", i, "]", 
						as.character(snp.obj$name), 
						as.character(snp.obj$gene), 
						snp.obj$chr, 
						snp.obj$loc, 
						snp.obj$maf, 
						ls$nmiss, 
						ls$rare, 
						ls$qv, 
						ls$pv, "\n");
	}	


	rs <- do.call( "rbind", rs.lst );
	ret <- data.frame( 
					index     = unlist(rs[,1,drop=F]), 
					snp.name  = unlist(rs.name), 
					gene.name = unlist(gene.name), 
					chr       = unlist(rs[,2,drop=F]), 
					pos       = unlist(rs[,3,drop=F]), 
					maf       = unlist(rs[,4,drop=F]), 
					nmiss     = unlist(rs[,5,drop=F]), 
					rare      = unlist(rs[,6,drop=F]), 
					q.lskat   = unlist(rs[,7,drop=F]), 
					p.lskat   = unlist(rs[,8,drop=F]));

	return(ret);
}

#public
longskat_snp_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, 
						file.phe.long, 
						file.phe.cov, 
						file.phe.time = NULL, 
						file.gene.set = NULL,
						snp.set = NULL,  
						options = list(),
						verbose = FALSE)
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
	
	chk.phe <- check_pheno_file( file.phe.long, file.phe.time, file.plink.fam ) 
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, file.plink.fam)
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	cat( "Starting to load all data files......\n");

	PF.gen <- read_gen_plink ( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, options$plink.path )

	PF.phe <- read_phe_cov ( file.phe.long, file.phe.time, 	file.phe.cov, PF.gen );

	PF.par <- list( file.plink.bed = file.plink.bed,
					file.plink.bim = file.plink.bim,
					file.plink.fam = file.plink.fam,
					file.phe.long  = file.phe.long,
					file.phe.time  = file.phe.time,
					file.phe.cov   = file.phe.cov,
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

	if (is.null(snp.set))
		snp.set<- c(1:NROW(PF.gen$snp$bim))

	if(is.character(snp.set))
	{
		snp.idx <- match(snp.set, PF.gen$snp$bim[,2]);
		if( sum(!is.na(snp.idx))==0)
			stop("All SNPs in the parameter 'snp.set' are available in current gene defintion file.\n");
			
		snp.set <- snp.idx[!is.na(snp.idx)];
	}

	if( length(which( snp.set > NROW(PF.gen$snp$bim[,2]) ))>0 )
	{
		warning("The snp range is out of data set.");
		snp.set <- snp.set[- which( snp.set > NROW(PF.gen$snp$bim[,2]) ) ];
	}

	cat("* SNP.RANGE =", min(snp.set),"-", max(snp.set), "[", length(snp.set), "]\n");
	
	cpu.fun<-function( sect )
	{
		g.range0 <- snp.set[((sect-1)*n.perjob+1):(sect*n.perjob)];
		if ( length(snp.set) < sect*n.perjob )
			g.range0 <- snp.set[((sect-1)*n.perjob+1):length(snp.set)];

		PF.gen.new <- clone_plink_refer(PF.gen)

		res.cluster <- longskat_snp_task( r.model,
						PF.gen.new,
						PF.phe,
						g.range0,
						options$weights.common,
						options$weights.rare,
						options$run.cpp,
						options$snp.impute,
						options$rare.cutoff,
						options$verbose );


		return(res.cluster);
	}

	tm.start <- proc.time();
	n.perjob <- 100;
	n.subjob <- ceiling( length(snp.set)/100 );
	res.cluster <- list();

	if( options$n.cpu>1 )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel=TRUE, cpus=options$n.cpu, type="SOCK")

		sfExport("n.perjob", "snp.set", "r.model", "PF.gen", "PF.phe", "options", "PF.par");
		ret.cluster <- sfClusterApplyLB( 1:n.subjob, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n"); 
		lskat.ret <- do.call("rbind", ret.cluster);
	}
	else
	{
		cat("Starting the LSKAT estimate for each SNP......\n");
		for(k in 1:n.subjob)
			res.cluster[[k]] <- cpu.fun(k);
	}
	
	lskat.ret <- do.call("rbind", res.cluster);

	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");


	ret = list( result=lskat.ret, par=PF.par, mle=r.model);
	class(ret) <- "LSKAT.snp.plink";

	return(ret);
}

#public:
summary.LSKAT.snp.plink<-function(r.lskat)
{
	return(r.lskat$result);
}

#public:
print.LSKAT.snp.plink<-function(object, useS4 = FALSE)
{
	cat("  LSKAT/SNP Summary\n");	
	cat("[1] Paramaters:\n");	
	cat("* PHE.LOG.FILE = ",   object$par$file.phe.long, "\n");	
	cat("* PHE.LOG.TIME = ",   object$par$file.phe.time, "\n");	
	cat("* PHE.COV.FILE = ",   object$par$file.phe.cov, "\n");	
	cat("* Covariate Count = ",object$par$cov.count, "\n");	
	cat("* Time Effect = ",    object$par$time.cov, "\n");	
	cat("* Weight Rare = ",    object$par$weights.rare, "\n");	
	cat("* Weight Common = ",  object$par$weights.common, "\n");	

	cat("[2] MLE Results:\n");	
	cat("* SIGMA_A =",         object$mle$par$sig.a, "\n");
	cat("* SIGMA_B =",         object$mle$par$sig.b, "\n");
	cat("* SIGMA_E =",         object$mle$par$sig.e, "\n");
	cat("* RHO =",             object$mle$par$rho, "\n");
	cat("* MU =",              object$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        object$mle$par$cov.effect, "\n");
	cat("* Beta(Time)=",       object$mle$par$time.effect, "\n");
	cat("* L(min) =",          object$mle$likelihood, "\n");

	cat("[3] Longitudinal SKAT Results:\n");
	cat("* SNP Count = ", NROW(object$result), "\n");	

	# 1st=index, 2nd=SNP.name, 3rd=Gene.name, 4th = Chr, 5th= Loc, 
	# 6th=MAF, 7th = NMISS, 8th=Rare, 9th=q.lskat, 10th= p.lskat,
	
	df.snp <- cbind( object$result[, c(2,3,4,5,6,7,10),drop=F], bonferroni=p.adjust( object$result[,10], method="bonferroni"));
	df.snp <- df.snp[ order(df.snp[,6]),,drop=F ]; 
	colnames(df.snp) <- c("SNP", "Gene", "Chr.", "Pos.", "MAF", "NMISS", "p-value", "Bonferroni");
	
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
	cat("* PHE.LOG.FILE = ", object$par$file.phe.long, "\n")	
	cat("* PHE.COV.FILE = ", object$par$file.phe.cov, "\n")	
	if (length(n.sig.i)>0)
	{
		cat("* Significant Genes (Bonferroni Correct):\n")	
		for(i in 1:length(n.sig.i))
			cat("* <= 10^(-", n.sig.i[i], ") Levels:", n.sig[i], "\n")	
	}
	
	cat("* Top 20 SNPs (Bonferroni Correct):\n");
	show(df.snp[c(1:ifelse(NROW(df.snp)>20, 20,NROW(df.snp))),]);

}

#public:
plot.LSKAT.snp.plink<-function( object, pdf.file=NA, title="",  y.max=NA )
{
	if(is.na(pdf.file)) pdf.file<-paste(object$par$file.phe.long, ".pdf", sep="");

	#CHR, POS, NAME, PV
	df.snp <- object$result[,c(4,5,2,10)];
	df.snp <- df.snp[ order(df.snp[,1], df.snp[,2]), ]; 
	df.snp[,4] <- df.snp[,4]*NROW(df.snp);

	pdf( pdf.file, width=6, height=3.5)

	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	draw_manhattan( df.snp[,c(1,4)], map.title=title, 0.0001/NROW(df.snp), 0.7, y.max= y.max );

	dev.off();

	cat("* Manhattan plot is save into ", pdf.file, "\n");	
}
  