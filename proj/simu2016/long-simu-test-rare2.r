
source("long-simu-run.r");

simu.snp.mat<- function(snp.file1, snp.file2, snp.range, n.sample, mat.count, rare.cutoff, rare.count.range, common.count.range, ncores)
{
	min.sect <- snp.range[1];
	max.sect <- snp.range[2];

	snp.hap <- read.table(snp.file1, header=F);
	snp.pos <- read.table(snp.file2, header=T);
	snp.maxpos <- max(snp.pos$CHROM_POS);
	snp.minpos <- min(snp.pos$CHROM_POS);

	if (is.null(rare.cutoff))
		rare.cutoff <- 1/sqrt(2*n.sample);

	mat.ret <- mclapply(1:mat.count, function(i)
	{
		set.seed( runif(1)*1000 + i*1000 + proc.time()[3] );
		snp.mat <- NULL;
		n.mat <- 1;

		while(n.mat <= 1)
		{
			snp.start <- as.integer(runif(1, min=snp.minpos, max=snp.maxpos));
			snp.sect  <- as.integer(runif(1, min=min.sect,   max=max.sect));

			p.sel <- which( snp.pos$CHROM_POS>snp.start & snp.pos$CHROM_POS<(snp.start+snp.sect) );

			snp.pair <- sample(1:(dim(snp.hap)[2]-2));
			snp.mat1 <- snp.hap[ snp.pair[1:n.sample], p.sel + 2 ];
			snp.mat2 <- snp.hap[ snp.pair[(n.sample+1):(2*n.sample)], p.sel + 2 ];
			snp.mat <- snp.mat1 + snp.mat2 -2;

			maf <- colMeans(snp.mat)/2;
			m.same <- which( maf==1 | maf==0 );
			if (length(m.same)>0)
				snp.mat <- snp.mat[, -m.same, drop=F ];

			if(NCOL(snp.mat) <= 1) next;

			#make sure of 2: Minor Allels, 0:Major Allels
			maf <- colMeans(snp.mat)/2;
			m.minor <- which( maf>0.5 );
			if (length(m.minor)>0)
				for(i in 1:length(m.minor))
					snp.mat[,m.minor[i]] <- 2 - snp.mat[,m.minor[i]];

			if(NCOL(snp.mat) <= 1) next;

			maf <- colMeans(snp.mat)/2;
			m.smallmaf <- which( maf < 5/n.sample);
			if (length(m.smallmaf)>0)
				snp.mat <- snp.mat[, -m.smallmaf,drop=F ];

			if(NCOL(snp.mat) <= 1) next;

			# at lease 1 RARE SNPs
			maf <- colMeans(snp.mat)/2;
			n.rare <- length(which(maf<rare.cutoff));

			#if (n.rare<=10)
			#	next;
			#
			#if ( n.rare>20 && n.rare>length(maf)*0.2 )
			#{
			#	n.rare0 <- as.integer(length(maf)*0.2)+1;
			#	if (n.rare0>1)
			#	{
			#		rare.rm <- which(maf<rare.cutoff)[c(n.rare0:n.rare)];
			#		if (length(rare.rm)>0) snp.mat <- snp.mat[, -rare.rm, drop=F ];
			#	}
			#}

			if ( n.rare < rare.count.range[1] || length(maf) - n.rare < common.count.range[1] )
				next;

			snp.keep <- c( which( maf < rare.cutoff )[c(1:rare.count.range[2])],
						   which( maf >= rare.cutoff )[c(1:common.count.range[2])] );
			snp.keep <- snp.keep[!is.na(snp.keep)]

			if (length(snp.keep)>0)
				snp.mat <- snp.mat[, snp.keep, drop=F ];

			n.mat <- n.mat + 1;
		}

		return(snp.mat);

	}, mc.cores=ncores);

	return(mat.ret);
}

simu.long.phe.power<-function( n.sample, snp.mat, par, f.simu )
{
	y.random <- f.simu(n.sample, par);

	par_cov <- cbind( rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	if(par$intercept)
		y.fixed  <- cbind(1, par_cov) %*% c(par$mu, par$a, par$b ) %*% t( rep(1, par$times))
	else
		y.fixed  <- par_cov %*% c(par$a, par$b ) %*% t( rep(1, par$times));

	y.time <- 0;
	if(!is.null(par$y.cov.time) && par$y.cov.time>0)
	{
		y.time <- par$par_t[1]*t(matrix(rep(c(1:par$times),n.sample), nrow=par$times))
		y.time <- y.time + par$par_t[2]*(t(matrix(rep(c(1:par$times),n.sample), nrow=par$times))**2)
	}

	rare.cutoff <- par$rare.cutoff;
	if (is.null(rare.cutoff)) rare.cutoff <- 1/sqrt(2*n.sample);

	maf <- colMeans(snp.mat)/2;
	n.snp <- length(maf);

	snp.comm <- 0;
	snp.rare <- 0;
	X.snp <- c();
	y.snp <- 0;
	for( i in 1:n.snp)
	{
		if (maf[i]<rare.cutoff)
		{
			#if (snp.rare < n.snp*0.2 )
			# the ratio of common causal snp is about 0.3 in 2015 test
			#if (snp.rare < 10 )
			if (snp.comm < par$max.rare.causal )
			{
				sign <- ifelse( runif(1) <= par$positive.effect, 1, -1 );
				#y.snp <- y.snp + sign * snp.mat[,i] * (rnorm(1, par$rare.c1, par$effect.sd) *abs(log10(maf[i])));
				y.snp <- y.snp + sign * snp.mat[,i] * par$rare.c1 *abs(log10(maf[i]));
				snp.rare <- snp.rare + 1;
				X.snp <-cbind(X.snp, snp.mat[,i]);
			}
		}
		else
		{
			# the ratio of common causal snp is about 0.3 in 2015 test
			#if (snp.comm < 4  )
			if ( snp.comm < par$max.common.causal )
			{
				sign <- ifelse( runif(1) <= par$positive.effect, 1, -1 );
				#y.snp <- y.snp + sign*snp.mat[,i] * rnorm(1, par$common.c1, par$effect.sd);
				y.snp <- y.snp + sign*snp.mat[,i] * par$common.c1;
				snp.comm <- snp.comm+1;
				X.snp <-cbind(X.snp, snp.mat[,i]);
			}
		}
	}

	y.snp <- y.snp%*%t(rep(1, par$times))

	y <- y.random + y.fixed + y.snp + y.time;

	h2 <- var(c(y.snp))/var(c(y.snp + y.random));
	get_stats<-function(y){return(c(mean(y, na.rm=T), min(y, na.rm=T), max(y, na.rm=T), sd(y, na.rm=T)))}
	ys <- rbind(get_stats(y), get_stats(y.fixed), get_stats(y.random), get_stats(y.snp), get_stats(y.time))
	rownames(ys) <- c("Y","Y.fixed","Y.random","Y.snp","Y.time");
	colnames(ys) <- c("Mean","Min","Max","SD");

show(round(ys,2));
cat( "h2=", round(h2,2),"Rare=", snp.rare, "\n");


	return(list(y=y, cov=par_cov, h2=h2, X.snp = X.snp))
}

power.test<-function(par, nloop, nsample, phe.dist, phe.cov, power.rdata, ncores=1 )
{
	cat(" [Power Test]\n")
	cat(" * loop=", nloop, "\n")
	cat(" * nsample=", nsample, "\n")
	cat(" * par=",  unlist(par), "\n")
	cat(" * phe.dist=", phe.dist, "\n")
	cat(" * phe.cov=", phe.cov, "\n")
	cat(" * ret.rdata=", power.rdata, "\n")

	r.1 <- r.2 <- r.3 <- NULL;
	f.simu <- f.mn.ar1;

	if ( phe.dist=="mn" && phe.cov=="ar1") f.simu <- f.mn.ar1;
	if ( phe.dist=="mn" && phe.cov=="sad") f.simu <- f.mn.sad;
	if ( phe.dist=="mn" && phe.cov=="cm" ) f.simu <- f.mn.cm;
	if ( phe.dist=="mt"  ) f.simu <- f.mt.ar1;
	if ( phe.dist=="msn" ) f.simu <- f.sn.ar1;
	if ( phe.dist=="mmn" ) f.simu <- f.mmn.ar1;

	# 50% positive effect
	par$positive.effect<- 0.5;
	r.3 <-check.power( par, nloop, nsample, f.simu, ncores )
	save(par, r.1, r.2, r.3, file=power.rdata);

	# 100% positive effect
	par$positive.effect<- 1;
	r.1 <- check.power( par, nloop, nsample, f.simu, ncores )
	save(par, r.1, r.2, r.3, file=power.rdata);

	# 80% positive effect
	par$positive.effect<- 0.8;
	r.2 <-check.power( par, nloop, nsample, f.simu, ncores )
	save(par, r.1, r.2, r.3, file=power.rdata);


	get_sig_level<-function( r.x, a.level )
	{
		return(c(
			length(which(r.x[,13]<a.level)),
			length(which(r.x[,14]<a.level)),
			length(which(r.x[,15]<a.level)),
			length(which(r.x[,18]<a.level)),
			length(which(r.x[,19]<a.level)),
			length(which(r.x[,20]<a.level)),
			length(which(r.x[,21]<a.level))));
	}
	a.level<-10^(-6);
	r.power<-c();
	r.power<-rbind( r.power, c(nsample, 1,  get_sig_level(r.1, a.level) ) );
	r.power<-rbind( r.power, c(nsample, 2,  get_sig_level(r.2, a.level) ) );
	r.power<-rbind( r.power, c(nsample, 3,  get_sig_level(r.3, a.level) ) );
	colnames(r.power) <- c("sample", "beta.plus", "lskat.pv", "lskat.burden.pv", "lskat.snp.pv", "skat.bl.pv", "burden.bl.pv", "skat.mu.pv", "burden.mu.pv");
	r.power6 <- r.power;

	a.level<-10^(-4);
	r.power<-c();
	r.power<-rbind( r.power, c(nsample, 1,  get_sig_level(r.1, a.level) ) );
	r.power<-rbind( r.power, c(nsample, 2,  get_sig_level(r.2, a.level) ) );
	r.power<-rbind( r.power, c(nsample, 3,  get_sig_level(r.3, a.level) ) );
	colnames(r.power) <- c("sample", "beta.plus", "lskat.pv", "lskat.burden.pv", "lskat.snp.pv", "skat.bl.pv", "burden.bl.pv", "skat.mu.pv", "burden.mu.pv");
	r.power4 <- r.power;

	a.level<-10^(-2);
	r.power<-c();
	r.power<-rbind( r.power, c(nsample, 1,  get_sig_level(r.1, a.level) ) );
	r.power<-rbind( r.power, c(nsample, 2,  get_sig_level(r.2, a.level) ) );
	r.power<-rbind( r.power, c(nsample, 3,  get_sig_level(r.3, a.level) ) );
	colnames(r.power) <- c("sample", "beta.plus", "lskat.pv", "lskat.burden.pv", "lskat.snp.pv", "skat.bl.pv", "burden.bl.pv", "skat.mu.pv", "burden.mu.pv");
	r.power2 <- r.power;

	save( r.power6, r.power4, r.power2, par, r.1, r.2, r.3, file=power.rdata );
}


testR <- function(power.test, par, ret.rdata, nloop=1000, nrep=1000, nsample=500, phe.dist="mn", phe.cov="ar1", ncores=1)
{
	par.default <- list(
		test.type   = "Joint",
	    mu          = 1,
		a           = 0.5,
		b           = 0.5,
		rho         = 0.75,
		sig_a       = 0.8,
		sig_b       = 0.8,
		sig_e       = 0.8,
		times       = 8,
		intercept   = F,
		y.cov.time  = 0,
		par_t       = c(0.2, -0.08),
	    snprange    = c(5*1000, 30*1000),
		a.level     = 10^(-6),
		w.common    = c(0.5, 0.5),
		w.rare      = c(1,25),
		common.c1   = 0.12,
		rare.c1     = 0.08,
		effect.sd   = 0.05,
		max.rare.causal = 10,
		max.common.causal = 4,
	    positive.effect = 1,
		rare.cutoff = 0.05 );

	par.default[names(par)] <- par;
	par <- par.default;

	#power.test <- T;
	#nloop      <- 2;
	#nrep       <- 1000;
	#phe.sample <- 500;
	#phe.dist   <- "msn"
	#phe.cov    <- "ar1"
	#phe.par    <- "0.7,0.5"
	#ret.rdata  <- "test-xxx.rdata"

	if ( power.test )
		power.test( par, nloop, nsample, phe.dist, phe.cov, ret.rdata, ncores )
	else
		type1.test( par, nloop, nrep, nsample, phe.dist, phe.cov, ret.rdata, ncores);
}
