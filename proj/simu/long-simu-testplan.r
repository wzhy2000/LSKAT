
source("long-simu-run.r");

simu.snp.mat<- function(snp.file1, snp.file2, snp.range, n.sample, mat.count, rare.cutoff,ncores)
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

			# at lease 50 RARE SNPs
			maf <- colMeans(snp.mat)/2;
			n.rare <- length(which(maf<rare.cutoff));
			if (n.rare<50)
				next;

			#if ( n.rare>=1 && n.rare>length(maf)*1/3 )
			#{
			#	n.rare0 <- as.integer(length(maf)*1/3)+1;
			#	if (n.rare0>1)
			#	{
			#		rare.rm <- which(maf<rare.cutoff)[c(n.rare0:n.rare)];
			#		if (length(rare.rm)>0) snp.mat <- snp.mat[, -rare.rm, drop=F ];
			#	}
			#}

			n.mat <- n.mat + 1;
		}
		
		return(snp.mat);
		
	}, mc.cores=ncores);

	return(mat.ret);
}

simu.long.phe.power<-function( n.sample, snp.mat, par, f.simu )
{
	par_cov <- cbind(1, rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	y <-  f.simu(n.sample, par) + par_cov%*%c(1, par$a, par$b )%*%t(rep(1, par$times));
	# 2 unobserved common SNP
	y <- y + round(runif(n.sample, -0.4, 2.4)) * 0.27;
	y <- y + round(runif(n.sample, -0.4, 2.4)) * 0.27;
	
	rare.cutoff <- par$rare.cutoff;
	if( is.null(rare.cutoff) ) rare.cutoff <- 1/sqrt(2*n.sample);
	
	maf    <- colMeans(snp.mat)/2;
	n.rare <- length(which(maf<rare.cutoff))
	n.snp  <- length(maf);
	
	rare.ratio <- cumsum( par$rare.ratio );
	snp.comm <- 0;
	snp.rare <- 0;
	for( i in 1:n.snp)
	{
		if (maf[i]<rare.cutoff)
		{
			effect <- par$rare.effect[ which.max(snp.rare/n.rare <= rare.ratio) ]
			direct <- par$rare.direct[ which.max(snp.rare/n.rare <= rare.ratio) ]
			if (is.na(effect)) effect <- par$rare.c1* direct * abs(log10(maf[i]));
			y <- y + snp.mat[,i] * effect;
			snp.rare <- snp.rare + 1;
		}
		else
		{
			if (snp.comm < par$common.cause )
			{
				y <- y + snp.mat[,i] * par$common.effect;
				snp.comm <- snp.comm+1;
			}
		}
	}
	
	return(list(y=y, cov=par_cov[,-1]))
}


testR <- function(power.test, ret.rdata, nloop=1000, nrep=1000, nsample=500, phe.dist="mn", phe.cov="ar1", 
				  a=NA, b=NA, rho=NA, sig_a=NA, sig_b=NA, sig_e=NA,
				  phe.par1=NULL, phe.par2=NULL, phe.par3=NULL, phe.par4=NULL, phe.par5=NULL, phe.par6=NULL, ncores=1)
{
	par <- list( 
		a           = 0.5, 
		b           = 0.5, 
		rho         = 0.75,
		sig_a       = 0.8, 
		sig_b       = 0.8, 
		sig_e       = 0.8, 
		beta.c1     = 0.4, 
		times       = 6, 
	    snprange    = c(30*1000, 30*1000),
		a.level     = 10^(-6), 
		w.common    = c(0.5, 0.5),
		w.rare      = c(1,25),
		common.effect = 0.27,
		common.cause = 3,
		rare.ratio  = c( 0.4, 0.6 ), 
	    rare.effect = c( NA, 0 ),
	    rare.direct = c( 1, 0 ),
	    rare.c1     = 0.01,
		rare.cutoff = 0.05 );

	if(!is.na(a))     par$a     <- a;
	if(!is.na(b))     par$b     <- b;
	if(!is.na(rho))   par$rho   <- rho;
	if(!is.na(sig_a)) par$sig_a <- sig_a;
	if(!is.na(sig_b)) par$sig_b <- sig_b;
	if(!is.na(sig_e)) par$sig_e <- sig_e;

	#power.test <- T;
	#nloop      <- 2;
	#nrep       <- 1000;
	#phe.sample <- 500;
	#phe.dist   <- "msn"
	#phe.cov    <- "ar1"
	#phe.par    <- "0.7,0.5"
	#ret.rdata  <- "test-xxx.rdata"

	par$par1 <- phe.par1;
	par$par2 <- phe.par2;
	par$par3 <- phe.par3;
	par$par4 <- phe.par4;
	par$par5 <- phe.par5;
	par$par6 <- phe.par6;

	if ( power.test )
		power.test( par, nloop, nsample, phe.dist, phe.cov, ret.rdata, ncores )
	else
		type1.test( par, nloop, nrep, nsample, phe.dist, phe.cov, ret.rdata, ncores);
}

#---------------------------------------------------
# Interface ()
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L1.mn.ar1.500.rdata sig.a=0.8 sig.e=0.8   < long-simu-testplan.r > power-L1-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.4 ret.rdata=power.L2.mn.ar1.500.rdata sig.a=0.4 sig.e=0.4   < long-simu-testplan.r > power-L2-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.2 ret.rdata=power.L3.mn.ar1.500.rdata sig.a=0.2 sig.e=0.2   < long-simu-testplan.r > power-L3-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L4.mn.ar1.500.rdata sig.a=0.2 sig.e=0.2   < long-simu-testplan.r > power-L4-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.2 ret.rdata=power.L5.mn.ar1.500.rdata sig.a=0.8 sig.e=0.8   < long-simu-testplan.r > power-L5-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mt phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L6.mt.ar1.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L6-mt-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mt phe.cov=ar1 phe.par=0.7,0.4 ret.rdata=power.L7.mt.ar1.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L7-mt-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=msn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L8.msn.ar1.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L8-msn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=msn phe.cov=ar1 phe.par=0.7,0.4 ret.rdata=power.L9.msn.ar1.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L9-msn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mmn phe.cov=ar1 phe.par=0.5,0.8,0.5,0.8,0.6 ret.rdata=power.L10.mmn.ar1.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L10-mmn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mmn phe.cov=ar1 phe.par=0.5,0.8,0.3,0.4,0.7 ret.rdata=power.L11.mmn.ar1.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L11-mmn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=sad phe.par=0.7,0.8 ret.rdata=power.L12.mn.sad.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L12-mn-sad-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=sad phe.par=0.7,0.4 ret.rdata=power.L13.mn.sad.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L13-mn-sad-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=cm phe.par=0.7,0.8 ret.rdata=power.L14.mn.cm.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L14-mn-cm-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=cm phe.par=0.7,0.4 ret.rdata=power.L15.mn.cm.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L15-mn-cm-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.87,1.48 ret.rdata=power.L16.mn.ar1.500.rdata sig.a=0.001 sig.e=1.67   < long-simu-testplan.r > power-L16-mn-ar1-500.out
# 
# $R  --vanilla --quite power=0 nloop=10 nrep=1000 phe.sample=500 phe.dist=mn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=type1.L1.mn.ar1.500.1.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > type1-L1-mn-ar1-500-1.out
#
# testCommand();
#---------------------------------------------------
testCommand<-function()
{
	par <- list( 
		a          = 0.5, 
		b          = 0.5, 
		rho        = 0.75,
		sig_a      = 0.8, 
		sig_b      = 0.8, 
		sig_e      = 0.8, 
		beta.c1    = 0.4, 
		times      = 6, 
	    snprange   = c(3*1000, 30*1000),
		a.level    = 10^(-6), 
		w.common   = c(0.5, 0.5),
		w.rare     = c(1,25),
		common.effect = 0.27,
		common.cause = 3,
		rare.ratio = c( 0.4, 0.6 ), 
	    rare.effect= c( NA, 0 ),
	    rare.direct= c( 1, 0 ),
		rare.cutoff= 0.05 );

	power.test  <- as.integer(get_con_param("power"))==1
	ret.rdata   <- as.character(get_con_param("ret.rdata"))
	nloop       <- as.integer(get_con_param("nloop"))
	nrep        <- as.integer(get_con_param("nrep"))
	nsample     <- as.integer(get_con_param("phe.sample"))
	phe.dist    <- as.character(get_con_param("phe.dist"))
	phe.cov     <- as.character(get_con_param("phe.cov"))
	phe.par     <- as.character(get_con_param("phe.par"))

	#power.test <- T;
	#nloop      <- 2;
	#nrep       <- 1000;
	#phe.sample <- 500;
	#phe.dist   <- "msn"
	#phe.cov    <- "ar1"
	#phe.par    <- "0.7,0.5"
	#ret.rdata  <- "test-xxx.rdata"

	pars <- strsplit(phe.par, ",");
	par$par1 <- as.numeric(pars[[1]][1]);
	par$par2 <- as.numeric(pars[[1]][2]);
	par$par3 <- as.numeric(pars[[1]][3]);
	par$par4 <- as.numeric(pars[[1]][4]);
	par$par5 <- as.numeric(pars[[1]][5]);
	par$par6 <- as.numeric(pars[[1]][6]);

	a     <- as.numeric( get_con_param("a") )
	b     <- as.numeric( get_con_param("b") )
	rho   <- as.numeric( get_con_param("rho") )
	sig_a <- as.numeric( get_con_param("sig.a") )
	sig_b <- as.numeric( get_con_param("sig.b") )
	sig_e <- as.numeric( get_con_param("sig.e") )	
	
	if(!is.na(a))     par$a     <- a;
	if(!is.na(b))     par$b     <- b;
	if(!is.na(rho))   par$rho   <- rho;
	if(!is.na(sig_a)) par$sig_a <- sig_a;
	if(!is.na(sig_b)) par$sig_b <- sig_b;
	if(!is.na(sig_e)) par$sig_e <- sig_e;

	if ( power.test )
		power.test( par, nloop, nsample, phe.dist, phe.cov, ret.rdata )
	else
		type1.test( par, nloop, nrep, nsample, phe.dist, phe.cov, ret.rdata);
}
