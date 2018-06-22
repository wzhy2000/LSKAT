longskat_gene_simulate<-function( power.test = TRUE, n.minsect=3000, n.maxsect=30000, n.sample=800, n.time=6, n.gene=10, 
    plink.format = FALSE,
    file.plink.prefix="LSKAT.plink.test",
    geno.miss=0.01,
    pheno.miss=0.10,
	pheno.dist="mn",
	pheno.cov="AR1",
	intercept=FALSE,
	par=list() )
{
	if(missing(par) || length(par)==0 )
		par <- list(b0=1, b1=0.5, b2=0.5, 
		    sig.a=0.8, sig.b=0.8, sig.e=0.8,
		    rho=0.7, 
		    cov.param   = c(0,7, 0.8, 0.2),
		    time.cov    = 0,
		    time.effect = c(0.2, -0.08),
		    max.common.causal  = 4,
		    coef.common.causal = 0.12,
		    max.rare.causal    = 10,
		    coef.rare.causal   = 0.08,
		    positive.ratio     = 1,
		    rare.cutoff        = 0.05 );

	if(!(tolower(pheno.dist)%in% c("mn", "mt", "msn", "mmn")))
		stop("The option 'pheno.dist' has four optional values: mn, mt, msn, mmn.");

	if(!(tolower(pheno.cov) %in% c("ar1", "sad", "cs")))
		stop("The option 'pheno.cov' has three optional values: ar1, sad, cs ");

	pheno.dist <- tolower(pheno.dist);
	pheno.cov <- tolower(pheno.cov);

	snp.mat <- simu_snp_mat( n.minsect, n.maxsect, n.sample, n.gene );

	f.simu <- f.mn.ar1;

	if (pheno.dist=="mn" && pheno.cov=="ar1") f.simu <- f.mn.ar1;
	if (pheno.dist=="mn" && pheno.cov=="sad") f.simu <- f.mn.sad;
	if (pheno.dist=="mn" && pheno.cov=="cs") f.simu  <- f.mn.cs;
	if (pheno.dist=="mt" ) f.simu <- f.mt.ar1;
	if (pheno.dist=="msn" ) f.simu <- f.sn.ar1;
	if (pheno.dist=="mmn" ) f.simu <- f.mmn.ar1;

	if( power.test )
		phe <- simu.long.phe.power( n.sample, n.time, par, f.simu, intercept, snp.mat[[1]] )
	else
		phe <- simu.long.phe.none( n.sample, n.time, par, f.simu, intercept );

	if(geno.miss!=0)
	{
		for(i in 1:n.gene)
		{
			snp.vec <- unlist(snp.mat[[i]]);
			idx.miss <- sample( length(snp.vec) )[ 1:round(runif(1, 0, geno.miss)*length(snp.vec)+1)];
			snp.vec[ idx.miss ] <- NA;
			snp.mat[[i]] <- array(snp.vec, dim=dim(snp.mat[[i]]));
			rownames(snp.mat[[i]]) <- paste("ID", 1:NROW(snp.mat[[i]]), sep="");
		}
	}

	if(pheno.miss!=0)
	{
		y.vec <- unlist(phe$y);
		idx.miss <- sample( length(y.vec) )[ 1:round(runif(1, 0, pheno.miss)*length(y.vec)+1)];
		y.vec[ idx.miss ] <- NA;
		phe$y <- array(y.vec, dim=dim(phe$y));
	}

	colnames(phe$y) <- paste("Y", 1:(NCOL(phe$y)), sep="") ;
	rownames(phe$y) <- paste("ID", 1:NROW(phe$y), sep="");

	colnames(phe$cov) <- paste("X", 1:(NCOL(phe$cov)), sep="") ;
	rownames(phe$cov) <- paste("ID", 1:NROW(phe$cov), sep="");

	if(!plink.format)
		return(list(phe.long = phe$y, phe.cov=phe$cov, snp.mat=snp.mat ))
	else
	{
		snp.chr10 <- c();
		snp.pos10 <- c();
		snp.mat10 <- c();
		for(i in 1:n.gene)
		{
			snp.chr10 <- c( snp.chr10, rep(i, NCOL(snp.mat[[i]])) );
			snp.pos10 <- c( snp.pos10, 1:NCOL(snp.mat[[i]]) );
			#smp.mat[[i]]: [N,P] -> [P,N]
			snp.mat10 <- rbind( snp.mat10, t(snp.mat[[i]]) );
		}

		df.gene<-c();
		for(i in 1:n.gene)
		{
			gene.name <- paste( "Gene", i, sep="");
			snp.name  <- paste( "SNP", i, 1:NCOL(snp.mat[[i]]), sep="-");
			df.gene   <- rbind( df.gene, data.frame(gene=gene.name, snp=snp.name));
		}

		snp.mat10[which(is.na(snp.mat10))] <- -1;
		snp.mat.ext <- cbind(chr=snp.chr10, pos=snp.pos10, snp.mat10);
		rownames(snp.mat.ext) <- df.gene$snp;
		plink <- convert_simpe_to_plink( snp.mat.ext, file.plink.prefix );

		plink$file.gene.set <- paste(file.plink.prefix, "-gene.tab", sep="");
		write.table(df.gene, file=plink$file.gene.set, quote=F, row.names=F, col.names=F, sep=" ");
		cat("Writing gene set file to", plink$file.gene.set,"\n");

		plink$file.phe.cov <- paste(file.plink.prefix, "-cov.csv", sep="");
		write.csv(phe$cov, file=plink$file.phe.cov, quote=F, row.names=T);
		cat("Writing covariate file to", plink$file.phe.cov,"\n");

		plink$file.phe.long <- paste(file.plink.prefix, "-long.csv", sep="")
		write.csv(phe$y, file=plink$file.phe.long, quote=F, row.names=T);
		cat("Writing longitudinal data file to", plink$file.phe.long,"\n");

		plink$phe.long   <- phe$y;
		plink$phe.cov    <- phe$cov;
		plink$snp.mat  <- snp.mat;

		return(plink);
	}
}

simu_snp_mat<-function( min.sect, max.sect, n.sample, mat.count)
{
	file.snp.hap1 <- system.file("extdata", "skat-test-1.hap.gz", package="LSKAT");
	file.snp.pos2 <- system.file("extdata", "skat-test-1.pos", package="LSKAT");

	snp.hap <- read.table(file.snp.hap1, header=F);
	snp.pos <- read.table(file.snp.pos2, header=T);
	snp.maxpos <- max(snp.pos$CHROM_POS);
	snp.minpos <- min(snp.pos$CHROM_POS);

	rare.cutoff <- 1/sqrt(2*n.sample);

	mat.ret <- list();
	n.mat <- 1;
	avoid_range <- c();
	while(n.mat<=mat.count)
	{
		snp.start <- as.integer(runif(1, min=snp.minpos, max=snp.maxpos));
		snp.sect  <- as.integer(runif(1, min=min.sect,   max=max.sect));
		p.sel<-which( snp.pos$CHROM_POS>snp.start & snp.pos$CHROM_POS<(snp.start+snp.sect) );

		if(n.mat==1)
		{
			avoid_range <- c( snp.start-1000, snp.start + snp.sect+1000);
		}
		else
		{
			if( length( intersect(avoid_range[1]:avoid_range[2], snp.start:(snp.start + snp.sect) )	)>0 )
				next;
		}

		snp.pair <- sample(1:(dim(snp.hap)[1]-2));
		snp.mat1 <- snp.hap[ snp.pair[1:n.sample], p.sel + 2 ];
		snp.mat2 <- snp.hap[ snp.pair[(n.sample+1):(2*n.sample)], p.sel + 2 ];
		snp.mat <- snp.mat1 + snp.mat2 -2;

		maf <- colMeans(snp.mat)/2;
		m.same <- which( maf==1 | maf==0 );
		if (length(m.same)>0)
			snp.mat <- snp.mat[, -m.same, drop=F ];

		#make sure of 2: Minor Allels, 0:Major Allels
		maf <- colMeans(snp.mat)/2;
		m.minor <- which( maf>0.5 );
		if (length(m.minor)>0)
			for(i in 1:length(m.minor))
				snp.mat[,m.minor[i]] <- 2 - snp.mat[,m.minor[i]];

		maf <- colMeans(snp.mat)/2;
		m.smallmaf <- which( maf < 5/n.sample);
		if (length(m.smallmaf)>0)
			snp.mat <- snp.mat[, -m.smallmaf,drop=F ];

		if (dim(snp.mat)[2]==1)
			next;

		maf <- colMeans(snp.mat)/2;
		n.rare <- length(which(maf<rare.cutoff));
		if (n.rare<1)
			next;

		if ( n.rare>=1 && n.rare>length(maf)*1/3 )
		{
			n.rare0 <- as.integer(length(maf)*1/3)+1;
			if (n.rare0>1)
			{
				rare.rm <- which(maf<rare.cutoff)[c(n.rare0:n.rare)];
				if (length(rare.rm)>0) snp.mat <- snp.mat[, -rare.rm, drop=F ];
			}
		}


		rownames(snp.mat) <- paste("ID", 1:NROW(snp.mat), sep="");
		mat.ret[[n.mat]] <- snp.mat;

		n.mat <- n.mat + 1;
	}

	return(mat.ret);
}

simu.long.phe.none<-function( n.sample, n.time, par, f.simu, intercept )
{
	cov.mat <- cbind( rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );
	if(intercept)
		y <-  f.simu(n.sample, n.time, par) + cbind(1, cov.mat )%*%c( par$b0, par$b1, par$b2 )%*%t(rep(1, n.time))
	else
		y <-  f.simu(n.sample, n.time, par) + cov.mat %*% c(  par$b1, par$b2 )%*%t(rep(1, n.time));

	return(list(y=y, cov=cov.mat));
}


simu.long.phe.power<-function( n.sample, n.time, par, f.simu, intercept, snp.mat)
{
	y.random <- f.simu(n.sample, n.time, par);

	cov.mat <- cbind( rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	if(intercept)
		y.fixed  <- cbind(1, cov.mat) %*% c(par$b0, par$b1, par$b2 ) %*% t( rep(1, n.time))
	else
		y.fixed  <- cov.mat %*% c(par$b1, par$b2 ) %*% t( rep(1, n.time));

	y.time <- 0;
	if(!is.null(par$time.cov) && par$time.cov>0)
	{
		y.time <- par$time.effect[1]*t(matrix(rep(c(1:n.time),n.sample), nrow=n.time))
		y.time <- y.time + par$time.effect[2]*(t(matrix(rep(c(1:n.time),n.sample), nrow=n.time))**2)
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
				sign <- ifelse( runif(1) <= par$positive.ratio, 1, -1 );
				y.snp <- y.snp + sign * snp.mat[,i] * par$coef.rare.causal *abs(log10(maf[i]));
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
				sign <- ifelse( runif(1) <= par$positive.ratio, 1, -1 );
				y.snp <- y.snp + sign*snp.mat[,i] * par$coef.common.causal;
				snp.comm <- snp.comm+1;
				X.snp <-cbind(X.snp, snp.mat[,i]);
			}
		}
	}

	y.snp <- y.snp%*%t(rep(1, n.time))

	y <- y.random + y.fixed + y.snp + y.time;

	h2 <- var(c(y.snp))/var(c(y.snp + y.random));
	get_stats<-function(y){return(c(mean(y, na.rm=T), min(y, na.rm=T), max(y, na.rm=T), sd(y, na.rm=T)))}
	ys <- rbind(get_stats(y), get_stats(y.fixed), get_stats(y.random), get_stats(y.snp), get_stats(y.time))
	rownames(ys) <- c("Y","Y.fixed","Y.random","Y.snp","Y.time");
	colnames(ys) <- c("Mean","Min","Max","SD");

	show(round(ys,2));
	cat( "h2=", round(h2,2),"Rare=", snp.rare, "\n");

	return(list(y=y, cov=cov.mat, h2=h2, X.snp = X.snp))
}

# time, rho, sig.a, sig.b, sig.e
f.mn.ar1<-function( sample, n.time, par)
{
	ncol <- n.time;

	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);

	sigma.a <- par$sig.a^2
	sigma.b <- par$sig.b^2*AR1;
	sigma.e <- diag(par$sig.e^2, ncol);

	r <- rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +
		 array(rnorm( sample*ncol,  0, par$sig.e ), dim=c(sample, ncol) ) ;

	sigma <- sigma.a + sigma.b + sigma.e;
#cat("MN.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");

	return(r);
}

f.mn.sad<-function( sample, n.time, par )
{
	ncol <- n.time;

	phi  <- par$rho;
	sad1 <- array(1, dim=c(ncol,ncol));
	for(i in 1:ncol)
		for(j in i:ncol)
		{
			sad1[i,j] <- phi^(j-i) * (1-phi^(2*i))/(1-phi^2)
			sad1[j,i] <- phi^(j-i) * (1-phi^(2*i))/(1-phi^2)
		}


	sigma.a <- par$sig.a^2
	sigma.b <- sad1* par$sig.b^2;
	sigma.e <- diag(par$sig.e^2, ncol);

	r <- rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +
		 array(rnorm( sample*ncol,  0, par$sig.e ), dim=c(sample, ncol) ) ;

	sigma <- sigma.a + sigma.b + sigma.e;
#cat("MN.SAD---",min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");

	return(r);
}

f.mn.cs<-function( sample, n.time, par )
{
	ncol <- n.time;
	CM1 <- array(par$rho,dim=c(ncol,ncol));
	for(i in 1:ncol)
		CM1[i,i] <- 1;

	sigma.a <- par$sig.a^2
	sigma.b <- par$sig.b^2 * CM1;
	sigma.e <- diag(par$sig.e^2, ncol);

	r <- rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +
		 array(rnorm( sample*ncol,  0, par$sig.e ), dim=c(sample, ncol) ) ;

	sigma <- sigma.a + sigma.b + sigma.e;
#cat("MN.CM---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");

	return(r);
}

f.mt.ar1<-function( sample, n.time, par )
{
	ncol <- n.time;
	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);

	sigma.a <- par$sig.a^2
	sigma.b <- par$sig.b^2*AR1;
	sigma.e <- diag(par$sig.e^2, ncol);

	r <- rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +
		 array(rt( sample*ncol,  df=10 ), dim=c(sample, ncol) ) ;

	sigma <- sigma.a + sigma.b + sigma.e;
#cat("MT.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");

	return(r);

}

f.sn.ar1<-function( sample, n.time, par )
{
	if(!requireNamespace("sn"))
		stop("Need to install package sn\n");

	ncol <- n.time;

	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);

	sigma.a <- par$sig.a^2
	sigma.b <- par$sig.b^2*AR1;
	sigma.e <- diag(par$sig.e^2, ncol);

	r <- rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +
		 array(sn::rsn( sample*ncol, omega=par$sig.e, alpha = 10 ), dim=c(sample, ncol) ) ;

	sigma <- sigma.a + sigma.b + sigma.e;
#cat("SN.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");
	return(r);
}

#parlist :  rho, sig.a, sig.b, sige_e, par1(sd1), par2(sd2), par3(ratio)

f.mmn.ar1<-function( sample, n.time, par )
{
	ncol <- n.time;

	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);

	sigma.a <- par$sig.a^2
	sigma.b <- par$sig.b^2*AR1;
	sigma.e <- diag(par$sig.e^2, ncol);

	r <- rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +
		 array( rnorm( sample*ncol, sd = par$cov.param[1] )*par$cov.param[3] +
		 	rnorm( sample*ncol, sd = par$cov.param[2] )*(1-par$cov.param[3]), dim=c(sample, ncol) ) ;

	sigma <- sigma.a + sigma.b + sigma.e;
#cat("MMN.AR1---", min(r), max(r), min(sigma), max(sigma),"\n");

	return(r);
}

convert_simpe_to_plink <- function( snp.mat, snp.file.base )
{
	chromosome <- snp.mat[,1];
	position <- snp.mat[,2];

	# PLINK raw data: 1/2/3==> AA,AB,BB, 0==>NA
	snp.mat <- snp.mat[,-c(1,2),drop=F] + 1;

	sub.name <- colnames(snp.mat);
	snp.name <- rownames(snp.mat);

	###snps
	dim.snps <- dim(snp.mat);

	snps <- as.raw( as.matrix(snp.mat ) );
	snps <- array(snps, dim=dim.snps);
	colnames(snps) <- sub.name;
	rownames(snps) <- snp.name;
	class(snps) <- "SnpMatrix";

	r <- write.plink( file.base=snp.file.base, snp.major = F, snps=t(snps),
	    	id=sub.name,
	    	father=rep(0,dim.snps[2]),
	    	mother=rep(0,dim.snps[2]),
	    	sex=rep(0,dim.snps[2]),
	    	phenotype=rep(-9,dim.snps[2]),
			chromosome=chromosome,
			genetic.distance=position,
			position= position,
			allele.1 = rep("A",dim.snps[1]),
			allele.2 = rep("B",dim.snps[1]),
			na.code=0);

	cat("Genotype files have been converted into PLINK binary format(bed/bim/fam)\n");

	return(list(file.plink.bed = paste(snp.file.base, ".bed", sep=""),
   	    	file.plink.bim = paste(snp.file.base, ".bim", sep=""),
   	    	file.plink.fam = paste(snp.file.base, ".fam", sep="")));
}
