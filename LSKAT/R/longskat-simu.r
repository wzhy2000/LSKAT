#public:

longskat_gene_simulate<-function(g.snp.hap1, g.snp.pos2, plink.format=F, file.plink.prefix=NA, n.minsect=3000, n.maxsect=30000, n.sample=800, 
	par=list(a=0.5, b=0.5, sig_a=0.8, sig_b=0.8, sig_e=0.8, 
	    beta.effect= c( 0.1, 0.1, 0, 0.8 ), 
	    c1 = c(0.2828, 0.2828, 0, 0), 
	    times = 6,
	    par1 = 0.7,
	    par2 = 0.4),
		covariance=list(phe.dist="mn", phe.cov="ar1"))
{
	snp.mat <- simu_snp_mat(g.snp.hap1, g.snp.pos2, n.minsect, n.maxsect, n.sample, 10);

	f.simu <- f.mn.ar1;
	if (covariance$phe.dist=="mn" && covariance$phe.cov=="ar1") f.simu <- f.mn.ar1; 
	if (covariance$phe.dist=="mn" && covariance$phe.cov=="sad") f.simu <- f.mn.sad; 
	if (covariance$phe.dist=="mn" && covariance$phe.cov=="cm") f.simu  <- f.mn.cm; 
	if (covariance$phe.dist=="mt" ) f.simu <- f.mt.ar1; 
	if (covariance$phe.dist=="msn" ) f.simu <- f.sn.ar1; 
	if (covariance$phe.dist=="mmn" ) f.simu <- f.mmn.ar1;

	#power test
	phe <- simu.long.phe.power( n.sample, snp.mat[[1]], par, f.simu);
	
	#type1 error test
	#phe <- simu.long.phe.none( n.sample, par, f.simu);

	if(!plink.format)
		return(list(phe.y=phe$y, phe.cov=phe$y.cov, snp=snp.mat[[1]]))
	else
	{	
		snp.chr10 <- c();
		snp.pos10 <- c();
		snp.mat10 <- c();
		for(i in 1:10)
		{
			snp.chr10 <- c( snp.chr10, rep(i, NCOL(snp.mat[[i]])) );
			snp.pos10 <- c( snp.pos10, 1:NCOL(snp.mat[[i]]) );
			#smp.mat[[i]]: [N,P] -> [P,N]
			snp.mat10 <- rbind( snp.mat10, t(snp.mat[[i]]) );
		}

		df.gene<-c();
		for(i in 1:10)
		{
			gene.name <- paste( "TGene", i, sep="");
			snp.name  <- paste( "G", i, 1:NCOL(snp.mat[[i]]), sep="-");
			df.gene   <- rbind( df.gene, data.frame(gene=gene.name, snp=snp.name));
		}	

		colnames(snp.mat10) <- 1:NCOL(snp.mat10);
		snp.mat.x <- cbind(chr=snp.chr10, pos=snp.pos10, snp.mat10);
		rownames(snp.mat.x) <- df.gene$snp;
		plink <- convert_simpe_to_plink( snp.mat.x, file.plink.prefix );
		
		plink$file.phe.cov <- paste(file.plink.prefix, "-cov.csv", sep="");
		write.csv(phe$y.cov, file=plink$file.phe.cov, quote=F, row.names=F);
		plink$file.phe.long <- paste(file.plink.prefix, "-long.csv", sep="")
		write.csv(phe$y, file=plink$file.phe.long, quote=F, row.names=F);
		plink$file.gene.set <- paste(file.plink.prefix, "-gene.tab", sep="");
		
		write.table(df.gene, file=plink$file.gene.set, quote=F, row.names=F, col.names=F, sep=" ");
		
		return(plink);
	}
}

simu_snp_mat<-function(snp.file1, snp.file2, min.sect, max.sect, n.sample, mat.count)
{
	snp.hap <- read.table(snp.file1, header=F);
	snp.pos <- read.table(snp.file2, header=T);
	snp.maxpos <- max(snp.pos$CHROM_POS);
	snp.minpos <- min(snp.pos$CHROM_POS);

	rare.cutoff <- 1/sqrt(2*n.sample);
	
	mat.ret <- list();	
	n.mat <- 1;
	while(n.mat<=mat.count)
	{
		snp.start <- as.integer(runif(1, min=snp.minpos, max=snp.maxpos));
		snp.sect  <- as.integer(runif(1, min=min.sect,   max=max.sect));
		
		p.sel<-which( snp.pos$CHROM_POS>snp.start & snp.pos$CHROM_POS<(snp.start+snp.sect) );

		snp.pair <- sample(1:(dim(snp.hap)[2]-2));
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
		
		mat.ret[[n.mat]] <- snp.mat;
		n.mat <- n.mat + 1;
	}
	
	return(mat.ret);
}

simu.long.phe.none<-function( n.sample, par, f.simu )
{
	par_cov <- cbind(1, rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	y <-  f.simu(n.sample, par) + par_cov%*%c( 1, par$a, par$b )%*%t(rep(1, par$times));

	return(list(y=y, cov=par_cov[,-1]))
}

simu.long.phe.power<-function( n.sample, snp.mat, par, f.simu )
{
	par_cov <- cbind(1, rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	y <-  f.simu(n.sample, par) + par_cov%*%c(1, par$a, par$b )%*%t(rep(1, par$times));

	rare.cutoff <- 1/sqrt(2*n.sample);
	maf <- colMeans(snp.mat)/2;
	n.snp <- length(maf);
	snp.rare.total <- length(which(maf<rare.cutoff))
	
	beta.ratio <- cumsum( par$beta.effect );
	snp.comm <- 0;
	snp.rare <- 0;
	for( i in 1:n.snp)
	{
		if (maf[i]<rare.cutoff)
		{
			beta <- par$c1[ which.max(snp.rare/snp.rare.total <= beta.ratio) ]
			if (is.na(beta)) beta<-runif(1, 0.15, 0.2461);
			y <- y + snp.mat[,i] * beta;
			snp.rare <- snp.rare + 1;
		}
		else
		{
			if (snp.comm < 2 )
			{
				y <- y + snp.mat[,i] * 0.27;
				snp.comm <- snp.comm+1;
			}
		}
	}

	#adding shareid columns
	y <- cbind(1:NROW(y), y);
	colnames(y) <- c("shareid", paste("Y", 1:(NCOL(y)-1), sep="") );
    y.cov <- cbind(1:NROW(y), par_cov[,-1]);
	colnames(y.cov) <- c("shareid", paste("COV", 1:(NCOL(y.cov)-1), sep=""));
	
	return(list(y=y, y.cov=y.cov));
}

# time, rho, sig_a, sig_b, sig_e
f.mn.ar1<-function( sample, par)
{
	ncol <- par$times;

	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);
		
	sigma.a <- par$sig_a^2
	sigma.b <- par$sig_b^2*AR1;
	sigma.e <- diag(par$sig_e^2, ncol);
	
	r <- rnorm( sample,  0, par$sig_a ) %*% array(1, dim=c(1,ncol))+    
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +  
		 array(rnorm( sample*ncol,  0, par$sig_e ), dim=c(sample, ncol) ) ;
	
	sigma <- sigma.a + sigma.b + sigma.e;
#cat("MN.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	
	
	return(r);
}

f.mn.sad<-function( sample, par )
{
	ncol <- par$times;

	phi  <- par$rho;
	sad1 <- array(1, dim=c(ncol,ncol));
	for(i in 1:ncol)
		for(j in i:ncol)
		{
			sad1[i,j] <- phi^(j-i) * (1-phi^(2*i))/(1-phi^2)				
			sad1[j,i] <- phi^(j-i) * (1-phi^(2*i))/(1-phi^2)				
		}


	sigma.a <- par$sig_a^2
	sigma.b <- sad1* par$sig_b^2;
	sigma.e <- diag(par$sig_e^2, ncol);

	r <- rnorm( sample,  0, par$sig_a ) %*% array(1, dim=c(1,ncol))+    
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +  
		 array(rnorm( sample*ncol,  0, par$sig_e ), dim=c(sample, ncol) ) ;
	
	sigma <- sigma.a + sigma.b + sigma.e;
cat("MN.SAD---",min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	

	return(r);
}

f.mn.cm<-function( sample, par )
{
	ncol <- par$times;
	CM1 <- array(par$rho,dim=c(ncol,ncol));
	for(i in 1:ncol)
		CM1[i,i] <- 1;
		
	sigma.a <- par$sig_a^2
	sigma.b <- par$sig_b^2 * CM1;
	sigma.e <- diag(par$sig_e^2, ncol);

	r <- rnorm( sample,  0, par$sig_a ) %*% array(1, dim=c(1,ncol))+    
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +  
		 array(rnorm( sample*ncol,  0, par$sig_e ), dim=c(sample, ncol) ) ;
	
	sigma <- sigma.a + sigma.b + sigma.e;
cat("MN.CM---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	

	return(r);
}

f.mt.ar1<-function( sample, par )
{
	ncol <- par$times;
	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);
		
	sigma.a <- par$sig_a^2
	sigma.b <- par$sig_b^2*AR1;
	sigma.e <- diag(par$sig_e^2, ncol);
	
	r <- rnorm( sample,  0, par$sig_a ) %*% array(1, dim=c(1,ncol))+    
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +  
		 array(rt( sample*ncol,  df=10 ), dim=c(sample, ncol) ) ;
	
	sigma <- sigma.a + sigma.b + sigma.e;
cat("MT.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	

	return(r);

}

f.sn.ar1<-function( sample, par )
{
	library(sn);
	
	ncol <- par$times;
	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);
		
	sigma.a <- par$sig_a^2
	sigma.b <- par$sig_b^2*AR1;
	sigma.e <- diag(par$sig_e^2, ncol);
	
	r <- rnorm( sample,  0, par$sig_a ) %*% array(1, dim=c(1,ncol))+    
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +  
		 array(rsn( sample*ncol, omega=par$sig_e, alpha = 10 ), dim=c(sample, ncol) ) ;
	
	sigma <- sigma.a + sigma.b + sigma.e;
cat("SN.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	
	return(r);
}

#parlist : times, rho, sig_a, sig_b, sige_e, par1(sd1), par2(sd2), par3(ratio) 

f.mmn.ar1<-function( sample, par )
{
	ncol <- par$times;

	AR1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1[i,j] <- par$rho^abs(i-j);

	sigma.a <- par$sig_a^2
	sigma.b <- par$sig_b^2*AR1;
	sigma.e <- diag(par$sig_e^2, ncol);
	
	r <- rnorm( sample,  0, par$sig_a ) %*% array(1, dim=c(1,ncol))+    
		 rmvnorm( sample,  rep(0, ncol), sigma.b ) +  
		 array( rnorm( sample*ncol, sd = par$par1 )*par$par3 + 
		 	rnorm( sample*ncol, sd = par$par2 )*(1-par$par3), dim=c(sample, ncol) ) ;

	sigma <- sigma.a + sigma.b + sigma.e;
cat("MMN.AR1---", min(r), max(r), min(sigma), max(sigma),"\n");		
	
	return(r);
}

