library(LSKAT);
library(SKAT);
library(mvtnorm);
library(parallel)
library(nlme)

g.snp.hap1  <- "skat-test-1.hap";
g.snp.pos2  <- "skat-test-1.pos"

get_con_param<-function(parm.id)
{
	for (e in commandArgs())
	{
		ta = strsplit(e,"=", fixed=TRUE);
		if(! is.na( ta[[1]][2]))
		{
			temp = ta[[1]][2];
			if( ta[[1]][1] == parm.id) {
				temp = as.character(temp);
				return (temp);
			}
		}
	}

	return(NA);
}


if(0)
{
 library(mvtnorm)
 sample<- 200;
 par <- list(times=6,  rho=0.7, sig_a=0.2, sig_b=0.8, sig_e=0.6, par1=0.2, par2=1.2, par3=0.8)
 dim( f.mn.ar1 ( sample, par ))
 dim( f.mn.sad ( sample, par ))
 dim( f.mn.cm  ( sample, par ))
 dim( f.mt.ar1 ( sample, par ))
 dim( f.sn.ar1 ( sample, par ))
 dim( f.mmn.ar1( sample, par ))
}

# time,s rho, sig_a, sig_b, sig_e
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
	
	#sigma <- sigma.a + sigma.b + sigma.e;
	#r <- rmvnorm( sample,  rep(0, ncol), sigma )  ;

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
		 array(rsn( sample*ncol, omega=par$sig_e, alpha = 40 ), dim=c(sample, ncol) ) ;
		 #array(rsn( sample*ncol, omega=par$sig_e, alpha = 1/4 ), dim=c(sample, ncol) ) ;
	
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


simu.long.phe.none<-function( n.sample, par, f.simu )
{

	par_cov <- cbind(rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	if(par$intercept)
		y.fixed  <- cbind(1, par_cov) %*% c(par$mu, par$a, par$b ) %*% t( rep(1, par$times))
	else
		y.fixed  <- par_cov %*% c(par$a, par$b ) %*% t( rep(1, par$times));

	y.random <- f.simu(n.sample, par);
	 
	y <- y.fixed + y.random;
	
	return(list(y=y, cov=par_cov))
}


check.type1err<-function( par, n.loop, n.rep, n.sample, f.simu, ncores=1)
{
	snp.mat  <- simu.snp.mat(g.snp.hap1, g.snp.pos2, par$snprange, n.sample, n.rep, par$rare.cutoff, par$rare.count.range, par$common.count.range,  ncores)

	rs <- mclapply(1:n.loop, function(i)
	{
		set.seed( runif(1)*1000 + i + proc.time()[3] );

		phe <- simu.long.phe.none( n.sample, par, f.simu);

		obj.h0 <- try( longskat_est_model( phe$y, phe$cov, y.cov.time=par$y.cov.time, g.maxiter=1, intercept=par$intercept ));
		
		if(class(obj.h0) != "try-error")
			rlist <- c( n.sample, obj.h0$par$sig_a, obj.h0$par$sig_b, obj.h0$par$sig_e, 
						obj.h0$par$rho,  
						obj.h0$par$mu, 
						obj.h0$par$par_cov, 
						obj.h0$likelihood)
		else
		{
			cat("--------------Failed to estimate the parameters in the model.\n");
			next;
		}
		
		rlist <- list();
		for(j in 1:n.rep)
		{
			r.lskat <- try( longskat_gene_run( obj.h0, snp.mat[[j]], 
					weights.rare = par$w.rare, 
					weights.common = par$w.common, 
					rare.cutoff = par$rare.cutoff,
					test.type = par$test.type,
					run.cpp = F,
					r.corr = 2) );

			if( is.null(r.lskat) || class(r.lskat)=="try-error" )
				next
			else
			{
				cat("TYPE1.i=", i, j, n.sample, par$times, r.lskat$p.lskat, r.lskat$p.burden, r.lskat$snp.total, r.lskat$snp.rare, 
							obj.h0$par$sig_a, obj.h0$par$sig_b, obj.h0$par$sig_e, obj.h0$par$rho, obj.h0$par$mu, obj.h0$par$par_cov, obj.h0$likelihood, "\n\n");

				rlist[[j]] <- c(i, j, n.sample, par$times, r.lskat$p.lskat, r.lskat$p.burden, r.lskat$snp.total, r.lskat$snp.rare, 
						obj.h0$par$sig_a, obj.h0$par$sig_b, obj.h0$par$sig_e, obj.h0$par$rho, obj.h0$par$mu, obj.h0$par$par_cov, obj.h0$likelihood );
			}
		}
		
		return( do.call(rbind, rlist) );
	}, mc.cores=ncores)
	
	rs <- do.call(rbind, rs);
	
	colnames(rs)<-c("i","j", "sample","times","p.lskat", "p.burden", "snp","rare", "sig_a", "sig_b","sig_e", "rho", "mu", "a", "b", "LR")
	
	save(rs, file=paste("simu-type1-", "0.rdata", sep="") );
	
	return(rs);
}

check.power<-function( par, n.loop, n.sample, f.simu, ncores=1 )
{
	snp.mat <- simu.snp.mat(g.snp.hap1, g.snp.pos2, par$snprange, n.sample, n.loop, par$rare.cutoff, par$rare.count.range, par$common.count.range, ncores)

	rs <- mclapply(1:n.loop, function(i)
	{
		rlist <- c();
		while(length(rlist)<=0)
		{
			phe <- simu.long.phe.power( n.sample, snp.mat[[i]], par, f.simu);

			obj.h0 <- try( longskat_est_model( phe$y, phe$cov, y.cov.time=par$y.cov.time, g.maxiter=1, intercept=par$intercept, method="REML" ));

			if(class(obj.h0) != "try-error"  && obj.h0$bSuccess)
				rlist <- c( n.sample, obj.h0$par$sig_a, obj.h0$par$sig_b, obj.h0$par$sig_e, 
						obj.h0$par$rho,  
						obj.h0$par$mu, 
						obj.h0$par$par_cov[c(1,2)], 
						obj.h0$likelihood)
			else
				next; 
		}

		if(par$intercept) 
			model <- "Y ~  1 + X1 + X2"
		else	
			model <- "Y ~  X1 + X2";

		if(!is.null(par$y.cov.time) && par$y.cov.time>0)
			model <- paste(model, "+ T + T2");

cat("formula=", model, "\n");

		y.phe <- cbind(subjects=seq(1, n.sample), phe$cov, Y=phe$y );
		y.df <- c();

		if(!is.null(par$y.cov.time) && par$y.cov.time>0)
		{		
			for(t in 1:par$times)
		        y.df <- rbind(y.df, cbind(y.phe[,c(1, 2,3, t+3)], t, t**2, t));
			colnames(y.df) <- c("subjects", "X1", "X2", "Y", "T", "T2", "time");	
		}
		else
		{		
			for(t in 1:par$times)
		        y.df <- rbind(y.df, cbind(y.phe[,c(1, 2,3, t+3)], t));
			colnames(y.df) <- c("subjects", "X1", "X2", "Y", "time");	
		}

		y.df <- as.data.frame(y.df);
        
		#ylm <- try ( lme(as.formula(model), random = ~ time| subjects + 1 | subjects,  correlation = corAR1(),  data=y.df) )
		ylm <- try ( lme(as.formula(model), random = ~ time| subjects,  correlation = corAR1(),  data=y.df) )
		if( class(ylm)!="try-error" )
			cat("***LME =", round(ylm$coefficients$fixed,3) , "\n")
		else
			cat("***LME =Unavalaible.\n");

		# LSKAT--gene 
		r.lskat <- try( longskat_gene_run( obj.h0, snp.mat[[i]], 
				weights.rare = par$w.rare, 
				weights.common = par$w.common, 
				rare.cutoff = par$rare.cutoff,
				test.type = par$test.type,
				run.cpp = F,
				r.corr = 2) );
		if( !is.null(r.lskat) && class(r.lskat)!="try-error" )
			rlist <- c(rlist, r.lskat$snp.total, r.lskat$snp.rare, r.lskat$p.lskat, r.lskat$p.burden )
		else
			rlist <- c(rlist, NA, NA, NA, NA)

		# LSKAT--SNP 
		if(0)
		{
			ls0.res <- c();
			for(k in 1:NCOL(snp.mat[[i]]))
			{
				snp.vec<- c(snp.mat[[i]][,k]);
				ls0 <- longskat_snp_run( obj.h0, 
								snp.vec,
								weights.rare=par$w.rare, 
								weights.common=par$w.common, 
								run.cpp=F,
								rare.cutoff = par$rare.cutoff) ;
				if(class(ls0)!="try-error")
					ls0.res <- rbind( ls0.res, c(k, ls0$maf, ls0$nmiss, ls0$snp.rare, ls0$pv));
			}


			if(NROW(ls0.res)>0)
			{
				snp.sig <- which.min(ls0.res[,5]);
				rlist <- c(rlist,  ls0.res[snp.sig, 5]*NCOL(snp.mat[[i]]), ls0.res[snp.sig, c(1,2) ] )
			}
			else
				rlist <- c(rlist,  NA, NA, NA )
		}
		else
			rlist <- c(rlist,  NA, NA, NA )
		
		# SKAT( Baseline ) 
		phe.bl  <- data.frame(phe$y[,1], phe$cov) ;
		colnames(phe.bl)<-c("Y","X1","X2");
		if(!is.null(par$y.cov.time) && par$y.cov.time>0)
			phe.bl = cbind(phe.bl, T=mean(1:par$times), T2=mean(c(1:par$times)**2) );

		obj.bl  <- SKAT_Null_Model(as.formula(model), data=phe.bl, out_type="C");

		# SKAT--Baseline--SKAT
		r.bl.skat <- try( SKAT_CommonRare( 
				as.matrix(snp.mat[[i]]), 
				obj.bl, 
				weights.beta.rare = par$w.rare, 
				weights.beta.common = par$w.common,
				method = "C", 
				r.corr.rare = 0, 
				r.corr.common = 0, 
				test.type = par$test.type,
				CommonRare_Cutoff = par$rare.cutoff ) );
		if(class(r.bl.skat)!="try-error")
			rlist <- c(rlist, r.bl.skat$p.value)
		else
			rlist <- c(rlist, NA)

		# SKAT--Baseline--Burden
		r.bl.burden <- try( SKAT_CommonRare( 
				as.matrix(snp.mat[[i]]), 
				obj.bl, 
				weights.beta.rare  = par$w.rare, 
				weights.beta.common = par$w.common,
				method = "C", 
				r.corr.rare = 1, 
				r.corr.common = 1, 
				test.type = par$test.type,
				CommonRare_Cutoff = par$rare.cutoff ));
		if(class(r.bl.burden)!="try-error")
			rlist <- c(rlist, r.bl.burden$p.value)
		else
			rlist <- c(rlist, NA )

		# SKAT( MEAN ) 
		mu      <- apply( phe$y, 1, mean );
		phe.mu  <- data.frame(mu, phe$cov) ;
		colnames(phe.mu)<-c("Y","X1","X2");
		if(!is.null(par$y.cov.time) && par$y.cov.time>0)
			phe.mu = cbind(phe.mu, T=mean(1:par$times), T2=mean(c(1:par$times)**2) )
		
		obj.mu  <- SKAT_Null_Model(as.formula(model), data=phe.mu, out_type="C");
		
		# SKAT--Mean--SKAT
		r.mu.skat <- try( SKAT_CommonRare( 
				as.matrix(snp.mat[[i]]), 
				obj.mu, 
				weights.beta.rare = par$w.rare, 
				weights.beta.common = par$w.common,
				method = "C", 
				r.corr.rare = 0, 
				r.corr.common = 0, 
				test.type = par$test.type,
				CommonRare_Cutoff = par$rare.cutoff ) )
		if(class(r.mu.skat)!="try-error")
			rlist <- c(rlist, r.mu.skat$p.value)
		else
			rlist <- c(rlist, NA)
			
		# SKAT--Mean--Burden
		r.mu.burden <- try( SKAT_CommonRare( 
				as.matrix(snp.mat[[i]]), 
				obj.mu, 
				weights.beta.rare  = par$w.rare, 
				weights.beta.common = par$w.common,
				method = "C", 
				r.corr.rare = 1, 
				r.corr.common = 1, 
				test.type = par$test.type,
				CommonRare_Cutoff = par$rare.cutoff ) );
		if(class(r.mu.burden)!="try-error")
			rlist <- c(rlist, r.mu.burden$p.value)
		else
			rlist <- c(rlist, NA)
				
		cat("POWER.i=", i, round(rlist[1:11],3), "\n\t", rlist[-c(1:11)], "\n\n");
		
		#rs <- rbind( rs, c(i, rlist ) );
		return(c(i, rlist ));
	}, mc.cores=ncores);

	rs <- do.call(rbind, rs);

	colnames(rs)<- c("idx", "sample", "sig_a", "sig_b", "sig_e", "rho",  "mu", "Xa", "Xb", "LR",
					"lskat.snp.total", "lskat.snp.rare", "lskat.pv", "lskat.burden.pv", 
					"lskat.snp.pv", "lskat.snp.loc", "lskat.snp.maf",
					"skat.bl.pv", "burden.bl.pv",  "skat.mu.pv", "burden.mu.pv");
	return(rs);
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

	# 40% causual RARE, 50% positive effect
	par$rare.ratio <- c( 0.4*0.5, 0.4*0.5, 0.6 ); 
	par$rare.effect<- c( NA,      NA,      0 );
	par$rare.direct<- c( 1,       -1,      0 );
	r.3 <-check.power( par, nloop, nsample, f.simu, ncores )

	save(par, r.1, r.2, r.3, file=power.rdata);
	
	# 40% causual RARE, 100% positive effect
	par$rare.ratio <- c( 0.4, 0.6 ); 
	par$rare.effect<- c( NA,  0   );
	par$rare.direct<- c( 1,   0 );
	r.1 <- check.power( par, nloop, nsample, f.simu, ncores )

	save(par, r.1, r.2, r.3, file=power.rdata);

	# 40% causual RARE, 80% positive effect
	par$rare.ratio <- c( 0.4*0.8, 0.4*0.2, 0.6 ); 
	par$rare.effect<- c( NA,      NA,      0 );
	par$rare.direct<- c( 1,       -1,      0 );
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

type1.test<-function(par, nloop, nrep, nsample, phe.dist, phe.cov, type1.rdata, ncores)
{
	f.simu <- f.mn.ar1;
	if (phe.dist=="mn" && phe.cov=="ar1") f.simu <- f.mn.ar1; 
	if (phe.dist=="mn" && phe.cov=="sad") f.simu <- f.mn.sad; 
	if (phe.dist=="mn" && phe.cov=="cm") f.simu  <- f.mn.cm; 
	if (phe.dist=="mt" ) f.simu <- f.mt.ar1; 
	if (phe.dist=="msn" ) f.simu <- f.sn.ar1; 
	if (phe.dist=="mmn" ) f.simu <- f.mmn.ar1; 

	cat(" [Type1 Test]\n")
	cat(" * loop=", nloop, "\n")
	cat(" * rep=", nrep, "\n")
	cat(" * nsample=", nsample, "\n")
	cat(" * par=",  unlist(par), "\n")
	cat(" * phe.dist=", phe.dist, "\n")
	cat(" * phe.cov=", phe.cov, "\n")
	cat(" * ret.rdata=", type1.rdata, "\n")
	
	r.1  <- check.type1err( par, nloop, nrep, nsample, f.simu, ncores )
	save(par, r.1, file=type1.rdata);

	a5.level<-10^(-5);
	a4.level<-10^(-4);
	a3.level<-10^(-3);
	a2.level<-10^(-2);

	r.type1<-c(nloop, nrep, nsample, 
			  mean(r.1[,7]),
			  mean(r.1[,8]),
	          length(which( r.1[,5]<a2.level)), 
	          length(which( r.1[,5]<a3.level)), 
	          length(which( r.1[,5]<a4.level)), 
	          length(which( r.1[,5]<a5.level)) );

	names(r.type1) <- c("loop", "repi", "sample", "snp.total", "snp.rare", "a2", "a3", "a4", "a5");
	save(r.type1, par, r.1, file=type1.rdata);
}
