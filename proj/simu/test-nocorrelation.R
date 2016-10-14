longskat_est_model<-function( y.long, y.cov, y.time = NULL, y.cov.time=0, intercept=FALSE, g.maxiter=20, par.init=list(), debug=F, method=c("REML", "ML") )
{
	if(missing(method)) method <- "REML";
	if(missing(debug)) debug <- F;

	na.idx <- which(is.na(rowSums(y.cov)));
	if(length(na.idx)>0)
	{
		y.long <- y.long[-na.idx,,drop=F];
		y.cov  <- y.cov[-na.idx,,drop=F]; 
		if(!is.null(y.time)) y.time <- y.time[-na.idx,,drop=F];
	}	

	ret <- longskat_est_ML( y.long, y.cov, y.time, y.cov.time, intercept, g.maxiter, par.init, debug );
		
	return(ret);	
}

longskat_est_ML<-function( y.long, y.cov0, y.time = NULL, y.cov.time=0, intercept=TRUE, g.maxiter=20, init.par=list(), debug=F )
{
	# check id matched before call here!
	
	y.cov <- y.cov0;
	
	if(is.data.frame(y.long))  y.long <- data.matrix(y.long)
	if(is.data.frame(y.cov))  y.cov <- data.matrix(y.cov)
	if(is.data.frame(y.time)) y.time<- data.matrix(y.time)

	
	ncol <- NCOL(y.long);
	nrow <- NROW(y.long);
	nCov <- NCOL(y.cov);

	if( is.null(y.time))
		y.time <- t(matrix(rep(c(1:ncol),nrow), nrow=ncol)) ;
	
	if(intercept) 
		y.cov <- cbind(1, y.cov)

	X <-  kronecker(y.cov, rep(1, ncol) )

	if( y.cov.time > 0 ) 
		for(i in 1:y.cov.time)
			X <- cbind(X, c(t(y.time))**i)	
#show(head(y.long));
#show(head(y.cov));
#show(head(y.time));

	get_par<-function(par, y, y.time, y.cov )
	{	
		par_rho=par[1];
		sig_a<-par[2];
		sig_b=par[3];
		sig_e<-par[4];

		AR.1 <- array(1:ncol, dim=c(ncol,ncol))
		AR.1 <- par_rho^abs(AR.1-t(AR.1));
		sigma <- diag(sig_e^2, ncol) + sig_a^2 + sig_b^2*AR.1;

		par_cov <- par[ c(5:(4+NCOL(y.cov))) ];
		y.delt <- y.long - t( array( X %*% par_cov, dim=c(NCOL(y.long), NROW(y.long))))
		
		y.mean <- mean(y.delt)
		y.delt - y.delt - y.mean;
	
		A <- 0;
		if(any(is.na(y.delt)))
		{
			for(i in 1:NROW(y.delt) )
			{
				t.sel <- which( !is.na(c(y.time[i,]) ) );
				sig <- sigma[t.sel, t.sel, drop=F];
				A <- A + sum( dmvnorm( y.delt[i,t.sel,drop=F], rep(0, length(t.sel)), sig, log=T ) );
			}
		}
		else
			A <- sum( dmvnorm( y.delt,  rep(0, ncol), sigma, log=T  ) );

		ret <- -A - log(abs(y.mean))* log(NCOL(y.long))

		return( ret );
	}
	
	est_par_cov<-function()
	{
		#1: rho, 2:sig_a, 3:sig_b, 4:sig_e
		par.init <- c( 0.5, rep( sd(y.long, na.rm=T)/3,3) );

		if(intercept)
			par.init <- c(par.init, mean(y.long, na.rm=T));

		for(i in 1:nCov)
			par.init <- c( par.init, 1/( mean(y.cov[,i], na.rm=T)^2 + 1) );

		if (y.cov.time>0)
		{
			par.init <- c( par.init, rep( 1/( mean(y.time, na.rm=T)^2 +1) , y.cov.time ) )
			if(is.null(y.time))
				y.time <- t( t(ifelse(is.na(y.long),NA, 1))*(rep(1:NROW(y.long))))
		}

		return(par.init);			
	}

	# par[1] = rho
	# par[2] = sig_a
	# par[3] = sig_b
	# par[4] = sig_e
	# par[5] = par_u
	# par[6, 5+nCov] = par_cov
	# par[6+nCov, (6:7+nCov) ] = par_cov_time
	# e.g.
	# par.init <- c(rho=0.75, sig_a=0.2, sig_b=0.3, sig_e=0.1, u=1, a=0.5, b=0.5); 
	
	tolerance <- 1;
	loop.n    <- 0;
	min.val   <- Inf;
	par.init  <- min.par <- est_par_cov();

	while( loop.n <= g.maxiter && tolerance > 1e-5 )
	{
		r0 <- try( optim( par.init, get_par, y = y.long, y.time=y.time, y.cov=y.cov, method = "BFGS", control=list(maxit=500) ), silent = F );

		if (class(r0)=="try-error")
		{
			if (min.val >= 1e8)
				loop.n <- loop.n + 0.2;
			par.init <- min.par*runif( length(min.par) );
			next;
		}

		loop.n <- loop.n+1;
		if ( r0$convergence==0 && r0$val<min.val)
		{
			if(debug) cat("  LOOP =", loop.n, "/", g.maxiter, " val=", r0$val,  "par=", r0$par, "\n");
			
			tolerance <- max( c( min.val, par.init) - c(r0$value, r0$par) );
			
			min.val <- r0$value;
			min.par <- r0$par;
			par.init<- min.par;
		}

		par.init <- par.init*runif( length(par.init), 0.8, 1.2 );
	}

	if(debug) cat("  Final =", loop.n, " min.val=", min.val, "par=", min.par, "\n");
	
	if( is.infinite( min.val) || is.na(min.val) || any(is.na(min.par))  )
		return( list(bSuccess=F) );
	

cat("SIG_A/B/E/R=", min.par[c(1:4)], "COV=", min.par[-c(1:4)], "LIKELIHOOD=", min.val, "\n"); 

	par_rho <- abs(min.par[1]);
	sig_a   <- abs(min.par[2]);;
	sig_b   <- abs(min.par[3]);;
	sig_e   <- abs(min.par[4]);
	par_cov <- min.par[-c(1:4)];
	par_mu  <- NA;
	par_t   <- NA;
	
	if(intercept)
	{
		par_mu  <- par_cov[1];
		par_cov <- par_cov[-1];
	}
	
	if( y.cov.time > 0)
	{
		par_t <- par_cov[-c(1:NCOL(y.cov))];
		par_cov <- par_cov[ c(1:NCOL(y.cov))];
	}

	y.delt <- y.long - t( array( X %*% min.par[-c(1:4)], dim=c(NCOL(y.long), NROW(y.long))))

cat("y.delt.mean=", mean(y.delt), "\n"); 

	pars <- list( intercept=intercept,
				  y.cov.time = y.cov.time,
			  	  mu     = par_mu, 
			  	  rho    = par_rho, 
			  	  sig_a  = sig_a, 
			  	  sig_b  = sig_b, 
			  	  sig_e  = sig_e, 
			  	  par_cov= par_cov,  
			  	  par_t  = par_t );  

	r.model <- list(par = pars, likelihood = min.val, y.delt=y.delt, y.time = y.time, y.cov = y.cov0, bSuccess=T );
	
	class(r.model) <- "LSKAT.null.model";

	return(r.model);
}

