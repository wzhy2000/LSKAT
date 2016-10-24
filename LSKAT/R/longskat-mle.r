#solve the problem 'Lapack routine dgesv: system is exactly singular: U[7,7] = 0'
Get_SKAT_Residuals.Get_X1 = function(X1){
	
	qr1<-qr(X1)
	q1<-ncol(X1)
	if(qr1$rank < q1){
		X1.svd<-svd(X1)
		X1 = X1.svd$u	
	} 

	return(X1)
}

#public

longskat_est_model<-function( phe.long, phe.cov, 
				phe.time = NULL, 
				time.cov = 0,
                intercept = FALSE, 
                method = c("REML", "ML"), 
                g.maxiter =  20, 
                par.init = list(), 
                verbose = F )
{
	y.long <- phe.long;
	y.cov  <- phe.cov;
	y.time <- phe.time;
	
	if(missing(method)) method <- "REML";
	if(missing(verbose)) verbose <- F;

	na.rm.id <- c();
	na.idx1 <- c(na.rm.id, which(is.na(rowSums(y.cov))) );
	if(length(na.idx1)>0)
	{
		na.rm.id <- c(na.rm.id, rownames(y.cov)[na.idx1] ); 
		y.long <- y.long[-na.idx1,,drop=F];
		y.cov  <- y.cov[-na.idx1,,drop=F]; 
		if(!is.null(y.time)) y.time <- y.time[-na.idx1,,drop=F];
	}	

	y.long.na <- y.long;
	y.long.na[!is.na(y.long.na)] <- 1;
	y.long.na[is.na(y.long.na)] <- 0;
	na.idx2 <- which(rowSums(y.long.na)==0);
	if(length(na.idx2)>0)
	{
		na.rm.id <- c(na.rm.id, rownames(y.long)[na.idx2]); 
		y.long <- y.long[-na.idx2,,drop=F];
		y.cov  <- y.cov[-na.idx2,,drop=F]; 
		if(!is.null(y.time)) y.time <- y.time[-na.idx2,,drop=F];
	}	

	if(method=="REML")
		ret <- longskat_est_REML( y.long, y.cov, y.time, time.cov, intercept, g.maxiter, par.init, verbose )
	else	
		ret <- longskat_est_ML( y.long, y.cov, y.time, time.cov, intercept, g.maxiter, par.init, verbose );
		
	ret$na.rm.id <- na.rm.id;
	
	if(ret$bSuccess)
		return(ret)
	else
		return(NULL);
}

# y.long format:  <shareid> trait1, ..., traitN
# y.time format: <shareid> time1, ..., timeN
# y.cov format:  <shareid> cov1, ..., covM

longskat_est_REML<-function( y.long, y.cov, y.time = NULL, time.cov = 0, intercept=FALSE, g.maxiter=20, par.init=list(), verbose=F )
{
	# check id matched before call here!

	if(is.data.frame(y.long))  y.long <- data.matrix(y.long)
	if(is.data.frame(y.cov))   y.cov <- data.matrix(y.cov)
	if(is.data.frame(y.time))  y.time<- data.matrix(y.time)
	
	ncol <- NCOL(y.long);
	nrow <- NROW(y.long);
	
	if( is.null(y.time))
		y.time <- t(matrix(rep(c(1:ncol),nrow), nrow=ncol)) ;

if(verbose)
{
	cat("y.cov", "ncol=", ncol, "intercept=", intercept,  "\n");
	show(head(y.cov));
}

	if(intercept) 
		X <-  kronecker(cbind(1, y.cov), rep(1, ncol) )
	else
		X <-  kronecker(y.cov, rep(1, ncol) )

	if( time.cov > 0 ) 
		for(i in 1:time.cov )
			X <- cbind(X, c(t(y.time))**i)	

	X.i <- lapply(1:nrow, function(i){ 
			YP<-t(y.long[i,,drop=F]); 
			XP<-X[(i*ncol - ncol+1):(i * ncol),,drop=F]; 
			return(XP[!is.na(YP),,drop=F])});
	Y.i <- lapply(1:nrow, function(i){ 
			YP<-t(y.long[i,,drop=F]); 
			return(YP[!is.na(YP),1])});
	U.i <- lapply(1:nrow, function(i){ 
			YP <- t(y.long[i,,drop=F]); 
			UD.i <- diag(1, ncol)[!is.na(YP), !is.na(YP), drop=F];
			return(as.matrix( cbind(1, UD.i) ) );
			});
	D.i <- lapply(1:nrow, function(i){
			return(diag(1,nrow(U.i[[i]])));} );

#show(head(y.long));
#show(head(y.cov));
#show(head(y.time));
	
	RMEL<-function(sig.a, sig.b, sig.e, par_rho)
	{
		AR.1 <- array(1:ncol, dim=c(ncol,ncol))
		AR.1 <- par_rho^abs(AR.1-t(AR.1));
	
		D.i <- lapply(1:nrow, function(i){
					D <- diag(1,ncol(U.i[[i]]));
					YP <- t(y.long[i,,drop=F]); 
					AR.i <- AR.1[!is.na(YP), !is.na(YP), drop=F];
					D[c(2:nrow(D)), c(2:nrow(D))] <- sig.b^2*AR.i; 
					D[1,1] <- sig.a^2;
					return(D);} );

		E.i  <- diag(sig.e^2, ncol);

		V_1 <- lapply(1:nrow, function(k) { kcol <- NROW(U.i[[k]]); solve( U.i[[k]] %*% D.i[[k]] %*% t(U.i[[k]]) + E.i[1:kcol, 1:kcol] ) } );
		XVX.i <- lapply(1:nrow, function(k){ t(X.i[[k]]) %*% V_1[[k]] %*% X.i[[k]]; });
		XVX <- 0;
		for(k in 1:nrow) XVX <- XVX + XVX.i[[k]];
		
		XVY.i <- lapply(1:nrow, function(k){ t(X.i[[k]]) %*% V_1[[k]] %*% Y.i[[k]]; });
		XVY <- 0;
		for(k in 1:nrow) XVY <- XVY + XVY.i[[k]];
		
		B <- solve(XVX) %*% XVY;
		YXB <- lapply(1:nrow, function(k){ Y.i[[k]] - X.i[[k]]%*%B; });

		Y.re <-  lapply(1:nrow, function(k){ t(YXB[[k]]) %*% V_1[[k]] %*% YXB[[k]]; }) ;		
		
		log.V.det <- lapply(1:nrow, function(k){ log(1/det(V_1[[k]]))});
		LR <-  -0.5*(sum( unlist(log.V.det )) + sum(unlist(Y.re))) - 0.5*log(abs(det(XVX)))
		
		return(list(LR=LR, B=B, YXB=YXB));
	}
		
	get_par<-function(par)
	{	
		par_rho<- par[1];
		if ( par_rho<0 || par_rho>0.99 )
			return(NaN);

		sig.a <- abs(par[2]);
		sig.b <- abs(par[3]);
		sig.e <- abs(par[4]);

		r <- RMEL(sig.a, sig.b, sig.e, par_rho);
		
		if(is.infinite(r$LR))
			r$LR <- .Machine$double.xmax

#cat( -r$LR, par, "\n");

		return( -r$LR );
	}

	# par[1] = rho
	# par[2] = sig.a
	# par[3] = sig.b
	# par[4] = sig.e
	# e.g.
	# par.init <- c(rho=0.75, sig.a=0.2, sig.b=0.3, sig.e=0.1 ); 
	
	if (is.null(par.init) || length(par.init)==0)
		par.init <- c( 0.5, rep(sd(y.long, na.rm=T)/3, 3) );

cat("REML parin=", par.init, "\n");

	sd.ref <- sd(y.long, na.rm=T);
	tolerance <- 1;
	loop.n <- 0;
	min.val <- Inf;
	min.par <- par.init;

	while( loop.n <= g.maxiter && tolerance > 1e-5 )
	{
		#r0 <- try( optim( par.init, get_par, method = "L-BFGS-B", #control=list(maxit=2500),
		#		lower=c(0.01, rep(sd.ref/100, 3 )) ,upper=c(0.99, rep(sd.ref*10, 3) ) ), silent = F );
		r0 <- try( optim( par.init, get_par, method = "BFGS", control=list(maxit=5000) ), silent = F );

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
			if(verbose) cat("  LOOP =", loop.n, "/", g.maxiter, " val=", r0$val,  "par=", r0$par, "\n");
			
			tolerance <- max( c( min.val, par.init) - c(r0$value, r0$par) );
			
			min.val <- r0$value;
			min.par <- r0$par;
			par.init<- min.par;
		}

		par.init <- par.init*runif( length(par.init), 0.8, 1.2 );
	}

	if(verbose) cat("  Final =", loop.n, " min.val=", min.val, "par=", min.par, "\n");
	
	if( is.infinite( min.val) || is.na(min.val) || any(is.na(min.par))  )
		return( list(bSuccess=F) );

	par_rho<- min.par[1]
	sig.a  <- abs(min.par[2]);
	sig.b  <- abs(min.par[3]);
	sig.e  <- abs(min.par[4]);

	LR <- RMEL(sig.a, sig.b, sig.e, par_rho  );
	
cat("SIG_A/B/E/R=", min.par, "COV=", c(LR$B), "DELT=", range(LR$YXB), "\n"); 

	cov.effect <- c(LR$B);
	par_mu  <- NA;
	time.effect   <- NA;
	
	if(intercept)
	{
		par_mu  <- cov.effect[1];
		cov.effect <- cov.effect[-1];
	}
	
	if( time.cov > 0)
	{
		time.effect <- cov.effect[-c(1:NCOL(y.cov))];
		cov.effect <- cov.effect[ c(1:NCOL(y.cov))];
	}


	y.delt <- y.long - t( array( X %*% LR$B, dim=c(NCOL(y.long), NROW(y.long))))

	pars <- list( intercept=intercept,
				  time.cov  = time.cov,
			  	  mu     = par_mu, 
			  	  rho    = par_rho, 
			  	  sig.a  = sig.a, 
			  	  sig.b  = sig.b, 
			  	  sig.e  = sig.e, 
			  	  cov.effect= cov.effect,  
			  	  time.effect  = time.effect);  

	r.model <- list(par = pars, likelihood = min.val, 
			phe.delt = y.delt, 
			phe.time = y.time, 
			phe.cov = y.cov, 
			bSuccess=T );
	
	class(r.model) <- "LSKAT.null.model";

	return(r.model);
}

#public
print.LSKAT.null.model <- function(r.model, useS4=FALSE)
{
	cat("  LSKAT/Gene Summary\n");	
	cat("[1] MLE Results:\n");	
	cat("* Intercept =",       r.model$par$intercept, "\n");
	cat("* SIGMA_A =", 		   r.model$par$sig.a, "\n");
	cat("* SIGMA_B =",         r.model$par$sig.b, "\n");
	cat("* SIGMA_E =",         r.model$par$sig.e, "\n");
	cat("* RHO =",             r.model$par$rho, "\n");
	cat("* MU =",              r.model$par$mu, "\n");
	cat("* Beta(Cov)=",        r.model$par$cov.effect, "\n");
	cat("* Beta(Time)=",       r.model$par$time.effect, "\n");
	cat("* L(min) =",          r.model$likelihood, "\n");
}


#private
get_ylog_list<-function(y.long)
{
	y.long.list <-list(); 
	y.una <-apply(y.long, 1, function(y.i){ length( which(!is.na(y.i))); } )

	for(i in 1:NCOL(y.long))
		y.long.list[[i]] <- y.long[which(y.una==i),c(1:i),drop=F]
	
	return(y.long.list)
}

longskat_est_ML<-function( y.long, y.cov, y.time = NULL, time.cov=0, intercept=TRUE, g.maxiter=20, init.par=list(), verbose=F )
{
	# check id matched before call here!
	
	if(is.data.frame(y.long))  y.long <- data.matrix(y.long)
	if(is.data.frame(y.cov))  y.cov <- data.matrix(y.cov)
	if(is.data.frame(y.time)) y.time<- data.matrix(y.time)

	ncol <- NCOL(y.long);
	nrow <- NROW(y.long);
	nCov <- NCOL(y.cov);

	if( is.null(y.time))
		y.time <- t(matrix(rep(c(1:ncol),nrow), nrow=ncol)) ;
	
	if(intercept) 
		X <-  kronecker(cbind(1, y.cov), rep(1, ncol) )
	else
		X <-  kronecker(y.cov, rep(1, ncol) )

	if( time.cov > 0 ) 
		for(i in 1:time.cov)
			X <- cbind(X, c(t(y.time))**i)	

#show(head(y.long));
#show(head(y.cov));
#show(head(y.time));

	get_par<-function(par, y, y.time, y.cov )
	{	
		par_rho<-par[1];
		if ( par_rho<0 || par_rho>=0.99 )
			return(NaN);

		sig.a<-par[2];
		sig.b<-par[3];
		sig.e<-par[4];

		AR.1 <- array(1:ncol, dim=c(ncol,ncol))
		AR.1 <- par_rho^abs(AR.1-t(AR.1));
		sigma <- diag(sig.e^2, ncol) + sig.a^2 + sig.b^2*AR.1;

		cov.effect <- par[ c(5:(4+NCOL(X))) ];
		y.delt <- y.long - t( array( X %*% cov.effect, dim=c(NCOL(y.long), NROW(y.long))))

		A <- 0;
		if(any(is.na(y.delt)))
		{
			for(i in 1:NROW(y.delt) )
			{
				t.sel <- which( !is.na(c(y.delt[i,]) ) & !is.na(c(y.time[i,]) ) );
				sig <- sigma[t.sel, t.sel, drop=F];
				A <- A + sum( dmvnorm( y.delt[i,t.sel,drop=F], rep(0, length(t.sel)), sig, log=T ) );
			}
		}
		else
			A <- sum( dmvnorm( y.delt,  rep(0, ncol), sigma, log=T  ) );

		return( -A );
	}
	
	est_par_cov<-function()
	{
		#1: rho, 2:sig.a, 3:sig.b, 4:sig.e
		par.init <- c( 0.5, sd(y.long, na.rm=T)/3, sd(y.long, na.rm=T)/3, sd(y.long, na.rm=T)/3 );

		if(intercept)
			par.init <- c(par.init, mean(y.long, na.rm=T));

		for(i in 1:nCov)
			par.init <- c( par.init, 1/( mean(y.cov[,i], na.rm=T)^2 + 1) );

		if (time.cov>0)
		{
			par.init <- c( par.init, rep( 1/( mean(y.time, na.rm=T)^2 +1) , time.cov ) )
			if(is.null(y.time))
				y.time <- t( t(ifelse(is.na(y.long),NA, 1))*(rep(1:NROW(y.long))))
		}

		return(par.init);			
	}

	# par[1] = rho
	# par[2] = sig.a
	# par[3] = sig.b
	# par[4] = sig.e
	# par[5] = par_u
	# par[6, 5+nCov] = cov.effect
	# par[6+nCov, (6:7+nCov) ] = time.effect
	# e.g.
	# par.init <- c(rho=0.75, sig.a=0.2, sig.b=0.3, sig.e=0.1, u=1, a=0.5, b=0.5); 
	
	tolerance <- 1;
	loop.n    <- 0;
	min.val   <- Inf;
	par.init  <- min.par <- est_par_cov();
	
	if(verbose)
		cat("  Initial parameters:", par.init, "\n");
	
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
			if(verbose) cat("  LOOP =", loop.n, "/", g.maxiter, " val=", r0$val,  "par=", r0$par, "\n");
			
			tolerance <- max( c( min.val, par.init) - c(r0$value, r0$par) );
			
			min.val <- r0$value;
			min.par <- r0$par;
			par.init<- min.par;
		}

		par.init <- par.init*runif( length(par.init), 0.8, 1.2 );
	}

	if(verbose) cat("  Final =", loop.n, " min.val=", min.val, "par=", min.par, "\n");
	
	if( is.infinite( min.val) || is.na(min.val) || any(is.na(min.par))  )
		return( list(bSuccess=F) );
	

cat("SIG_A/B/E/R=", min.par[c(1:4)], "COV=", min.par[-c(1:4)], "LIKELIHOOD=", min.val, "\n"); 

	par_rho <- min.par[1];
	sig.a   <- abs(min.par[2]);
	sig.b   <- abs(min.par[3]);
	sig.e   <- abs(min.par[4]);
	cov.effect <- min.par[-c(1:4)];
	par_mu  <- NA;
	time.effect   <- NA;
	
	if(intercept)
	{
		par_mu  <- cov.effect[1];
		cov.effect <- cov.effect[-1];
	}
	
	if( time.cov > 0)
	{
		time.effect <- cov.effect[-c(1:NCOL(y.cov))];
		cov.effect <- cov.effect[ c(1:NCOL(y.cov))];
	}

	y.delt <- y.long - t( array( X %*% min.par[-c(1:4)], dim=c(NCOL(y.long), NROW(y.long))))

	pars <- list( intercept=intercept,
				  time.cov= time.cov,
			  	  mu      = par_mu, 
			  	  rho     = par_rho, 
			  	  sig.a   = sig.a, 
			  	  sig.b   = sig.b, 
			  	  sig.e   = sig.e, 
			  	  cov.effect  = cov.effect,  
			  	  time.effect = time.effect );  

	r.model <- list(par = pars, likelihood = min.val, 
			phe.delt = y.delt, 
			phe.time = y.time, 
			phe.cov = y.cov, 
			bSuccess=T );
	
	class(r.model) <- "LSKAT.null.model";
	return(r.model);
}

