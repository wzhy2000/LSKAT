library(mvtnorm)

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

	if(method=="REML")
		ret <- longskat_est_REML( y.long, y.cov, y.time, y.cov.time, intercept, g.maxiter, par.init, debug )
	else	
		ret <- longskat_est_ML( y.long, y.cov, y.time, y.cov.time, intercept, g.maxiter, par.init, debug );
		
	return(ret);	
}

# y.long format:  <shareid> trait1, ..., traitN
# y.time format: <shareid> time1, ..., timeN
# y.cov format:  <shareid> cov1, ..., covM

longskat_est_REML<-function( y.long, y.cov, y.time = NULL, y.cov.time=0, intercept=FALSE, g.maxiter=20, par.init=list(), debug=F )
{
	# check id matched before call here!

	if(is.data.frame(y.long))  y.long <- data.matrix(y.long)
	if(is.data.frame(y.cov))   y.cov <- data.matrix(y.cov)
	if(is.data.frame(y.time))  y.time<- data.matrix(y.time)
	
	ncol <- NCOL(y.long);
	nrow <- NROW(y.long);
	
	if( is.null(y.time))
		y.time <- t(matrix(rep(c(1:ncol),nrow), nrow=ncol)) ;

if(debug)
{
	cat("y.cov", "ncol=", ncol, "intercept=", intercept,  "\n");
	show(head(y.cov));
}

	if(intercept) 
		X <-  kronecker(cbind(1, y.cov), rep(1, ncol) )
	else
		X <-  kronecker(y.cov, rep(1, ncol) )

	if( y.cov.time > 0 ) 
		for(i in 1:y.cov.time)
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
	
	RMEL<-function(sig_a, sig_b, sig_e, par_rho)
	{
		AR.1 <- array(1:ncol, dim=c(ncol,ncol))
		AR.1 <- par_rho^abs(AR.1-t(AR.1));
	
		D.i <- lapply(1:nrow, function(i){
					D <- diag(1,ncol(U.i[[i]]));
					YP <- t(y.long[i,,drop=F]); 
					AR.i <- AR.1[!is.na(YP), !is.na(YP), drop=F];
					D[c(2:nrow(D)), c(2:nrow(D))] <- sig_b^2*AR.i; 
					D[1,1] <- sig_a^2;
					return(D);} );

		E.i  <- diag(sig_e^2, ncol);

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

		sig_a <- abs(par[2]);
		sig_b <- abs(par[3]);
		sig_e <- abs(par[4]);

		r <- RMEL(sig_a, sig_b, sig_e, par_rho);
		
		if(is.infinite(r$LR))
			r$LR <- .Machine$double.xmax

#cat( -r$LR, par, "\n");

		return( -r$LR );
	}

	# par[1] = rho
	# par[2] = sig_a
	# par[3] = sig_b
	# par[4] = sig_e
	# e.g.
	# par.init <- c(rho=0.75, sig_a=0.2, sig_b=0.3, sig_e=0.1 ); 
	
	if (is.null(par.init) || length(par.init)==0)
		par.init <- c( 0.5, sd(y.long, na.rm=T), sd(y.long, na.rm=T), sd(y.long, na.rm=T));

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

	par_rho<- min.par[1]
	sig_a  <- abs(min.par[2]);
	sig_b  <- abs(min.par[3]);
	sig_e  <- abs(min.par[4]);

	LR <- RMEL(sig_a, sig_b, sig_e, par_rho  );
	
cat("SIG_A/B/E/R=", min.par, "COV=", c(LR$B), "DELT=", range(LR$YXB), "\n"); 

	par_cov <- c(LR$B);
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


	y.delt <- y.long - t( array( X %*% LR$B, dim=c(NCOL(y.long), NROW(y.long))))

	pars <- list( intercept=intercept,
				  y.cov.time = y.cov.time,
			  	  mu     = par_mu, 
			  	  rho    = par_rho, 
			  	  sig_a  = sig_a, 
			  	  sig_b  = sig_b, 
			  	  sig_e  = sig_e, 
			  	  par_cov= par_cov,  
			  	  par_t  = par_t );

	r.model <- list(par = pars, likelihood = min.val, y.delt=y.delt, y.time = y.time, y.cov = y.cov, bSuccess=T );
	
	class(r.model) <- "LSKAT.null.model";

	return(r.model);
}

#public
print.LSKAT.null.model <- function(r.model, useS4=FALSE)
{
	cat("  LSKAT/Gene Summary\n");	
	cat("[1] MLE Results:\n");	
	cat("* Intercept =",       r.model$par$intercept, "\n");
	cat("* SIGMA_A =", 		   r.model$par$sig_a, "\n");
	cat("* SIGMA_B =",         r.model$par$sig_b, "\n");
	cat("* SIGMA_E =",         r.model$par$sig_e, "\n");
	cat("* RHO =",             r.model$par$rho, "\n");
	cat("* MU =",              r.model$par$mu, "\n");
	cat("* Beta(Cov)=",        r.model$par$par_cov, "\n");
	cat("* Beta(Time)=",       r.model$par$par_t, "\n");
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

longskat_est_ML<-function( y.long, y.cov, y.time = NULL, y.cov.time=0, intercept=TRUE, g.maxiter=20, init.par=list(), debug=F )
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

#show(head(y.long));
#show(head(y.cov));
#show(head(y.time));

	get_par<-function(par, y, y.time, y.cov )
	{	
		par_rho<-par[1]
		if ( par_rho<0 || par_rho>=0.99 )
			return(NaN);

		sig_a<-par[2]
		sig_b<-par[3]
		sig_e<-par[4]
		par_u<-par[5]
		par_cov <- par[ c(6:(5+nCov)) ];
		
		#y.ind <- par[c((6+nCov):(6+nCov+nrow-1))]

		AR.1 <- array(1:ncol, dim=c(ncol,ncol))
		AR.1 <- par_rho^abs(AR.1-t(AR.1));
		sigma <- diag(sig_e^2, ncol) + sig_a^2 + sig_b^2*AR.1;

		y.delt <- y - cbind(1, y.cov)%*%(c(par_u, par_cov))%*%array(1, dim=c(1,ncol))

		par_t   <- c();
		if( y.cov.time >0 )
		{
			y.time.max <- max( y.time, na.rm=T);	
			y.time.min <- min( y.time, na.rm=T);	

			y.time <- (y.time - y.time.min)/(y.time.max - y.time.min)*2 - 1;
			par_t  <-  par[5 + nCov + c(1:y.cov.time)];
			
			if(y.cov.time >= 1)
				y.delt <- y.delt - (y.time) * par_t[1];
			if(y.cov.time >= 2)
				y.delt <- y.delt - 0.5 * (3*y.time^2 - 1) * par_t[2];
			if(y.cov.time >= 3)
				y.delt <- y.delt - 0.5 * (5*y.time^3 - 3*y.time) * par_t[3];
			if(y.cov.time >= 4)
				y.delt <- y.delt - 1/8 * (35*y.time^4 - 30*y.time^2 + 3 ) * par_t[4];
		}
		
		# y.temp <- c(y.delt)
		# if ( length(which(is.na(y.temp)))>0) y.temp[ is.na(y.temp) ]<-0;
		# y.delt <- array(y.temp, dim=c(nrow, ncol));
		# A <- -sum( dmvnorm( y.delt,  rep(0, ncol), sigma, log=T ) ) 
		
		#y.list <- get_ylog_list(y.delt);
		#A <- 0;
		#for(i in 1:length(y.list) )
		#{
		#	if(!is.null(y.list[[i]]) && NROW(y.list[[i]])>0)
		#		A <- A - sum( dmvnorm( y.list[[i]], rep(0, i), sigma[1:i, 1:i, drop=F], log=T ) );

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
		{
			A <- sum( dmvnorm( y.delt,  rep(0, ncol), sigma, log=T  ) );
			
			#y.mean <- y.ind %*% t(rep(1, ncol));
			#A <- sum( dmvnorm( y.delt-y.mean,  rep(0, ncol), sigma, log=T  ) );
		}	
		return( -A );
	}
	
	est_par_cov<-function()
	{
		par.init.cov <- c(mean(y.long, na.rm=T));
		for(i in 1:nCov)
			par.init.cov <- c( par.init.cov, 1/( mean(y.cov[,i], na.rm=T)^2 + 1) );

		return(par.init.cov);			
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

	par.init.cov <- est_par_cov();
	par.init <- c( 0.5, sd(y.long, na.rm=T), sd(y.long, na.rm=T), sd(y.long, na.rm=T), par.init.cov );

	if (y.cov.time>0)
	{
		par.init <- c( par.init, rep( 1/( mean(y.time, na.rm=T)^2 +1) , y.cov.time ) )
		if(is.null(y.time))
			y.time <- t( t(ifelse(is.na(y.long),NA, 1))*(rep(1:NROW(y.long))))
	}
	
	tolerance <- 1;
	loop.n <- 0;
	min.val <- Inf;
	min.par <- par.init;

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

#min.par <- c( 0.9632793, 0.0008507287, 4.260085, 0.8689705, 27.82508, -1.522521, -6.556081, 16.20164, 5.955584, -32.51609, 8.446485, 2.068376, -2.723876, -0.9025116, -0.3033274);
#min.par <- c( 0.7, 0.8, 0.8, 0.8, 1, 0.5, 0.5);
#min.val <- 37391.62;
#loop.n  <- 10.2;
	
	if(debug) cat("  Final =", loop.n, " min.val=", min.val, "par=", min.par, "\n");
	
	if( is.infinite( min.val) || is.na(min.val) || any(is.na(min.par))  )
		return( list(bSuccess=F) );

	par_u   <- min.par[5];
	par_cov <- min.par[c(6:(5+nCov))];

	if(!is.null(init.par$mu)) 	par_u <- init.par$mu;
	if(!is.null(init.par$cov)) 	par_cov <- init.par$cov;

	y.delt  <- y.long - cbind(1, y.cov)%*%(c( par_u, par_cov))%*%array(1, dim=c(1,ncol))

	par_t <- c();	
	if( y.cov.time>0 )
	{
		#par_t  <- min.par[ 5 + nCov +c(1:y.cov.time)];
		#for( i in 1:length(par_t))
		#	y.delt <- y.delt - (y.time^i) * par_t[i];	

		y.time.max <- max( y.time, na.rm=T);	
		y.time.min <- min( y.time, na.rm=T);	

		y.time <- (y.time - y.time.min)/(y.time.max - y.time.min)*2 - 1;
		par_t  <-  min.par[5 + nCov + c(1:y.cov.time)];

		if(y.cov.time >= 1)
			y.delt <- y.delt - (y.time) * par_t[1];
		if(y.cov.time >= 2)
			y.delt <- y.delt - 0.5 * (3*y.time^2 - 1) * par_t[2];
		if(y.cov.time >= 3)
			y.delt <- y.delt - 0.5 * (5*y.time^3 - 3*y.time) * par_t[3];
		if(y.cov.time >= 4)
			y.delt <- y.delt - 1/8 * (35*y.time^4 - 30*y.time^2 + 3 ) * par_t[4];
	}

cat(par_u, par_cov, par_t, "DELT=", range(y.delt), "\n"); 
	
	pars <- list( intercept = TRUE,
				  y.cov.time = y.cov.time,
				  mu     = par_u, 
			  	  rho    = min.par[1], 
			  	  sig_a  = abs(min.par[2]), 
			  	  sig_b  = abs(min.par[3]), 
			  	  sig_e  = abs(min.par[4]), 
		  	      par_cov= par_cov, 
		  	      par_t  = par_t );
	r.model <- list(par = pars, likelihood = min.val, y.delt=y.delt, y.time = y.time, y.cov = y.cov );
	class(r.model) <- "LSKAT.null.model";

	return(r.model);
}

longskat_est_model_LR_random<-function( y.long, y.cov, y.time = NULL, y.cov.time=0, g.maxiter=20, random.effect=NULL, init.par=list(), debug=F )
{
	# check id matched before call here!
	
	if(is.data.frame(y.long))  y.long <- data.matrix(y.long)
	if(is.data.frame(y.cov))  y.cov <- data.matrix(y.cov)
	if(is.data.frame(y.time)) y.time<- data.matrix(y.time)
	
	ncol <- NCOL(y.long);
	nrow <- NROW(y.long);
	nCov <- NCOL(y.cov);

	if(is.null(random.effect))
		u.mat <- array(0, dim=c(nrow, ncol))
	else
		u.mat <- as.matrix(random.effect);
	
	if( is.null(y.time))
		y.time <- t(matrix(rep(c(1:ncol),nrow), nrow=ncol)) ;

	X <-  kronecker(cbind(1, y.cov), rep(1, ncol) )
	if( y.cov.time > 0 ) 
		for(i in 1:y.cov.time)
			X <- cbind(X, t(c(y.time))**i)	

	X.i <- lapply(1:nrow, function(i){ 
			YP<-t(y.long[i,,drop=F]); 
			XP<-X[(i*ncol - ncol+1):(i * ncol),,drop=F]; 
			return(XP[!is.na(YP),,drop=F])});
	Y.i <- lapply(1:nrow, function(i){ 
			YP<-t(y.long[i,,drop=F]); 
			return(YP[!is.na(YP),1])});
	U.i <- lapply(1:nrow, function(i){ 
			YP <- t(y.long[i,,drop=F]); 
			UP <- kronecker(u.mat[i,,drop=F], rep(1, ncol) )
			return(UP[!is.na(YP),,drop=F])});
	D.i <- lapply(1:nrow, function(i){ 
			D <- array(1, dim=c(ncol(u.mat),ncol(u.mat)));
			return(D);});
			
	RMEL<-function(sig_a, sig_b, sig_e, par_rho)
	{
		AR.1 <- array(1:ncol, dim=c(ncol,ncol))
		AR.1 <- par_rho^abs(AR.1-t(AR.1));
		E.i  <- diag(sig_e^2, ncol) + sig_a^2 + sig_b^2*AR.1;
		
		V_1 <- lapply(1:nrow, function(k) { solve( U.i[[k]] %*% D.i[[k]] %*% t(U.i[[k]]) + E.i )} );
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
		sig_a <- abs(par[2]);
		sig_b <- abs(par[3]);
		sig_e <- abs(par[4]);
		par_rho<- par[1];
		if ( par_rho<0 || par_rho>=0.99 )
			return(NaN);

		r <- RMEL(sig_a, sig_b, sig_e, par_rho);
		
		return( -r$LR );
	}

	# par[1] = rho
	# par[2] = sig_a
	# par[3] = sig_b
	# par[4] = sig_e
	# e.g.
	# par.init <- c(rho=0.75, sig_a=0.2, sig_b=0.3, sig_e=0.1 ); 

	par.init <- c( 0.5, sd(y.long, na.rm=T), sd(y.long, na.rm=T), sd(y.long, na.rm=T) );

	tolerance <- 1;
	loop.n <- 0;
	min.val <- Inf;
	min.par <- par.init;

	while( loop.n <= g.maxiter && tolerance > 1e-5 )
	{
		r0 <- try( optim( par.init, get_par, method = "BFGS", control=list(maxit=500) ), silent = F );

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

	par_rho<- min.par[1]
	sig_a  <- abs(min.par[2]);
	sig_b  <- abs(min.par[3]);
	sig_e  <- abs(min.par[4]);

	LR <- RMEL(sig_a, sig_b, sig_e, par_rho);
	
cat("COV=", c(LR$B), "DELT=", range(LR$YXB), "\n"); 

	par_cov <- c(LR$B);
	par_mu  <- par_cov[1];
	par_cov <- par_cov[-1];
	par_t <- par_cov[-c(1:nCov)];

	y.delt <- y.long - cbind(1,y.cov) %*% LR$B %*% rep(1, ncol);
	
	pars <- list( mu     = par_mu, 
			  	  rho    = min.par[1], 
			  	  sig_a  = sig_a, 
			  	  sig_b  = sig_b, 
			  	  sig_e  = sig_e, 
		  	      par_cov= par_cov, 
		  	      par_t  = par_t );

	r.model <- list(par = pars, likelihood = min.val, y.cov.time=y.cov.time, y.delt=y.delt, y.time = y.time, y.cov = y.cov );
	
	class(r.model) <- "LSKAT.null.model";

	return(r.model);
}
