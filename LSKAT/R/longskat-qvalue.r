
#########################################################
#
#	Functions from CompQuadForm package
#			date: 12/05/2011

SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001)
{
	r <- length(lambda)
	if (length(h) != r) stop("lambda and h should have the same length!")
	if (length(delta) != r) stop("lambda and delta should have the same length!")

	out <- .C("qfc",lambdas = as.double(lambda),
				noncentral  = as.double(delta),
				df          = as.integer(h),
				r           = as.integer(r),
				sigma       = as.double(sigma),
				q           = as.double(q),
				lim         = as.integer(lim),
				acc         = as.double(acc),
				trace       = as.double(rep(0,7)),
				ifault      = as.integer(0),
				res         = as.double(0),PACKAGE="LSKAT")

	out$res <- 1 - out$res

	return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
}

SKAT_liu <- function(q, lambda, h = rep(1,length(lambda)), delta = rep(0,length(lambda)))
{
	r <- length(lambda)
	if (length(h) != r) stop("lambda and h should have the same length!")
	if (length(delta) != r) stop("lambda and delta should have the same length!")

	c1 <- sum(lambda*h) + sum(lambda*delta)

	c2 <- sum(lambda^2*h) + 2*sum(lambda^2*delta)

	c3 <- sum(lambda^3*h) + 3*sum(lambda^3*delta)

	c4 <- sum(lambda^4*h) + 4*sum(lambda^4*delta)

	s1 <- c3/(c2^(3/2))

	s2 <- c4/c2^2

	muQ <- c1

	sigmaQ <- sqrt(2*c2)

	tstar <- (q-muQ)/sigmaQ

	if (s1^2>s2) {
		a <- 1/(s1-sqrt(s1^2-s2))
		delta <- s1*a^3-a^2
		l <- a^2-2*delta
	} else {
		a <- 1/s1
		delta <- 0
		l <- c2^3/c3^2
	}

	muX <- l+delta
	sigmaX <- sqrt(2)*a
	Qq <- pchisq(tstar*sigmaX+muX,df=l,ncp=delta,lower.tail=FALSE)
	return(Qq)
}


SKAT_liu.MOD <- function(q, lambda, h = rep(1,length(lambda)), delta = rep(0,length(lambda)))
{
  	r <- length(lambda)
  	if (length(h) != r) stop("lambda and h should have the same length!")
  	if (length(delta) != r) stop("lambda and delta should have the same length!")

  	c1 <- sum(lambda*h) + sum(lambda*delta)
  	c2 <- sum(lambda^2*h) + 2*sum(lambda^2*delta)
  	c3 <- sum(lambda^3*h) + 3*sum(lambda^3*delta)
  	c4 <- sum(lambda^4*h) + 4*sum(lambda^4*delta)

  	s1 <- c3/(c2^(3/2))
  	s2 <- c4/c2^2

  	muQ <- c1
  	sigmaQ <- sqrt(2*c2)
  	tstar <- (q-muQ)/sigmaQ

  	if (s1^2>s2) {
		a <- 1/(s1-sqrt(s1^2-s2))
	   	delta <- s1*a^3-a^2
		l <- a^2-2*delta
	} else {
		delta <- 0
		l = 1/s2
		a = sqrt(l)
	}

	muX <- l+delta
	sigmaX <- sqrt(2)*a
	Qq <- pchisq(tstar*sigmaX+muX,df=l,ncp=delta,lower.tail=FALSE)

	return(Qq)
}

Get_PValue.Lambda<-function(lambda,Q)
{
	#print(lambda)
	n1<-length(Q)

	p.val<-rep(0,n1)
	p.val.liu<-rep(0,n1)
	is_converge<-rep(0,n1)
	p.val.liu<-SKAT_liu.MOD(Q, lambda)

	for(i in 1:n1){
		out<-SKAT_davies(Q[i],lambda,acc=10^(-6))

		p.val[i]<-out$Qq
		is_converge[i]<-1

		# check convergence
		if(length(lambda) == 1){
			p.val[i]<-p.val.liu[i]
		} else if(out$ifault != 0){
			is_converge[i]<-0
		}

		# check p-value
		if(p.val[i] > 1 || p.val[i] <= 0 ){
			is_converge[i]<-0
			p.val[i]<-p.val.liu[i]
		}
	}

	p.val.msg = NULL
	#cat(p.val[1])
	if(p.val[1] == 0)
	{
		#param<-Get_Liu_Params_Mod_Lambda(lambda)
		#p.val.msg<-Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ, param$muX, param$sigmaQ, param$sigmaX, param$l, param$d)
		p.val.msg<-"p.val is close to ZERO."

	}

	return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge, pval.zero.msg=p.val.msg))
}

Get_Lambda<-function(K)
{
	out.s <- try(eigen(K,symmetric=TRUE, only.values = TRUE))
	#print(out.s$values)

	if (class(out.s)=="try-error") { show(K); browser(); }

	#out.s1<-eigen(K,symmetric=TRUE)
	#print(out.s1$values)

	lambda1<-out.s$values
	IDX1<-which(lambda1 >= 0)

	# eigenvalue bigger than sum(eigenvalues)/1000
	IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)

	if(length(IDX2) == 0){
		stop("No Eigenvalue is bigger than 0!!")
	}

	lambda<-lambda1[IDX2]
	return(lambda)
}

get_Qv_pvalue<-function(Q, W)
{
	lambda<-Get_Lambda(W)
	re<-Get_PValue.Lambda(lambda,Q)
	return(re)
}

get_Qu_pvalue<-function(Q, W)
{
	lambda<-Get_Lambda( matrix(1, nrow=NROW(W), ncol=NCOL(W) ) %*% W )
	re<-Get_PValue.Lambda(lambda,Q)
	return(re)
}

SKAT_Scale_Genotypes <- function(X1, Z, weights.common=c(1,1), weights.rare=c(1,25), weights=NULL, rare.cutoff=NULL, test.type="Joint", r.corr.common=0, r.corr.rare=0)
{
	n <- NROW(Z);
	if ( is.null(rare.cutoff) )
		rare.cutoff <- 1/sqrt(2*n);

	Z.maf <- colMeans(Z)/2;

	## NO common SNP, but Common.Only
	if ( length ( which(Z.maf > rare.cutoff) )==0 && toupper(test.type) == toupper("Common.Only"))
		return(list(new=NULL, maf=NULL, rare=0))

	## NO rare SNP, but Rare.Only
	if ( length ( which(Z.maf <= rare.cutoff) )==0 && toupper(test.type) == toupper("Rare.Only"))
		return(list(new=NULL, maf=NULL, rare=0))

	## only have common SNPs
	if ( length ( which(Z.maf <= rare.cutoff) )==0 || toupper(test.type) == toupper("Common.Only"))
	{
		Z <- Z[, Z.maf > rare.cutoff, drop=F];
		Z.maf <- colMeans(Z)/2;

		wr <- get_weights(Z.maf, NROW(Z), weights.common, weights.rare, rare.cutoff);
		Z <- t( t(Z) * wr )
		return(list(new=Z, maf=Z.maf, rare=0))
	}


	## only have RARE SNPs
	if ( length ( which(Z.maf > rare.cutoff) )==0 || toupper(test.type) == toupper("Rare.Only"))
	{
		Z <- Z[, Z.maf <= rare.cutoff, drop=F];
		Z.maf <- colMeans(Z)/2;

		wr <- get_weights(Z.maf, NROW(Z), weights.common, weights.rare,rare.cutoff);
		Z <- t(t(Z) * wr );
		return(list(new=Z, maf=Z.maf, rare=NCOL(Z)))
	}

	Z1 <- Z[ ,which(Z.maf <= rare.cutoff),drop=F]
	Z2 <- Z[ ,which(Z.maf > rare.cutoff),drop=F]

	pi_1 = NULL
	p.m1 <- NCOL(Z1)
	p.m2 <- NCOL(Z2)


	MAF1<- colMeans(Z1)/2
	if(is.null(weights))
	{
		weights<-SKAT:::Beta.Weights(MAF1,weights.rare)
	}

	# Weights for rare variants, but no weights for common variants
	Z1 = t(t(Z1) * (weights))

	MAF2 <- colMeans(Z2)/2
	weights2<-SKAT:::Beta.Weights(MAF2,weights.common)
	Z2 = t(t(Z2) * (weights2))


   	# r.corr
   	if(r.corr.rare == 1){
  		Z1<-cbind(rowSums(Z1))
   	} else if(r.corr.rare > 0){

   		p.m<-dim(Z1)[2]
		R.M<-diag(rep(1-r.corr.rare,p.m)) + matrix(rep( r.corr.rare, p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z1<- Z1 %*% t(L)
   	}

   	if(r.corr.common == 1){
  		Z2<-cbind(rowSums(Z2))
   	} else if(r.corr.common > 0){

   		p.m<-dim(Z2)[2]
		R.M<-diag(rep(1-r.corr.common,p.m)) + matrix(rep(r.corr.common,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z2<- Z2 %*% t(L)
   	}

	Z1.1<-Z1 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z1)
	Z2.1<-Z2 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z2)

	temp1<-t(Z1.1) %*% Z1.1
	temp2<-t(Z2.1) %*% Z2.1
	z1.mean<-sum(diag(temp1))
	z2.mean<-sum(diag(temp2))

	z1.var<-sum(temp1 * temp1)
	z2.var<-sum(temp2 * temp2)


	# sqrt-sqrt because temp1  is ^2
	Z1.new<-Z1/sqrt(sqrt(z1.var))
	Z2.new<-Z2/sqrt(sqrt(z2.var))

	return(list(new=cbind(Z1.new, Z2.new),rare=length(MAF1), maf=c(MAF1, MAF2), Z1=Z1.new, Z2=Z2.new, z1.mean=z1.mean, z2.mean=z2.mean, z1.var=z1.var, z2.var=z2.var, Z1.org=Z1, Z2.org=Z2))
}
