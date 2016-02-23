/* lskat_R.cpp  -	LSKAT GENE application
 *	Copyright (C) 2014
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "lskat_R.h"

#include "fm_rlogger.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_err.h"
#include "fm_new.h"

/*
est_gen_Q.R<-function(y.delt, maf, Z, X, par_null)
{
	sig_a2  <- par_null[1]^2
	sig_b2  <- par_null[2]^2
	sig_e2  <- par_null[3]^2
	par_rho <- par_null[4]

	n <- dim(y.delt)[1];
	m <- dim(y.delt)[2];
	k <- length(maf);

	Q.v<-0;
	Q.w<-c();

	AR.1 <- array(0,dim=c(m,m));
	for(i in 1:m)
	for(j in 1:m)
		AR.1[i,j] <- par_rho^abs(i-j);

	V.j <- array(1, dim=c(m,m)) * sig_a2 +  AR.1 * sig_b2 + diag(1, m) * sig_e2
	V.j_1 <- solve(V.j);

	V.j_x <- list();
	Y.j_x <- list();
	M.j_x <- c();

	for(i in 1:n)
	{
		y.i <- y.delt[i,];
		y.j_i <- V.j[!is.na(y.i), !is.na(y.i)];
		V.j_x[[i]] <- solve(y.j_i);
		Y.j_x[[i]] <- y.i[!is.na(y.i)];
		M.j_x[i] <- length(Y.j_x[[i]]);
	}

	for(i in 1:k)
	{
		Q.i <- 0;
		for(j in 1:n)
			Q.i <- Q.i + Z[j,i]*t(rep(1, M.j_x[j] )) %*% V.j_x[[j]] %*% (Y.j_x[[j]])
		Q.v <- Q.v + (Q.i)^2;
	}

	Q.v <- Q.v/2

	n.x <- dim(X)[2];
	W0 <- array(0, dim=c(k,k));
	W1 <- array(0, dim=c(k,n.x));
	W2 <- array(0, dim=c(n.x,n.x));
	W3 <- array(0, dim=c(n.x,k));
	for(i in 1:n)
	{
		m.i <- M.j_x[i]
		kr.Z<- kronecker( Z[i,,drop=F], array(1, dim=c(m.i,1)) )
		kr.X<- kronecker( X[i,,drop=F], array(1, dim=c(m.i,1)) )

		W0 <- W0 + t(kr.Z) %*% V.j_x[[i]] %*% kr.Z;
		W1 <- W1 + t(kr.Z) %*% V.j_x[[i]] %*% kr.X;
		W2 <- W2 + t(kr.X) %*% V.j_x[[i]] %*% kr.X;
		W3 <- W3 + t(kr.X) %*% V.j_x[[i]] %*% kr.Z;
	}

	Q.w <- W0 - W1 %*% solve(W2) %*% W3;

	return(list(v=Q.v, w=Q.w/2));
}

*/

int kronecker_vm( CFmVector& A, CFmMatrix& B, CFmMatrix* pRet )
{
	int NB = B.GetNumCols();
	pRet->Resize( 1 * B.GetNumRows(), A.GetLength() * NB );

	for(int i=0; i<A.GetLength(); i++)
	{
		for(int j=0; j<B.GetNumRows(); j++)
		for(int k=0; k<NB; k++)
			pRet->Set(j, k+i*NB, A[i]*B.Get(j,k) );
	}

	return(0);
}

SEXP Xest_gen_Q_C( CFmMatrix* pFmYDelt, CFmMatrix* pFmZ, CFmMatrix* pFmX, CFmVector* pFmMaf, CFmVector* pFmParNull)
{
	CFmNewTemp fmRef;

	double sig_a2  = pow(pFmParNull->Get(0), 2);
	double sig_b2  = pow(pFmParNull->Get(1), 2);
	double sig_e2  = pow(pFmParNull->Get(2), 2);
	double par_rho = pFmParNull->Get(3);

	int N = pFmYDelt->GetNumRows();
	int M = pFmYDelt->GetNumCols();
	int K = pFmMaf->GetLength();

//Rprintf("N=%d, M=%d, K=%d a2=%f b2=%f e2=%f rho=%f\n", N, M, K, sig_a2, sig_b2, sig_e2, par_rho);

	CFmMatrix fmAR1( M, M );
	for(int i=0; i<M; i++)
	for(int j=0; j<M; j++)
		fmAR1.Set( i, j, pow( par_rho, abs( i - j ) ) );

	CFmMatrix fmV_j( M, M );
	CFmMatrix fmV_j1( M, M );
	//CFmMatrix fmDiag( M, true, 1.0 ); //???HERE
	CFmMatrix fmDiag( M, M );
	for(int i=0; i<M; i++) fmDiag.Set(i, i, 1.0);

	fmV_j = (fmV_j + 1.0) * sig_a2 +  fmAR1 * sig_b2 + fmDiag * sig_e2;
	fmV_j1 = fmV_j.GetInverted();

	CFmVector fmVectMj_x(N, 0.0);
	CFmVector fmVecTmp (M, 0.0);
	CFmVector fmVecTmp2(M, 0.0);
	CFmMatrix fmVj_i(M, M);

	CFmMatrix** ppVj = Calloc(N, CFmMatrix*);
	CFmVector** ppYj = Calloc( N, CFmVector*);

	for(int i=0; i<N ;i++)
	{
		fmVecTmp = pFmYDelt->GetRow(i);
		fmVecTmp2.Resize(0);
		for(int j=0; j<fmVecTmp.GetLength(); j++)
		{
			if (!isnan(fmVecTmp[j]))
				fmVecTmp2.Put(j);
		}

		int NonNA = fmVecTmp2.GetLength();

		ppVj[i] = new (fmRef) CFmMatrix(NonNA, NonNA);
		ppYj[i] = new (fmRef) CFmVector(NonNA, 0.0);

		fmVj_i.Resize(NonNA, NonNA);
		if (NonNA>0)
		{
			for( int k=0; k<NonNA; k++)
			for( int l=0; l<NonNA; l++)
				fmVj_i.Set(k, l,fmV_j.Get( (int)fmVecTmp2[k], (int)fmVecTmp2[l] ) );
			*(ppVj[i]) = fmVj_i.GetInverted( );

			fmVecTmp.RemoveNan();
			*(ppYj[i]) = fmVecTmp;

			fmVectMj_x[i] = fmVecTmp.GetLength();
		}

	}

	CFmVector fmQi(1, 0.0);
	CFmMatrix fmTrans(1, N );
	double fQv=0.0;

	for(int i=0; i<K ;i++)
	{
		fmQi.Resize(1);
		for(int j=0; j<N ;j++)
		{
			fmTrans.Resize(1, fmVectMj_x[j]);
			for( int l=0;l<(int)(fmVectMj_x[j]);l++)
				fmTrans.Set(0, l, 1.0);
			fmQi = (fmTrans * (*(ppVj[j])) * (*(ppYj[j])) * pFmZ->Get(j,i)).GetRow(0) + fmQi;
		}

		fmQi = fmQi*fmQi;
		fQv = fQv + fmQi.Sum();
	}

	fQv = fQv/2.0;

	int NX = pFmX->GetNumCols();
	CFmMatrix fmW0( K, K );
	CFmMatrix fmW1( K, NX );
	CFmMatrix fmW2( NX,NX );
	CFmMatrix fmW3( NX,K );
	CFmMatrix fmKrZ( 0,0);
	CFmMatrix fmKrX( 0,0);
	CFmMatrix fmKron( N, 1 );

	for(int i=0; i<N; i++)
	{
		fmKron.Resize( fmVectMj_x[i],1 );
		for(int k=0;k<fmVectMj_x[i]; k++) fmKron.Set(k, 0, 1.0);

		fmVecTmp = pFmZ->GetRow(i);
		kronecker_vm( fmVecTmp, fmKron, &fmKrZ );
		fmVecTmp = pFmX->GetRow(i);
		kronecker_vm( fmVecTmp, fmKron, &fmKrX );

		fmW0 = fmW0 + fmKrZ.GetTransposed() * (*(ppVj[i])) * fmKrZ;
		fmW1 = fmW1 + fmKrZ.GetTransposed() * (*(ppVj[i])) * fmKrX;
		fmW2 = fmW2 + fmKrX.GetTransposed() * (*(ppVj[i])) * fmKrX;
		fmW3 = fmW3 + fmKrX.GetTransposed() * (*(ppVj[i])) * fmKrZ;
	}

	CFmMatrix fmQw( K, K );

	fmQw = fmW0 - fmW1 * fmW2.GetInverted() * fmW3;
	fmQw = fmQw / 2.0;

	for(int i=0; i<N; i++) { destroy( ppVj[i] );}
	for(int i=0; i<N; i++) { destroy( ppYj[i] );}
	Free(ppVj);
	Free(ppYj);

	//double fQv = 0.5;
	//CFmMatrix fmQw( K, K );
	//for(int i=0; i<K; i++)  fmQw.Set(i, i, i+1);

	SEXP sRet, t;
   	PROTECT(sRet = t = allocList(2));

	SEXP expVS = GetSEXP(&fmQw);
	SETCAR( t, expVS );
	SET_TAG(t, install("w") );
	t = CDR(t);

	CFmVector frmQv(1, fQv);
	SEXP expVS1 = GetSEXP(&frmQv);
	SETCAR( t, expVS1 );
	SET_TAG(t, install("v") );
	t = CDR(t);

	UNPROTECT(1);

    return(sRet);
}

CFmMatrix* getMatrixData(SEXP pMat)
{
	SEXP Rdim = getAttrib(pMat, R_DimSymbol);
	int nrow = INTEGER(Rdim)[0];
	int ncol = INTEGER(Rdim)[1];

    double* p0 = REAL(pMat) ;

 	CFmNewTemp fmRef;
 	CFmMatrix* p = new (fmRef) CFmMatrix(nrow, ncol) ;

    int i,j;
    for( i=0; i<nrow; i++)
    	for( j=0; j<ncol; j++)
               p->Set(i, j, p0[i+nrow*j] ) ;

	return(p);
}

CFmVector* getVectorData(SEXP pVec)
{
	int nlen = length(pVec);

    double* p0 = REAL(pVec) ;

 	CFmNewTemp fmRef;
    CFmVector* p = new (fmRef) CFmVector( nlen, 0.0) ;

    int i,j;
    for( i=0; i<nlen; i++)
               p->Set(i,  p0[i]);

	return(p);
}

SEXP _est_gen_Q_C( SEXP spYdelt,
  		   	SEXP spZ,
  		   	SEXP spX,
  		   	SEXP spMaf,
		   	SEXP spParNull)
{
	// int nUsed0, nTotal0;
	// CFmVector::StatCache( &nTotal0, &nUsed0 );
	// int nUsed1, nTotal1;
	// CFmMatrix::StatCache( &nTotal1, &nUsed1 );
	// Rprintf( "Enter C Range, Vec.count=%d, Mat.count=%d\n", nTotal0, nTotal1);

	CFmMatrix* pFmYDelt = getMatrixData(spYdelt);
	CFmMatrix* pFmZ     = getMatrixData(spZ);
	CFmMatrix* pFmX     = getMatrixData(spX);
	CFmVector* pFmMaf   = getVectorData(spMaf);
	CFmVector* pFmParNull = getVectorData(spParNull);

	SEXP ret;
	try
	{
		ret = Xest_gen_Q_C( pFmYDelt, pFmZ, pFmX, pFmMaf, pFmParNull);
	}
    catch(const char* str)
    {
        _log_error( _HI_, "Exception=%s", str);
        return( R_NilValue );
    }

	destroy( pFmYDelt );
	destroy( pFmZ );
	destroy( pFmX );
	destroy( pFmMaf );
	destroy( pFmParNull );

	// CFmVector::StatCache( &nTotal0, &nUsed0 );
	// CFmMatrix::StatCache( &nTotal1, &nUsed1 );
	// Rprintf( "Leave C Range, Vec.count=%d, Mat.count=%d\n", nTotal0, nTotal1);

	return(ret);
}
