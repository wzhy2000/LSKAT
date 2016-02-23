/* bls_RI.cpp  -	BLS application
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "lskat_R.h"

#define BOOLEAN_ELT(x,__i__)	LOGICAL(x)[__i__]
#define INTEGER_ELT(x,__i__)	INTEGER(x)[__i__]
#define NUMERIC_ELT(x,__i__)	REAL(x)[__i__]

SEXP est_gen_Q_C( SEXP spYdelt,
  		   	SEXP spZ,
  		   	SEXP spX,
  		   	SEXP spMaf,
		   	SEXP spParNull)
{
	SEXP ret = _est_gen_Q_C( spYdelt, spZ, spX, spMaf, spParNull);

	return(ret);
}
