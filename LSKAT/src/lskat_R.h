// lskat_R.h: interface for the LSKAT.
//
//////////////////////////////////////////////////////////////////////

#if !defined(LSKAT_R_H__INCLUDED_)
#define LSKAT_R_H__INCLUDED_

#include <stdbool.h>
#include "fm_linux.h"

#ifdef __cplusplus
extern "C" {
#endif

SEXP _est_gen_Q_C( SEXP spYdelt, SEXP spZ, SEXP spX, SEXP spMaf, SEXP spParNull);

#ifdef __cplusplus
}
#endif

#endif


