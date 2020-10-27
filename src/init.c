#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _OMR_EM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _OMR_NR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _OMR_rcpparma_bothproducts(SEXP);
extern SEXP _OMR_rcpparma_hello_world();
extern SEXP _OMR_rcpparma_innerproduct(SEXP);
extern SEXP _OMR_rcpparma_outerproduct(SEXP);

static const R_CallMethodDef CallEntries[] = {
	{ "_OMR_EM",                    (DL_FUNC)&_OMR_EM,                    10 },
{ "_OMR_NR",                    (DL_FUNC)&_OMR_NR,                    10 },
{ "_OMR_rcpparma_bothproducts", (DL_FUNC)&_OMR_rcpparma_bothproducts,  1 },
{ "_OMR_rcpparma_hello_world",  (DL_FUNC)&_OMR_rcpparma_hello_world,   0 },
{ "_OMR_rcpparma_innerproduct", (DL_FUNC)&_OMR_rcpparma_innerproduct,  1 },
{ "_OMR_rcpparma_outerproduct", (DL_FUNC)&_OMR_rcpparma_outerproduct,  1 },
{ NULL, NULL, 0 }
};

void R_init_OMR(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
