#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP c_quantile2(SEXP, SEXP);
extern SEXP running_quantile_windowed(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP running_quantile_windowed_bootstrap(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"c_quantile2",                         (DL_FUNC) &c_quantile2,                         2},
    {"running_quantile_windowed",           (DL_FUNC) &running_quantile_windowed,           5},
    {"running_quantile_windowed_bootstrap", (DL_FUNC) &running_quantile_windowed_bootstrap, 5},
    {NULL, NULL, 0}
};

void R_init_climdex_pcic(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
