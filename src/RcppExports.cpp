// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dummy1_cpp
void dummy1_cpp();
RcppExport SEXP _bobFunctions_dummy1_cpp() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    dummy1_cpp();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bobFunctions_dummy1_cpp", (DL_FUNC) &_bobFunctions_dummy1_cpp, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bobFunctions(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
