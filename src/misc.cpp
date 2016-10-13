
#include <Rcpp.h>
#include "misc.h"
using namespace Rcpp;

//------------------------------------------------
// define very small number for catching underflow problems
// DEFINED IN HEADER

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// mean of vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB) {
    if (logA-logB > 100) {
        return(logA);
    } else if (logB-logA > 100) {
        return(logB);
    }
    double output = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
    return(output);
}

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
// DEFINED IN HEADER

//------------------------------------------------
// force exits R, which has the advantage that windows will appear automatically on re-opening
// [[Rcpp::export]]
void exit_cpp() {
    exit(0);
}