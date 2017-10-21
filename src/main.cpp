
#include <Rcpp.h>
#include "misc.h"

using namespace std;

//------------------------------------------------
// dummy function
// [[Rcpp::export]]
void dummy1_cpp() {
    print("Rcpp function working");
}
