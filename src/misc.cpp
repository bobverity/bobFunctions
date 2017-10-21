
#include <Rcpp.h>
#include <string>
#include "misc.h"

using namespace std;

//------------------------------------------------
// global variable very large/small number for catching over/underflow problems
// OVERFLO (defined in header)
// UNDERFLO (defined in header)

//------------------------------------------------
// basic sum over elements in a vector, templated for different data types
// sum (defined in header)

//------------------------------------------------
// mean of vector, templated for different data types
// mean (defined in header)

//------------------------------------------------
// min of vector, templated for different data types
// min (defined in header)

//------------------------------------------------
// max of vector, templated for different data types
// max (defined in header)

//------------------------------------------------
// push back multiple values to vector, templated for different data types
// push_back_multiple (defined in header)

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB) {
    if (logA-logB > 100) {
        return(logA);
    } else if (logB-logA > 100) {
        return(logB);
    }
    double output = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
    return output;
}

//------------------------------------------------
// helper function for printing a single value or series of values, templated for different data types
// print (defined in header)

//------------------------------------------------
// helper function for printing contents of a vector, templated for different data types
// printVector (defined in header)

//------------------------------------------------
// helper function for printing contents of a matrix, templated for different data types
// printMatrix (defined in header)

//------------------------------------------------
// helper function for printing contents of a 3D array, templated for different data types
// printArray (defined in header)

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void printStars(int n, string title) {
    Rcpp::Rcout << title;
    for (int i=0; i<n; i++) {
        Rcpp::Rcout << "*";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// analogue of R function seq() for integers
vector<int> seq_int(int from, int to, int by) {
    int n = floor((to-from)/double(by)) + 1;
    vector<int> ret(n,from);
    for (int i=1; i<n; i++) {
        from += by;
        ret[i] = from;
    }
    return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to bool format.
int Rcpp_to_bool(SEXP x) {
    return Rcpp::as<bool>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int Rcpp_to_int(SEXP x) {
    return Rcpp::as<int>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double Rcpp_to_double(SEXP x) {
    return Rcpp::as<double>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to string format.
string Rcpp_to_string(SEXP x) {
    return Rcpp::as<string>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<bool> format.
vector<bool> Rcpp_to_vector_bool(SEXP x) {
    return Rcpp::as<vector<bool> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
vector<int> Rcpp_to_vector_int(SEXP x) {
    return Rcpp::as<vector<int> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
vector<double> Rcpp_to_vector_double(SEXP x) {
    return Rcpp::as<vector<double> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
vector<string> Rcpp_to_vector_string(SEXP x) {
    return Rcpp::as<vector<string> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<bool>> format.
vector< vector<bool> > Rcpp_to_mat_bool(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<bool> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<bool> >(x[i]);
    }
    return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
vector< vector<int> > Rcpp_to_mat_int(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<int> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<int> >(x[i]);
    }
    return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
vector< vector<double> > Rcpp_to_mat_double(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<double> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<double> >(x[i]);
    }
    return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<string>> format.
vector< vector<string> > Rcpp_to_mat_string(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<string> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<string> >(x[i]);
    }
    return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<bool>>> format.
vector< vector< vector<bool> > > Rcpp_to_array_bool(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<bool> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<bool> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<bool> >(x_i[j]);
        }
    }
    return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
vector< vector< vector<double> > > Rcpp_to_array_double(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<double> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<double> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<double> >(x_i[j]);
        }
    }
    return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
vector< vector< vector<int> > > Rcpp_to_array_int(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<int> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<int> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<int> >(x_i[j]);
        }
    }
    return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<string>>> format.
vector< vector< vector<string> > > Rcpp_to_array_string(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<string> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<string> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<string> >(x_i[j]);
        }
    }
    return ret;
}

