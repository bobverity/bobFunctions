
#ifndef __bobFunctions__misc__
#define __bobFunctions__misc__

//------------------------------------------------
// define very small number for catching underflow problems
#define UNDERFLO   1e-100

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types).
template<class TYPE>
TYPE sum(std::vector<TYPE> &x) {
    TYPE output = 0;
    for (int i=0; i<int(x.size()); i++)
        output += x[i];
    return(output);
}

//------------------------------------------------
// mean of vector (templated for different data types)
template<class TYPE>
double mean(std::vector<TYPE> &x) {
    return(sum(x)/double(x.size()));
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double logSum(double logA, double logB);

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
template<class TYPE>
void print(TYPE x) {
    std::cout << x << "\n";
}

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
template<class TYPE>
void printVector(std::vector<TYPE> &x) {
    for (int i=0; i<x.size(); i++) {
        std::cout << x[i] << " ";
    }
    std::cout << "\n";
}

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
template<class TYPE>
void printMatrix(std::vector< std::vector<TYPE> > &M) {
    for (int i=0; i<M.size(); i++) {
        for (int j=0; j<M[i].size(); j++) {
            std::cout << M[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
template<class TYPE>
void printArray(std::vector< std::vector< std::vector<TYPE> > > &x) {
    for (int i=0; i<x.size(); i++) {
        std::cout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<x[i].size(); j++) {
            for (int k=0; k<x[i][j].size(); k++) {
                std::cout << x[i][j][k] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


#endif