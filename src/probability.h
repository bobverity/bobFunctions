
#ifndef __bobFunctions__probability__
#define __bobFunctions__probability__

//------------------------------------------------
// draw from uniform(a,b) distribution, where a <= x < b
double runif1(double a=0, double b=1.0);

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum=1);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability
int sample2(int b, int a=1);

//------------------------------------------------
// sample without replacement from vector x
std::vector<int> sample3(std::vector<int> x, int n);

//------------------------------------------------
// draw from Bernoulli(p) distribution
int rbernoulli1(double p);

//------------------------------------------------
// draw from binomial(n,p) distribution
int rbinom1(int n, double p);

//------------------------------------------------
// draw from Poisson(rate) distribution
int rpois1(double rate);

//------------------------------------------------
// draw from hypergeometric distribution. Returns number of white balls given N total balls, K total white balls and n samples (without replacement). Chooses between two different sampling methods depending on parameter values.
int rhyper1(int N, int K, int n);
int rhyper1_method1(int N, int K, int n);
int rhyper1_method2(int N, int K, int n);

//------------------------------------------------
// pmf of hypergeometric distribution
double dhyper1(int N, int K, int n, int x);

//------------------------------------------------
// draw from exponential(rate) distribution (expectation = 1/rate). Found to be faster than the Rf_rexp(1/rate) method.
double rexp1(double rate);

//------------------------------------------------
// sample from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd);

#endif