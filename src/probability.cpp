
#include <Rcpp.h>
#include <random>
#include "probability.h"
#include "misc.h"

using namespace std;

//-- set random seed --
random_device rd;
default_random_engine generator(rd());

uniform_real_distribution<double> uniform_0_1(0.0,1.0);

//------------------------------------------------
// draw from uniform(a,b) distribution, where a <= x < b
double runif1(double a, double b) {
    uniform_real_distribution<double> uniform_a_b(a,b);
    return(uniform_a_b(generator));
}

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(vector<double> &p, double pSum) {
    double rand = pSum*uniform_0_1(generator);
    double z = 0;
    for (int i=0; i<int(p.size()); i++) {
        z += p[i];
        if (rand<z)
            return i+1;
    }
    return(0);
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal probability
int sample2(int a, int b) {
    int z = floor(runif1(a, b+1));
    return(z);
}

//------------------------------------------------
// sample without replacement from vector x
vector<int> sample3(vector<int> x, int n) {
    vector<int> output; // faster if initialised at length n?
    int rand1;
    for (int i=0; i<n; i++) {
        rand1 = floor(runif1(0,x.size())); //replace with sample2 once written?
        output.push_back(x[rand1]);
        x.erase(x.begin()+rand1);
    }
    return(output);
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
int rbernoulli1(double p) {
    bernoulli_distribution dist_bernoulli(p);
    return(dist_bernoulli(generator));
}

//------------------------------------------------
// draw from binomial(n,p) distribution
int rbinom1(int n, double p) {
    binomial_distribution<int> dist_binomial(n,p);
    return(dist_binomial(generator));
}

//------------------------------------------------
// draw from Poisson(rate) distribution
int rpois1(double rate) {
    poisson_distribution<int> dist_poisson(rate);
    return(dist_poisson(generator));
}

//------------------------------------------------
// draw from hypergeometric distribution. Returns number of white balls given N total balls, K total white balls and n samples (without replacement). Chooses between two different sampling methods depending on parameter values.
int rhyper1(int N, int K, int n) {
    int x = (n<10*(1+log(n))) ? rhyper1_method1(N,K,n) : rhyper1_method2(N,K,n);
    return(x);
}
// method1 simply draws from the sampling process
int rhyper1_method1(int N, int K, int n) {
    int x = 0;
    for (int i=0; i<n; i++) {
        if (runif1(0,1)<(K-x)/double(N-i))
            x++;
    }
    return(x);
}
// method2 uses inverse CDF method. Starts searching at the mode and explores outwards in either direction until a hit is obtained
int rhyper1_method2(int N, int K, int n) {
    // draw uniform value
    double u = runif1(0,1);
    // start x at mode
    int x = (n+1)*(K+1)/(N+2);
    // calculate hypergeometric density at this point
    double p = dhyper1(N,K,n,x);
    // if hit then done
    if (u<p)
        return(x);
    
    // otherwise define left and right versions and store cumulative density s
    int x_left = x;
    int x_right = x;
    double p_left = p;
    double p_right = p;
    double s = p;
    bool goLeft = (x>0);
    
    // search alternating left then right. pmf in either case is calculated from previous values using a simple recursion
    for (int i=0; i<(N+1); i++) {
        if (goLeft) {
            x_left--;
            p_left /= (K-x_left)*(n-x_left)/double((x_left+1)*(N-K-n+x_left+1));
            s += p_left;
            // if hit then done
            if (u<s)
                return(x_left);
            if (x_right<N)
                goLeft = false;
        } else {
            p_right *= (K-x_right)*(n-x_right)/double((x_right+1)*(N-K-n+x_right+1));
            x_right++;
            s += p_right;
            // if hit then done
            if (u<s)
                return(x_right);
            if (x_left>0)
                goLeft = true;
        }
    }
    // return -1 if some kind of error
    return(-1);
}

//------------------------------------------------
// pmf of hypergeometric distribution
double dhyper1(int N, int K, int n, int x) {
    double z1 = lgamma(K+1)-lgamma(x+1)-lgamma(K-x+1);
    double z2 = lgamma(N-K+1)-lgamma(n-x+1)-lgamma(N-K-(n-x)+1);
    double z3 = lgamma(N+1)-lgamma(n+1)-lgamma(N-n+1);
    return(exp(z1+z2-z3));
}

//------------------------------------------------
// draw from exponential(rate) distribution (expectation = 1/rate). Found to be faster than the Rf_rexp(1/rate) method.
double rexp1(double rate) {
    exponential_distribution<double> dist_exponential(rate);
    return(dist_exponential(generator));
}

//------------------------------------------------
// draw from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate) {
    gamma_distribution<double> gamma_dist(shape,1.0/rate);
    double x = gamma_dist(generator);
    
    // check for zero or infinite values (corrects for bug in some compilers)
    while (x==0 || (1.0/x)==0)
        x = gamma_dist(generator);
    
    return(x);
}

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta) {
    double X1 = rgamma1(alpha,1.0);
    double X2 = rgamma1(beta,1.0);
    return(X1/(X1+X2));
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
    normal_distribution<double> normal_dist(mean,sd);
    return(normal_dist(generator));
}
