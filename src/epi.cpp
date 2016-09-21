
#include <Rcpp.h>
#include "misc.h"
#include "probability.h"
#include "epi.h"

using namespace std;

//------------------------------------------------
// draw from asynchronous stochastic SIS model. Return infectives at all time points at which new infection or recovery occurs. Stop when maxIterations is reached.
// [[Rcpp::export]]
Rcpp::List SIS_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double r = Rcpp::as<double>(args["r"]);
    int I_start = Rcpp::as<int>(args["I_start"]);
    int N = Rcpp::as<int>(args["N"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int I = I_start;
    double N_inv = 1/double(N);
    vector<double> t_vec(maxIterations);
    vector<double> I_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    I_vec[0] = I_start;
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rate1 = beta*I*(1-I*N_inv);
        rate2 = r*I;
        rateTotal = rate1+rate2;
        
        // draw new time
        t += rexp1(rateTotal);
        t_vec[i] = t;
        
        // draw event
        rand1 = unif_rand();
        if (rand1<(rate1/rateTotal)) {
            I++;
        } else {
            I--;
        }
        
        // store values
        I_vec[i] = I;
        
        // abort on extinction
        if (I==0)
            break;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("t")=t_vec, Rcpp::Named("I")=I_vec);
    
}

//------------------------------------------------
// draw from synchronous stochastic SIS model. Return infectives at defined time points.
// [[Rcpp::export]]
Rcpp::List SIS_stochastic_sync_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double r = Rcpp::as<double>(args["r"]);
    int I_start = Rcpp::as<int>(args["I_start"]);
    int N = Rcpp::as<int>(args["N"]);
    vector<double> t_vec = Rcpp::as<vector<double> >(args["t_vec"]);
    
    // setup some initial parameters
    int I = I_start;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> I_vec(t_size);
    I_vec[0] = I;
    
    // carry out simulation loop
    double rate1, rate2, prob1, prob2, delta_t;
    int rand1, rand2;
    for (int i=1; i<t_size; i++) {
        
        // calculate rates of all events
        rate1 = beta*I*N_inv;
        rate2 = r;
        
        // convert to probabilities
        delta_t = t_vec[i]-t_vec[i-1];
        prob1 = 1 - exp(-rate1*delta_t);
        prob2 = 1 - exp(-rate2*delta_t);
        
        // draw events
        rand1 = rbinom1(N-I, prob1); // new infections
        rand2 = rbinom1(I, prob2); // recoveries
        I += rand1 - rand2;
        
        // store values
        I_vec[i] = I;
        
        // abort on extinction
        if (I==0)
            break;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("I")=I_vec);
    
}

//------------------------------------------------
// draw from stochastic SIS model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
// [[Rcpp::export]]
Rcpp::List SIS_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double r = Rcpp::as<double>(args["r"]);
    int I_start = Rcpp::as<int>(args["I_start"]);
    int N = Rcpp::as<int>(args["N"]);
    vector<double> t_vec = Rcpp::as<vector<double> >(args["t_vec"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int I = I_start;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> I_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    int i=0, j=0;
    while (i<maxIterations) { // while loop means i has value even after loop completes
        i++;
        
        // calculate rates of all events
        rate1 = beta*I*(1-I*N_inv);
        rate2 = r*I;
        rateTotal = rate1+rate2;
        
        // draw new time
        t += rexp1(rateTotal);
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size)
                break;
            I_vec[j] = I;
            j++;
        }
        
        // draw event
        rand1 = unif_rand();
        if (rand1<(rate1/rateTotal)) {
            I++;
        } else {
            I--;
        }
        
        // abort on end of t_vec
        if (j==t_size)
            break;
        
    } // end of simulation loop
    
    // report if maxIterations reached
    if (i==maxIterations)
        cout << "maxIterations reached\n";
    
    // return values
    return Rcpp::List::create(Rcpp::Named("I")=I_vec);
    
}

//------------------------------------------------
// draw from asynchronous stochastic SIR model. Return state of the system at all time points at which an event occurs. Stop when maxIterations is reached. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
// [[Rcpp::export]]
Rcpp::List SIR_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double r = Rcpp::as<double>(args["r"]);
    double mu = Rcpp::as<double>(args["mu"]);
    int I_init = Rcpp::as<int>(args["I_init"]);
    int R_init = Rcpp::as<int>(args["R_init"]);
    int N = Rcpp::as<int>(args["N"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    vector<double> t_vec(maxIterations);
    vector<double> S_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    vector<double> I_vec(maxIterations, -1);
    vector<double> R_vec(maxIterations, -1);
    S_vec[0] = S;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rateTotal;
    vector<double> rates(5);
    int rand1;
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rates[0] = beta*S*I*N_inv; // infection
        rates[1] = r*I; // recovery
        rates[2] = mu*S; // natural death in S
        rates[3] = mu*I; // natural death in I
        rates[4] = mu*R; // natural death in R
        rateTotal = sum(rates);
        
        // abort if rateTotal==0 (system has reached a stationary point)
        if (rateTotal<UNDERFLO)
            break;
        
        // draw new time
        t += rexp1(rateTotal);
        t_vec[i] = t;
        
        // draw event
        rand1 = sample1(rates, rateTotal);
        switch(rand1) {
            case 1: // infection
                S --;
                I ++;
                break;
            case 2: // recovery
                I --;
                R ++;
                break;
            case 3: // natural death in S
                // (death exactly matched by new birth)
                break;
            case 4: // natural death in I
                I --;
                S ++;
                break;
            case 5: // natural death in R
                R --;
                S ++;
        }
        
        // store values
        S_vec[i] = S;
        I_vec[i] = I;
        R_vec[i] = R;
        
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("t")=t_vec,
                              Rcpp::Named("S")=S_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
    
}

//------------------------------------------------
// draw from synchronous stochastic SIR model. Return state of the system at all time points at which an event occurs. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
// [[Rcpp::export]]
Rcpp::List SIR_stochastic_sync_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double r = Rcpp::as<double>(args["r"]);
    double mu = Rcpp::as<double>(args["mu"]);
    int I_init = Rcpp::as<int>(args["I_init"]);
    int R_init = Rcpp::as<int>(args["R_init"]);
    int N = Rcpp::as<int>(args["N"]);
    vector<double> t_vec = Rcpp::as<vector<double> >(args["t_vec"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> S_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> I_vec(t_size, -1);
    vector<double> R_vec(t_size, -1);
    S_vec[0] = S;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double delta_t;
    double rate_inf;
    vector<double> probs(5);
    vector<int> rands(6);
    for (int i=1; i<t_size; i++) {
        
        // calculate rates of all events
        rate_inf = beta*I*N_inv; // infection
        
        // convert to probabilities, allowing for competing hazards
        delta_t = t_vec[i]-t_vec[i-1];
        probs[0] = 1 - exp(-(rate_inf + mu)*delta_t); // infection or death in S
        probs[1] = rate_inf/(rate_inf+mu); // infection
        probs[2] = 1 - exp(-(r + mu)*delta_t); // recovery or death in I
        probs[3] = r/(r+mu); // recovery
        probs[4] = 1 - exp(-mu*delta_t); // death in R
        
        // draw events
        rands[0] = rbinom1(S, probs[0]); // infection or death in S
        rands[1] = rbinom1(rands[0], probs[1]); // infection
        rands[2] = rbinom1(I, probs[2]); // recovery or death in I
        rands[3] = rbinom1(rands[2], probs[3]); // recovery
        rands[4] = rands[2] - rands[3]; // death in I
        rands[5] = rbinom1(R, probs[4]); // death in R
        
        S += -rands[1] + rands[4] + rands[5];
        I += rands[1] - rands[3] - rands[4];
        R += rands[3] - rands[5];
        
        // store values
        S_vec[i] = S;
        I_vec[i] = I;
        R_vec[i] = R;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("S")=S_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
    
}

//------------------------------------------------
// Draw from stochastic SIR model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
// [[Rcpp::export]]
Rcpp::List SIR_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double r = Rcpp::as<double>(args["r"]);
    double mu = Rcpp::as<double>(args["mu"]);
    int I_init = Rcpp::as<int>(args["I_init"]);
    int R_init = Rcpp::as<int>(args["R_init"]);
    int N = Rcpp::as<int>(args["N"]);
    vector<double> t_vec = Rcpp::as<vector<double> >(args["t_vec"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> S_vec(t_size, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    vector<double> I_vec(t_size, -1);
    vector<double> R_vec(t_size, -1);
    S_vec[0] = S;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rateTotal;
    vector<double> rates(5);
    int rand1;
    int i=0, j=0;
    while (i<maxIterations) { // while loop means i has value even after loop completes
        i++;
        
        // calculate rates of all events
        rates[0] = beta*S*I*N_inv; // infection
        rates[1] = r*I; // recovery
        rates[2] = mu*S; // natural death in S
        rates[3] = mu*I; // natural death in I
        rates[4] = mu*R; // natural death in R
        rateTotal = sum(rates);
        
        // abort if rateTotal==0 (system has reached a stationary point)
        if (rateTotal<UNDERFLO)
            break;
        
        // draw new time
        t += rexp1(rateTotal);
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size)
                break;
            S_vec[j] = S;
            I_vec[j] = I;
            R_vec[j] = R;
            j++;
        }
        
        // draw event
        rand1 = sample1(rates, rateTotal);
        switch(rand1) {
            case 1: // infection
                S --;
                I ++;
                break;
            case 2: // recovery
                I --;
                R ++;
                break;
            case 3: // natural death in S
                // (death exactly matched by new birth)
                break;
            case 4: // natural death in I
                I --;
                S ++;
                break;
            case 5: // natural death in R
                R --;
                S ++;
        }
        
        // abort on end of t_vec
        if (j==t_size)
            break;
        
        
    } // end of simulation loop
    
    // report if maxIterations reached
    if (i==maxIterations)
        cout << "maxIterations reached\n";
    
    // return values
    return Rcpp::List::create(Rcpp::Named("S")=S_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
    
}