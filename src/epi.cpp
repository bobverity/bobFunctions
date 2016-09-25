
#include <Rcpp.h>
#include "misc.h"
#include "probability.h"
#include "epi.h"

using namespace std;

//------------------------------------------------
// Draw from asynchronous stochastic SIS model. Return infectives at all time points at which new infection or recovery occurs. Stop when maxIterations is reached.
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
// Draw from stochastic SIS model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
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
// Draw from synchronous stochastic SIS model. Return infectives at defined time points.
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
// Draw from asynchronous stochastic SIR model. Return state of the system at all time points at which an event occurs. Stop when maxIterations is reached. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
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
    vector<double> S_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
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

//------------------------------------------------
// Draw from synchronous stochastic SIR model. Return state of the system at all time points at which an event occurs. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
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
// Draw from asynchronous stochastic SLIR model, where L is an incubation (lag) stage of defined length. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
// [[Rcpp::export]]
Rcpp::List SLIR_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double dur_lag = Rcpp::as<double>(args["dur_lag"]);
    double r = Rcpp::as<double>(args["r"]);
    int I_init = Rcpp::as<int>(args["I_init"]);
    int R_init = Rcpp::as<int>(args["R_init"]);
    int N = Rcpp::as<int>(args["N"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int L = 0;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    vector<double> t_vec(maxIterations);
    vector<double> S_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    vector<double> L_vec(maxIterations, -1);
    vector<double> I_vec(maxIterations, -1);
    vector<double> R_vec(maxIterations, -1);
    S_vec[0] = S;
    L_vec[0] = 0;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    vector<double> lagList;
    int lagList_size = 0;
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rate1 = beta*S*I*N_inv; // infection (move to lag state)
        rate2 = r*I; // recovery
        rateTotal = rate1 + rate2;
        
        // draw new time
        rand1 = rexp1(rateTotal);
        if (lagList_size>0) {
            if ((t+rand1)<lagList[0]) {
                t += rand1;
            } else {
                t = lagList[0];
                lagList.erase (lagList.begin());
                lagList_size --;
                L--;
                I++;
                goto store;
            }
        } else {
            t += rand1;
        }
        
        // draw event
        rand1 = runif1();
        if (rand1<(rate1/rateTotal)) {
            S--;
            L++;
            lagList.push_back(t + dur_lag);
            lagList_size ++;
        } else {
            I--;
            R++;
        }
        
        // store values
        store:
        t_vec[i] = t;
        S_vec[i] = S;
        L_vec[i] = L;
        I_vec[i] = I;
        R_vec[i] = R;
        
        // abort if system has reached a stationary point
        if (L==0 && I==0)
            break;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("t")=t_vec,
                              Rcpp::Named("S")=S_vec,
                              Rcpp::Named("L")=L_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
    
}

//------------------------------------------------
// Draw from stochastic SLIR model, where L is an incubation (lag) stage of defined length, using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
// [[Rcpp::export]]
Rcpp::List SLIR_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double beta = Rcpp::as<double>(args["beta"]);
    double dur_lag = Rcpp::as<double>(args["dur_lag"]);
    double r = Rcpp::as<double>(args["r"]);
    int I_init = Rcpp::as<int>(args["I_init"]);
    int R_init = Rcpp::as<int>(args["R_init"]);
    int N = Rcpp::as<int>(args["N"]);
    vector<double> t_vec = Rcpp::as<vector<double> >(args["t_vec"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int S = N-I_init-R_init;
    int L = 0;
    int I = I_init;
    int R = R_init;
    double N_inv = 1/double(N);
    int t_size = int(t_vec.size());
    vector<double> S_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> L_vec(t_size, -1);
    vector<double> I_vec(t_size, -1);
    vector<double> R_vec(t_size, -1);
    S_vec[0] = S;
    L_vec[0] = L;
    I_vec[0] = I;
    R_vec[0] = R;
    
    // carry out simulation loop
    double t=0, rate1, rate2, rateTotal, rand1;
    vector<double> lagList;
    int lagList_size = 0;
    bool new_event;
    int i=0, j=0;
    while (i<maxIterations) { // while loop means i has value even after loop completes
        i++;
        
        // calculate rates of all events
        rate1 = beta*S*I*N_inv; // infection (move to lag state)
        rate2 = r*I; // recovery
        rateTotal = rate1 + rate2;
        
        // draw new time
        rand1 = rexp1(rateTotal);
        new_event = true;
        if (lagList_size>0) {
            if ((t+rand1)<lagList[0]) {
                t += rand1;
            } else {
                t = lagList[0];
                lagList.erase (lagList.begin());
                lagList_size --;
                new_event = false;
            }
        } else {
            t += rand1;
        }
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size)
                break;
            S_vec[j] = S;
            L_vec[j] = L;
            I_vec[j] = I;
            R_vec[j] = R;
            j++;
        }
        
        // draw event
        if (new_event) {
            rand1 = runif1();
            if (rand1<(rate1/rateTotal)) {
                S--;
                L++;
                lagList.push_back(t + dur_lag);
                lagList_size ++;
            } else {
                I--;
                R++;
            }
        } else {
            L--;
            I++;
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
                              Rcpp::Named("L")=L_vec,
                              Rcpp::Named("I")=I_vec,
                              Rcpp::Named("R")=R_vec);
    
}

//------------------------------------------------
// Draw from asynchronous stochastic Ross-Macdonald model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
/*
 a : human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
 p : mosquito probability of surviving one day.
 g : mosquito instantaneous death rate. g = -log(p) unless specified.
 u : intrinsic incubation period. The number of days from infection to infectiousness in a human host.
 v : extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
 r : daily recovery rate.
 b : probability a human becomes infected after being bitten by an infected mosquito.
 c : probability a mosquito becomes infected after biting an infected human.
 E_h_init : initial number of infected but not infectious humans.
 I_h_init : initial number of infectious humans.
 H : human population size.
 E_m_init : initial number of infected but not yet infectious mosquitoes.
 I_m_init : initial number of infectious mosquitoes.
 m : ratio of adult female mosquitoes to humans. Population density of adult female mosquitoes is equal to M = m*H.
 maxIterations : exit if this number of iterations is reached.
*/
// [[Rcpp::export]]
Rcpp::List RM1_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double a = Rcpp::as<double>(args["a"]);
    double g = Rcpp::as<double>(args["g"]);
    double u = Rcpp::as<double>(args["u"]);
    double v = Rcpp::as<double>(args["v"]);
    double r = Rcpp::as<double>(args["r"]);
    double b = Rcpp::as<double>(args["b"]);
    double c = Rcpp::as<double>(args["c"]);
    int E_h_init = Rcpp::as<int>(args["E_h_init"]);
    int I_h_init = Rcpp::as<int>(args["I_h_init"]);
    int H = Rcpp::as<int>(args["H"]);
    int E_m_init = Rcpp::as<int>(args["E_m_init"]);
    int I_m_init = Rcpp::as<int>(args["I_m_init"]);
    double m = Rcpp::as<double>(args["m"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int S_h = H-E_h_init-I_h_init;
    int E_h = E_h_init;
    int I_h = I_h_init;
    int M = round(m*H);
    int S_m = M-E_m_init-I_m_init;
    int E_m = E_m_init;
    int I_m = I_m_init;
    double H_inv = 1/double(H);
    vector<double> t_vec(maxIterations);
    vector<double> S_h_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    vector<double> E_h_vec(maxIterations, -1);
    vector<double> I_h_vec(maxIterations, -1);
    vector<double> S_m_vec(maxIterations, -1);
    vector<double> E_m_vec(maxIterations, -1);
    vector<double> I_m_vec(maxIterations, -1);
    vector<string> event_vec(maxIterations);
    S_h_vec[0] = S_h;
    E_h_vec[0] = E_h;
    I_h_vec[0] = I_h;
    S_m_vec[0] = S_m;
    E_m_vec[0] = E_m;
    I_m_vec[0] = I_m;
    event_vec[0] = "start";
    
    // if there are members of E_h and E_m states then populate lists
    vector<double> E_h_list, E_m_list;
    for (int i=0; i<E_h; i++) {
        E_h_list.push_back(u);
    }
    for (int i=0; i<E_m; i++) {
        E_m_list.push_back(v);
    }
    
    // carry out simulation loop
    double t=0, rateTotal, rand1;
    vector<double> rates(5);
    int event_type, randint1, randint2;
    bool lag1, lag2;
    string event_name = "";
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rates[0] = a*b*S_h*I_m*H_inv; // human infection (move to E_h state)
        rates[1] = r*I_h; // human recovery
        
        rates[2] = a*c*I_h*S_m/H; // mosquito infection (move to E_m state)
        rates[3] = g*E_m; // mosquito death in E_m state
        rates[4] = g*I_m; // mosquito death in I_m state
        rateTotal = sum(rates);
        
        // draw new time
        rand1 = rexp1(rateTotal);
        
        // work out whether new time goes past next human or mosquito incubation time
        lag1 = false;
        lag2 = false;
        if (E_h>0)
            lag1 = (t+rand1)>=E_h_list[0];
        if (E_m>0)
            lag2 = (t+rand1)>=E_m_list[0];
        
        // choose event type (1=new event, 2=E_h->I_h, 3=E_m->I_m)
        event_type = 1;
        if (lag1 && lag2) {
            if (E_h_list[0]<E_m_list[0]) {
                event_type = 2;
            } else {
                event_type = 3;
            }
        } else if (lag1 && !lag2) {
            event_type = 2;
        } else if (!lag1 && lag2) {
            event_type = 3;
        }
        
        // draw event
        if (event_type==1) {
            t += rand1;
            randint1 = sample1(rates, rateTotal);
            switch(randint1) {
                case 1: // human infection (move to E_h state)
                    S_h --;
                    E_h ++;
                    E_h_list.push_back(t + u);
                    event_name = "human infection";
                    break;
                case 2: // human recovery
                    I_h --;
                    S_h ++;
                    event_name = "human recovery";
                    break;
                case 3: // mosquito infection (move to E_m state)
                    S_m --;
                    E_m ++;
                    E_m_list.push_back(t + v);
                    event_name = "mosquito infection";
                    break;
                case 4: // mosquito death in E_m state
                    randint2 = sample2(E_m);
                    E_m_list.erase (E_m_list.begin()+randint2-1);
                    E_m --;
                    S_m ++; // (new birth)
                    event_name = "mosquito death in E_m";
                    break;
                case 5: // mosquito death in I_m state
                    I_m --;
                    S_m ++; // (new birth)
                    event_name = "mosquito death in I_m";
            }
        } else if (event_type==2) {
            t = E_h_list[0];
            E_h --;
            I_h ++;
            E_h_list.erase (E_h_list.begin());
            event_name = "human progression to I_h";
        } else {
            t = E_m_list[0];
            E_m --;
            I_m ++;
            E_m_list.erase (E_m_list.begin());
            event_name = "mosquito progression to I_m";
        }
        
        // abort if system has reached a stationary point
        if (E_h==0 && I_h==0 && E_m==0 && I_m==0)
            break;
        
        // store values
        t_vec[i] = t;
        S_h_vec[i] = S_h;
        E_h_vec[i] = E_h;
        I_h_vec[i] = I_h;
        S_m_vec[i] = S_m;
        E_m_vec[i] = E_m;
        I_m_vec[i] = I_m;
        event_vec[i] = event_name;
        
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("t")=t_vec,
                              Rcpp::Named("S_h")=S_h_vec,
                              Rcpp::Named("E_h")=E_h_vec,
                              Rcpp::Named("I_h")=I_h_vec,
                              Rcpp::Named("S_m")=S_m_vec,
                              Rcpp::Named("E_m")=E_m_vec,
                              Rcpp::Named("I_m")=I_m_vec,
                              Rcpp::Named("event")=event_vec
                              );
    
}

//------------------------------------------------
// Draw from stochastic Ross-Macdonald model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
/*
 a : human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
 p : mosquito probability of surviving one day.
 g : mosquito instantaneous death rate. g = -log(p) unless specified.
 u : intrinsic incubation period. The number of days from infection to infectiousness in a human host.
 v : extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
 r : daily recovery rate.
 b : probability a human becomes infected after being bitten by an infected mosquito.
 c : probability a mosquito becomes infected after biting an infected human.
 E_h_init : initial number of infected but not infectious humans.
 I_h_init : initial number of infectious humans.
 H : human population size.
 E_m_init : initial number of infected but not yet infectious mosquitoes.
 I_m_init : initial number of infectious mosquitoes.
 m : ratio of adult female mosquitoes to humans. Population density of adult female mosquitoes is equal to M = m*H.
 times : vector of times at which output should be returned.
 maxIterations : exit if this number of iterations is reached.
 */
// [[Rcpp::export]]
Rcpp::List RM1_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double a = Rcpp::as<double>(args["a"]);
    double g = Rcpp::as<double>(args["g"]);
    double u = Rcpp::as<double>(args["u"]);
    double v = Rcpp::as<double>(args["v"]);
    double r = Rcpp::as<double>(args["r"]);
    double b = Rcpp::as<double>(args["b"]);
    double c = Rcpp::as<double>(args["c"]);
    int E_h_init = Rcpp::as<int>(args["E_h_init"]);
    int I_h_init = Rcpp::as<int>(args["I_h_init"]);
    int H = Rcpp::as<int>(args["H"]);
    int E_m_init = Rcpp::as<int>(args["E_m_init"]);
    int I_m_init = Rcpp::as<int>(args["I_m_init"]);
    double m = Rcpp::as<double>(args["m"]);
    vector<double> t_vec = Rcpp::as<vector<double> >(args["t_vec"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int S_h = H-E_h_init-I_h_init;
    int E_h = E_h_init;
    int I_h = I_h_init;
    int M = round(m*H);
    int S_m = M-E_m_init-I_m_init;
    int E_m = E_m_init;
    int I_m = I_m_init;
    double H_inv = 1/double(H);
    int t_size = int(t_vec.size());
    vector<double> S_h_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> E_h_vec(t_size, -1);
    vector<double> I_h_vec(t_size, -1);
    vector<double> S_m_vec(t_size, -1);
    vector<double> E_m_vec(t_size, -1);
    vector<double> I_m_vec(t_size, -1);
    S_h_vec[0] = S_h;
    E_h_vec[0] = E_h;
    I_h_vec[0] = I_h;
    S_m_vec[0] = S_m;
    E_m_vec[0] = E_m;
    I_m_vec[0] = I_m;
    
    // if there are members of E_h and E_m states then populate lists
    vector<double> E_h_list, E_m_list;
    for (int i=0; i<E_h; i++) {
        E_h_list.push_back(u);
    }
    for (int i=0; i<E_m; i++) {
        E_m_list.push_back(v);
    }
    
    // carry out simulation loop
    double t=0, rateTotal, rand1;
    vector<double> rates(5);
    int event_type, randint1, randint2;
    bool lag1, lag2;
    int i=0, j=0;
    while (i<maxIterations) { // while loop means i has value even after loop completes
        i++;
        
        // calculate rates of all events
        rates[0] = a*b*S_h*I_m*H_inv; // human infection (move to E_h state)
        rates[1] = r*I_h; // human recovery
        
        rates[2] = a*c*I_h*S_m/H; // mosquito infection (move to E_m state)
        rates[3] = g*E_m; // mosquito death in E_m state
        rates[4] = g*I_m; // mosquito death in I_m state
        rateTotal = sum(rates);
        
        // draw new time
        rand1 = rexp1(rateTotal);
        
        // work out whether new time goes past next human or mosquito incubation time
        lag1 = false;
        lag2 = false;
        if (E_h>0)
            lag1 = (t+rand1)>=E_h_list[0];
        if (E_m>0)
            lag2 = (t+rand1)>=E_m_list[0];
        
        // draw new time and event type (1=new event, 2=E_h->I_h, 3=E_m->I_m)
        if (lag1 && lag2) {
            if (E_h_list[0]<E_m_list[0]) {
                event_type = 2;
                t = E_h_list[0];
            } else {
                event_type = 3;
                t = E_m_list[0];
            }
        } else if (lag1 && !lag2) {
            event_type = 2;
            t = E_h_list[0];
        } else if (!lag1 && lag2) {
            event_type = 3;
            t = E_m_list[0];
        } else {
            event_type = 1;
            t += rand1;
        }
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size)
                break;
            S_h_vec[j] = S_h;
            E_h_vec[j] = E_h;
            I_h_vec[j] = I_h;
            S_m_vec[j] = S_m;
            E_m_vec[j] = E_m;
            I_m_vec[j] = I_m;
            j++;
        }
        
        // draw event
        if (event_type==1) {
            randint1 = sample1(rates, rateTotal);
            switch(randint1) {
                case 1: // human infection (move to E_h state)
                    S_h --;
                    E_h ++;
                    E_h_list.push_back(t + u);
                    break;
                case 2: // human recovery
                    I_h --;
                    S_h ++;
                    break;
                case 3: // mosquito infection (move to E_m state)
                    S_m --;
                    E_m ++;
                    E_m_list.push_back(t + v);
                    break;
                case 4: // mosquito death in E_m state
                    randint2 = sample2(E_m);
                    E_m_list.erase (E_m_list.begin()+randint2-1);
                    E_m --;
                    S_m ++; // (new birth)
                    break;
                case 5: // mosquito death in I_m state
                    I_m --;
                    S_m ++; // (new birth)
            }
        } else if (event_type==2) {
            E_h --;
            I_h ++;
            E_h_list.erase (E_h_list.begin());
        } else {
            E_m --;
            I_m ++;
            E_m_list.erase (E_m_list.begin());
        }
        
        // abort on end of t_vec
        if (j==t_size)
            break;
        
    } // end of simulation loop
    
    // report if maxIterations reached
    if (i==maxIterations)
        cout << "maxIterations reached\n";
    
    // return values
    return Rcpp::List::create(Rcpp::Named("S_h")=S_h_vec,
                              Rcpp::Named("E_h")=E_h_vec,
                              Rcpp::Named("I_h")=I_h_vec,
                              Rcpp::Named("S_m")=S_m_vec,
                              Rcpp::Named("E_m")=E_m_vec,
                              Rcpp::Named("I_m")=I_m_vec
                              );
    
}