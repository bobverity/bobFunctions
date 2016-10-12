
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
 mu : mosquito instantaneous death rate. mu = -log(p) unless specified.
 u : intrinsic incubation period. The number of days from infection to infectiousness in a human host.
 v : extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
 r : daily recovery rate.
 b : probability a human becomes infected after being bitten by an infected mosquito.
 c : probability a mosquito becomes infected after biting an infected human.
 Eh_init : initial number of infected but not infectious humans.
 Ih_init : initial number of infectious humans.
 Em_init : initial number of infected but not yet infectious mosquitoes.
 Im_init : initial number of infectious mosquitoes.
 H : human population size.
 M : mosquito population size (number of adult female mosquitoes).
 maxIterations : exit if this number of iterations is reached.
*/
// [[Rcpp::export]]
Rcpp::List RM1_stochastic_async_cpp(Rcpp::List args) {
    
    // convert input format
    double a = Rcpp::as<double>(args["a"]);
    double mu = Rcpp::as<double>(args["mu"]);
    double u = Rcpp::as<double>(args["u"]);
    double v = Rcpp::as<double>(args["v"]);
    double r = Rcpp::as<double>(args["r"]);
    double b = Rcpp::as<double>(args["b"]);
    double c = Rcpp::as<double>(args["c"]);
    int Eh_init = Rcpp::as<int>(args["Eh_init"]);
    int Ih_init = Rcpp::as<int>(args["Ih_init"]);
    int Em_init = Rcpp::as<int>(args["Em_init"]);
    int Im_init = Rcpp::as<int>(args["Im_init"]);
    int H = Rcpp::as<int>(args["H"]);
    int M = Rcpp::as<int>(args["M"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int Sh = H - Eh_init - Ih_init;
    int Eh = Eh_init;
    int Ih = Ih_init;
    int Sm = M - Em_init - Im_init;
    int Em = Em_init;
    int Im = Im_init;
    double H_inv = 1/double(H);
    vector<double> t_vec(maxIterations);
    vector<double> Sh_vec(maxIterations, -1); // -1 acts as an indicator that these values should be trimmed from the final output in the R function
    vector<double> Eh_vec(maxIterations, -1);
    vector<double> Ih_vec(maxIterations, -1);
    vector<double> Sm_vec(maxIterations, -1);
    vector<double> Em_vec(maxIterations, -1);
    vector<double> Im_vec(maxIterations, -1);
    vector<string> event_vec(maxIterations);
    Sh_vec[0] = Sh;
    Eh_vec[0] = Eh;
    Ih_vec[0] = Ih;
    Sm_vec[0] = Sm;
    Em_vec[0] = Em;
    Im_vec[0] = Im;
    event_vec[0] = "start";
    
    // if there are members of Eh and Em states then populate lists
    vector<double> Eh_list, Em_list;
    for (int i=0; i<Eh; i++) {
        Eh_list.push_back(u);
    }
    for (int i=0; i<Em; i++) {
        Em_list.push_back(v);
    }
    
    // carry out simulation loop
    double t=0, rateTotal, rand1;
    vector<double> rates(5);
    int event_type, randint1, randint2;
    bool lag1, lag2;
    string event_name = "";
    for (int i=1; i<maxIterations; i++) {
        
        // calculate rates of all events
        rates[0] = a*b*Sh*Im*H_inv; // human infection (move to Eh state)
        rates[1] = r*Ih; // human recovery
        
        rates[2] = a*c*Ih*Sm/H; // mosquito infection (move to Em state)
        rates[3] = mu*Em; // mosquito death in Em state
        rates[4] = mu*Im; // mosquito death in Im state
        rateTotal = sum(rates);
        
        // draw new time
        rand1 = rexp1(rateTotal);
        
        // work out whether new time goes past next human or mosquito incubation time
        lag1 = false;
        lag2 = false;
        if (Eh>0)
            lag1 = (t+rand1)>=Eh_list[0];
        if (Em>0)
            lag2 = (t+rand1)>=Em_list[0];
        
        // choose event type (1=new event, 2=Eh->Ih, 3=Em->Im)
        event_type = 1;
        if (lag1 && lag2) {
            if (Eh_list[0]<Em_list[0]) {
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
                case 1: // human infection (move to Eh state)
                    Sh --;
                    Eh ++;
                    Eh_list.push_back(t + u);
                    event_name = "human infection";
                    break;
                case 2: // human recovery
                    Ih --;
                    Sh ++;
                    event_name = "human recovery";
                    break;
                case 3: // mosquito infection (move to Em state)
                    Sm --;
                    Em ++;
                    Em_list.push_back(t + v);
                    event_name = "mosquito infection";
                    break;
                case 4: // mosquito death in Em state
                    randint2 = sample2(Em);
                    Em_list.erase(Em_list.begin()+randint2-1);
                    Em --;
                    Sm ++; // (new birth)
                    event_name = "mosquito death in Em";
                    break;
                case 5: // mosquito death in Im state
                    Im --;
                    Sm ++; // (new birth)
                    event_name = "mosquito death in Im";
            }
        } else if (event_type==2) {
            t = Eh_list[0];
            Eh --;
            Ih ++;
            Eh_list.erase(Eh_list.begin());
            event_name = "human progression to Ih";
        } else {
            t = Em_list[0];
            Em --;
            Im ++;
            Em_list.erase(Em_list.begin());
            event_name = "mosquito progression to Im";
        }
        
        // abort if system has reached a stationary point
        if (Eh==0 && Ih==0 && Em==0 && Im==0)
            break;
        
        // store values
        t_vec[i] = t;
        Sh_vec[i] = Sh;
        Eh_vec[i] = Eh;
        Ih_vec[i] = Ih;
        Sm_vec[i] = Sm;
        Em_vec[i] = Em;
        Im_vec[i] = Im;
        event_vec[i] = event_name;
        
    } // end of simulation loop
    
    // return values
    return Rcpp::List::create(Rcpp::Named("time")=t_vec,
                              Rcpp::Named("Sh")=Sh_vec,
                              Rcpp::Named("Eh")=Eh_vec,
                              Rcpp::Named("Ih")=Ih_vec,
                              Rcpp::Named("Sm")=Sm_vec,
                              Rcpp::Named("Em")=Em_vec,
                              Rcpp::Named("Im")=Im_vec,
                              Rcpp::Named("event")=event_vec
                              );
    
}

//------------------------------------------------
// Draw from stochastic Ross-Macdonald model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
/*
 a : human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
 p : mosquito probability of surviving one day.
 mu : mosquito instantaneous death rate. mu = -log(p) unless specified.
 u : intrinsic incubation period. The number of days from infection to infectiousness in a human host.
 v : extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
 r : daily recovery rate.
 b : probability a human becomes infected after being bitten by an infected mosquito.
 c : probability a mosquito becomes infected after biting an infected human.
 Eh_init : initial number of infected but not infectious humans.
 Ih_init : initial number of infectious humans.
 Em_init : initial number of infected but not yet infectious mosquitoes.
 Im_init : initial number of infectious mosquitoes.
 H : human population size.
 M : mosquito population size (number of adult female mosquitoes).
 times : vector of times at which output should be returned.
 maxIterations : exit if this number of iterations is reached.
 */
// [[Rcpp::export]]
Rcpp::List RM1_stochastic_hybrid_cpp(Rcpp::List args) {
    
    // convert input format
    double a = Rcpp::as<double>(args["a"]);
    double g = Rcpp::as<double>(args["mu"]);
    double u = Rcpp::as<double>(args["u"]);
    double v = Rcpp::as<double>(args["v"]);
    double r = Rcpp::as<double>(args["r"]);
    double b = Rcpp::as<double>(args["b"]);
    double c = Rcpp::as<double>(args["c"]);
    int Eh_init = Rcpp::as<int>(args["Eh_init"]);
    int Ih_init = Rcpp::as<int>(args["Ih_init"]);
    int Em_init = Rcpp::as<int>(args["Em_init"]);
    int Im_init = Rcpp::as<int>(args["Im_init"]);
    int H = Rcpp::as<int>(args["H"]);
    int M = Rcpp::as<int>(args["M"]);
    vector<double> t_vec = Rcpp::as<vector<double> >(args["t_vec"]);
    int maxIterations = Rcpp::as<int>(args["maxIterations"]);
    
    // setup some initial parameters
    int Sh = H - Eh_init - Ih_init;
    int Eh = Eh_init;
    int Ih = Ih_init;
    int Sm = M - Em_init - Im_init;
    int Em = Em_init;
    int Im = Im_init;
    double H_inv = 1/double(H);
    int t_size = int(t_vec.size());
    vector<double> Sh_vec(t_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> Eh_vec(t_size, -1);
    vector<double> Ih_vec(t_size, -1);
    vector<double> Sm_vec(t_size, -1);
    vector<double> Em_vec(t_size, -1);
    vector<double> Im_vec(t_size, -1);
    Sh_vec[0] = Sh;
    Eh_vec[0] = Eh;
    Ih_vec[0] = Ih;
    Sm_vec[0] = Sm;
    Em_vec[0] = Em;
    Im_vec[0] = Im;
    
    // if there are members of Eh and Em states then populate lists
    vector<double> Eh_list, Em_list;
    for (int i=0; i<Eh; i++) {
        Eh_list.push_back(u);
    }
    for (int i=0; i<Em; i++) {
        Em_list.push_back(v);
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
        rates[0] = a*b*Sh*Im*H_inv; // human infection (move to Eh state)
        rates[1] = r*Ih; // human recovery
        
        rates[2] = a*c*Ih*Sm/H; // mosquito infection (move to Em state)
        rates[3] = g*Em; // mosquito death in Em state
        rates[4] = g*Im; // mosquito death in Im state
        rateTotal = sum(rates);
        
        // draw new time
        rand1 = rexp1(rateTotal);
        
        // work out whether new time goes past next human or mosquito incubation time
        lag1 = false;
        lag2 = false;
        if (Eh>0)
            lag1 = (t+rand1)>=Eh_list[0];
        if (Em>0)
            lag2 = (t+rand1)>=Em_list[0];
        
        // draw new time and event type (1=new event, 2=Eh->Ih, 3=Em->Im)
        if (lag1 && lag2) {
            if (Eh_list[0]<Em_list[0]) {
                event_type = 2;
                t = Eh_list[0];
            } else {
                event_type = 3;
                t = Em_list[0];
            }
        } else if (lag1 && !lag2) {
            event_type = 2;
            t = Eh_list[0];
        } else if (!lag1 && lag2) {
            event_type = 3;
            t = Em_list[0];
        } else {
            event_type = 1;
            t += rand1;
        }
        
        // fill in up to next value of t_vec
        while (t>=t_vec[j]) {
            if (j==t_size)
                break;
            Sh_vec[j] = Sh;
            Eh_vec[j] = Eh;
            Ih_vec[j] = Ih;
            Sm_vec[j] = Sm;
            Em_vec[j] = Em;
            Im_vec[j] = Im;
            j++;
        }
        
        // draw event
        if (event_type==1) {
            randint1 = sample1(rates, rateTotal);
            switch(randint1) {
                case 1: // human infection (move to Eh state)
                    Sh --;
                    Eh ++;
                    Eh_list.push_back(t + u);
                    break;
                case 2: // human recovery
                    Ih --;
                    Sh ++;
                    break;
                case 3: // mosquito infection (move to Em state)
                    Sm --;
                    Em ++;
                    Em_list.push_back(t + v);
                    break;
                case 4: // mosquito death in Em state
                    randint2 = sample2(Em);
                    Em_list.erase(Em_list.begin()+randint2-1);
                    Em --;
                    Sm ++; // (new birth)
                    break;
                case 5: // mosquito death in Im state
                    Im --;
                    Sm ++; // (new birth)
            }
        } else if (event_type==2) {
            Eh --;
            Ih ++;
            Eh_list.erase(Eh_list.begin());
        } else {
            Em --;
            Im ++;
            Em_list.erase(Em_list.begin());
        }
        
        // abort on end of t_vec
        if (j==t_size)
            break;
        
    } // end of simulation loop
    
    // report if maxIterations reached
    if (i==maxIterations)
        cout << "maxIterations reached\n";
    
    // return values
    return Rcpp::List::create(Rcpp::Named("Sh")=Sh_vec,
                              Rcpp::Named("Eh")=Eh_vec,
                              Rcpp::Named("Ih")=Ih_vec,
                              Rcpp::Named("Sm")=Sm_vec,
                              Rcpp::Named("Em")=Em_vec,
                              Rcpp::Named("Im")=Im_vec
                              );
    
}

//------------------------------------------------
// Draw from synchronous stochastic Ross-Macdonald model. Return state of the system at known time points. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
/*
 a : human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
 p : mosquito probability of surviving one day.
 mu : mosquito instantaneous death rate. mu = -log(p) unless specified.
 u : intrinsic incubation period. The number of days from infection to infectiousness in a human host.
 v : extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
 r : daily recovery rate.
 b : probability a human becomes infected after being bitten by an infected mosquito.
 c : probability a mosquito becomes infected after biting an infected human.
 Eh_init : initial number of infected but not infectious humans.
 Ih_init : initial number of infectious humans.
 Em_init : initial number of infected but not yet infectious mosquitoes.
 Im_init : initial number of infectious mosquitoes.
 H : human population size.
 M : mosquito population size (number of adult female mosquitoes).
 times : vector of times at which output should be returned.
 */
// [[Rcpp::export]]
Rcpp::List RM1_stochastic_sync_cpp(Rcpp::List args) {
    
    // convert input format
    double a = Rcpp::as<double>(args["a"]);
    double mu = Rcpp::as<double>(args["mu"]);
    double u = Rcpp::as<double>(args["u"]);
    double v = Rcpp::as<double>(args["v"]);
    double r = Rcpp::as<double>(args["r"]);
    double b = Rcpp::as<double>(args["b"]);
    double c = Rcpp::as<double>(args["c"]);
    int Eh_init = Rcpp::as<int>(args["Eh_init"]);
    int Ih_init = Rcpp::as<int>(args["Ih_init"]);
    int Em_init = Rcpp::as<int>(args["Em_init"]);
    int Im_init = Rcpp::as<int>(args["Im_init"]);
    int H = Rcpp::as<int>(args["H"]);
    int M = Rcpp::as<int>(args["M"]);
    vector<double> times = Rcpp::as<vector<double> >(args["times"]);
    
    // setup some initial parameters
    int Sh = H-Eh_init-Ih_init;
    int Eh = Eh_init;
    int Ih = Ih_init;
    int Sm = M-Em_init-Im_init;
    int Em = Em_init;
    int Im = Im_init;
    double H_inv = 1/double(H);
    int times_size = int(times.size());
    vector<double> Sh_vec(times_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> Eh_vec(times_size, -1);
    vector<double> Ih_vec(times_size, -1);
    vector<double> Sm_vec(times_size, -1);
    vector<double> Em_vec(times_size, -1);
    vector<double> Im_vec(times_size, -1);
    Sh_vec[0] = Sh;
    Eh_vec[0] = Eh;
    Ih_vec[0] = Ih;
    Sm_vec[0] = Sm;
    Em_vec[0] = Em;
    Im_vec[0] = Im;
    
    // convert u and v lag times to integer steps
    double delta_t = times[1]-times[0];
    int u_step = round(u/delta_t);
    int v_step = round(v/delta_t);
    
    // create vectors to store lag states
    vector<double> Eh_list(u_step+1);
    vector<double> Em_list(v_step+1);
    Eh_list[0] = Eh;
    Em_list[0] = Em;
    int uBuffer_index = 0;
    int uBuffer_index_delay = 1; // note that uBuffer_index_delay is actually u_step steps BEHIND uBuffer_index, but due to the ring-buffer looping round this actually places it one step in front of uBuffer_index at all times.
    int vBuffer_index = 0;
    
    // carry out simulation loop
    double rate_h_infection, rate_m_infection;
    double prob_h_infection, prob_m_infection_or_death, prob_m_infection;
    double prob_h_recovery = 1 - exp(-r*delta_t);
    double prob_m_death = 1 - exp(-mu*delta_t);
    int h_infection, h_recovery, m_infection_or_death, m_infection, m_death_Em, m_death_Im;
    int j2;
    for (int i=1; i<times_size; i++) {
        
        // calculate rates of events
        rate_h_infection = a*b*Im*H_inv; // human infection (move to Eh state)
        rate_m_infection = a*c*Ih*H_inv; // mosquito infection (move to Em state)
        
        // convert to probabilities, allowing for competing hazards
        prob_h_infection = 1 - exp(-rate_h_infection*delta_t); // human infection (move to Eh state)
        prob_m_infection_or_death = 1 - exp(-(rate_m_infection + mu)*delta_t); // mosquito infection or death in Sm state (competing hazards)
        prob_m_infection = rate_m_infection/(rate_m_infection + mu); // mosquito infection
        
        // update ring buffer indices
        uBuffer_index = uBuffer_index==u_step ? 0 : uBuffer_index+1;
        uBuffer_index_delay = uBuffer_index_delay==u_step ? 0 : uBuffer_index_delay+1;
        vBuffer_index = vBuffer_index==v_step ? 0 : vBuffer_index+1;
        
        // human events
        h_infection = rbinom1(Sh, prob_h_infection); // human infection (move to Eh state)
        h_recovery = rbinom1(Ih, prob_h_recovery); // human recovery
        Sh += -h_infection + h_recovery; // update Sh
        
        Eh_list[uBuffer_index] = h_infection; // add infecteds to Eh list
        Eh += Eh_list[uBuffer_index] - Eh_list[uBuffer_index_delay]; // update Eh
        
        Ih += Eh_list[uBuffer_index_delay] - h_recovery; // update Ih
        
        // mosquito events
        m_infection_or_death = rbinom1(Sm, prob_m_infection_or_death); // mosquito infection or death in Sm state (competing hazards)
        m_infection = rbinom1(m_infection_or_death, prob_m_infection); // mosquito infection
        m_death_Im = rbinom1(Im, prob_m_death); // mosquito death in Im state
        Sm += -m_infection + m_death_Im; // update Sm
        
        Em_list[vBuffer_index] = m_infection; // add infecteds to Em list
        Em += m_infection; // add new infecteds to Em
        
        j2 = vBuffer_index;
        for (int j=0; j<v_step; j++) { // loop through all previous entries in Em list
            j2 = j2==0 ? v_step : j2-1;
            if (Em_list[j2]>0) {
                m_death_Em = rbinom1(Em_list[j2], prob_m_death); // mosquito death in this Em state
                Em_list[j2] -= m_death_Em; // subtract mosquito deaths from this element
                Em -= m_death_Em; // subtract mosquito deaths from Em counter
                Sm += m_death_Em; // respawn mosquitoes in Sm state
            }
        }
        Em -= Em_list[j2]; // at this stage j2 is behind vBuffer_index by v_step steps. Subtract final Em list entries from Em counter as they transition from lag state
        
        Im += Em_list[j2] - m_death_Im; // update Im
        
        // store values
        Sh_vec[i] = Sh;
        Eh_vec[i] = Eh;
        Ih_vec[i] = Ih;
        Sm_vec[i] = Sm;
        Em_vec[i] = Em;
        Im_vec[i] = Im;
        
    }
    
    // return values
    return Rcpp::List::create(Rcpp::Named("Sh")=Sh_vec,
                              Rcpp::Named("Eh")=Eh_vec,
                              Rcpp::Named("Ih")=Ih_vec,
                              Rcpp::Named("Sm")=Sm_vec,
                              Rcpp::Named("Em")=Em_vec,
                              Rcpp::Named("Im")=Im_vec
                              );
}

//------------------------------------------------
// Draw from a synchronous stochastic version of a particular Ross-Macdonald-style model (see \code{?RM2_deterministic} for details of the model). Return state of the system at known time points. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
/*
 a : human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
 mu : mosquito instantaneous  death rate.
 lambda : mosquito instantaneous birth rate.
 u : intrinsic incubation period. The number of days from infection to infectiousness in a human host.
 v : extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
 r : daily recovery rate.
 b : probability a human becomes infected after being bitten by an infected mosquito.
 c : probability a mosquito becomes infected after biting an infected human.
 Eh_init : initial number of infected but not infectious humans.
 Ih_init : initial number of infectious humans.
 Em_init : initial number of infected but not yet infectious mosquitoes.
 Im_init : initial number of infectious mosquitoes.
 H : human population size.
 M_init : initial mosquito population size.
 times : vector of times at which output should be returned.
 Ktimes : vector of times at which carrying capacity is defined.
 Kvalues : vector of carrying capacities that come into action at \code{Ktimes}.
 */
// [[Rcpp::export]]
Rcpp::List RM2_stochastic_sync_cpp(Rcpp::List args) {
    
    // convert input format
    double a = Rcpp::as<double>(args["a"]);
    double mu = Rcpp::as<double>(args["mu"]);
    double lambda = Rcpp::as<double>(args["lambda"]);
    double u = Rcpp::as<double>(args["u"]);
    double v = Rcpp::as<double>(args["v"]);
    double r = Rcpp::as<double>(args["r"]);
    double b = Rcpp::as<double>(args["b"]);
    double c = Rcpp::as<double>(args["c"]);
    int Eh_init = Rcpp::as<int>(args["Eh_init"]);
    int Ih_init = Rcpp::as<int>(args["Ih_init"]);
    int Em_init = Rcpp::as<int>(args["Em_init"]);
    int Im_init = Rcpp::as<int>(args["Im_init"]);
    int H = Rcpp::as<int>(args["H"]);
    int M_init = Rcpp::as<int>(args["M_init"]);
    vector<double> times = Rcpp::as<vector<double> >(args["times"]);
    vector<double> Ktimes = Rcpp::as<vector<double> >(args["Ktimes"]);
    vector<int> Kvalues = Rcpp::as<vector<int> >(args["Kvalues"]);
    
    // setup some initial parameters
    int Sh = H-Eh_init-Ih_init;
    int Eh = Eh_init;
    int Ih = Ih_init;
    int K = Kvalues[0];
    int M = M_init;
    int Sm = M-Em_init-Im_init;
    int Em = Em_init;
    int Im = Im_init;
    double H_inv = 1/double(H);
    int times_size = int(times.size());
    vector<double> Sh_vec(times_size, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
    vector<double> Eh_vec(times_size, -1);
    vector<double> Ih_vec(times_size, -1);
    vector<double> Sm_vec(times_size, -1);
    vector<double> Em_vec(times_size, -1);
    vector<double> Im_vec(times_size, -1);
    vector<double> M_vec(times_size, -1);
    vector<double> K_vec(times_size, -1);
    Sh_vec[0] = Sh;
    Eh_vec[0] = Eh;
    Ih_vec[0] = Ih;
    Sm_vec[0] = Sm;
    Em_vec[0] = Em;
    Im_vec[0] = Im;
    M_vec[0] = M;
    K_vec[0] = K;
    
    // convert u and v lag times to integer steps
    double delta_t = times[1]-times[0];
    int u_step = round(u/delta_t);
    int v_step = round(v/delta_t);
    
    // create vectors to store lag states
    vector<double> Eh_list(u_step+1);
    vector<double> Em_list(v_step+1);
    Eh_list[0] = Eh;
    Em_list[0] = Em;
    int uBuffer_index = 0;
    int uBuffer_index_delay = 1; // note that uBuffer_index_delay is actually u_step steps BEHIND uBuffer_index, but due to the ring-buffer looping round this actually places it one step in front of uBuffer_index at all times.
    int vBuffer_index = 0;
    
    // carry out simulation loop
    double rate_h_infection, rate_m_infection, rate_m_birth;
    double prob_h_infection, prob_Sm_infection_death, prob_Sm_infection, prob_Sm_birth, prob_EmIm_birth_death, prob_EmIm_birth;
    double prob_h_recovery = 1 - exp(-r*delta_t);
    int h_infection, h_recovery, Sm_infection_death, Sm_infection, Sm_birth, Sm_death, Em_birth_death, Em_birth, Em_death, Im_birth_death, Im_birth, Im_death;
    int j2, Kindex=0;
    int Ktimes_size = int(Ktimes.size());
    for (int i=1; i<times_size; i++) {
        
        // change carrying capacity if needed
        if (Kindex<Ktimes_size) {
            if (times[i]>Ktimes[Kindex]) {
                K = Kvalues[Kindex];
                Kindex++;
            }
        }
        
        // calculate rates of events
        rate_h_infection = a*b*Im*H_inv; // human infection (move to Eh state)
        rate_m_infection = a*c*Ih*H_inv; // mosquito infection (move to Em state)
        rate_m_birth = lambda*(1-M/double(K)); // mosquito birth in Sm state
        rate_m_birth = rate_m_birth<0 ? 0 : rate_m_birth; // don't let birth rate go below zero
        
        // convert rates to probabilities, allowing for competing hazards
        prob_h_infection = 1 - exp(-rate_h_infection*delta_t); // human infection (move to Eh state)
        prob_Sm_infection_death = 1 - exp(-(rate_m_infection + mu)*delta_t); // mosquito infection or death in Sm state (competing hazards)
        prob_Sm_infection = rate_m_infection/(rate_m_infection + mu); // mosquito infection
        prob_Sm_birth = 1 - exp(-rate_m_birth*delta_t); // mosquito birth in Sm state
        prob_EmIm_birth_death = 1 - exp(-(rate_m_birth + mu)*delta_t); // mosquito birth or death in Em or Im states (competing hazards)
        prob_EmIm_birth = rate_m_birth/(rate_m_birth + mu); // mosquito birth in Em or Im states
        
        // update ring buffer indices
        uBuffer_index = uBuffer_index==u_step ? 0 : uBuffer_index+1;
        uBuffer_index_delay = uBuffer_index_delay==u_step ? 0 : uBuffer_index_delay+1;
        vBuffer_index = vBuffer_index==v_step ? 0 : vBuffer_index+1;
        
        // human events
        h_infection = rbinom1(Sh, prob_h_infection); // human infection (move to Eh state)
        h_recovery = rbinom1(Ih, prob_h_recovery); // human recovery
        Sh += -h_infection + h_recovery; // update Sh
        
        Eh_list[uBuffer_index] = h_infection; // add infecteds to Eh list
        Eh += Eh_list[uBuffer_index] - Eh_list[uBuffer_index_delay]; // update Eh
        
        Ih += Eh_list[uBuffer_index_delay] - h_recovery; // update Ih
        
        // mosquito events
        Sm_infection_death = rbinom1(Sm, prob_Sm_infection_death); // mosquito infection or death in Sm state (competing hazards)
        Sm_infection = rbinom1(Sm_infection_death, prob_Sm_infection); // mosquito infection
        Sm_death = Sm_infection_death - Sm_infection; // mosquito death in Sm state
        Sm_birth = rbinom1(Sm-Sm_death, prob_Sm_birth); // mosquito birth in Sm state (applies to all remaining mosquitoes, inclusing those that are infected)
        
        Sm += -Sm_infection + Sm_birth - Sm_death; // update Sm
        M += Sm_birth - Sm_death; // update M with births and deaths
        
        Em_list[vBuffer_index] = Sm_infection; // add infecteds to Em list
        Em += Sm_infection; // add new infecteds to Em
        
        j2 = vBuffer_index;
        for (int j=0; j<v_step; j++) { // loop through all previous entries in Em list
            j2 = j2==0 ? v_step : j2-1;
            if (Em_list[j2]>0) {
                Em_birth_death = rbinom1(Em_list[j2], prob_EmIm_birth_death); // mosquito birth or death in this Em state
                Em_birth = rbinom1(Em_birth_death, prob_EmIm_birth); // mosquito birth by members of this Em state
                Em_death = Em_birth_death - Em_birth; // mosquito death in this Em state
                
                Em_list[j2] -= Em_death; // subtract mosquito deaths from this element
                Em -= Em_death; // subtract mosquito deaths from Em counter
                Sm += Em_birth; // add new births to Sm
                M += Em_birth - Em_death; // update M with births and deaths
            }
        }
        Em -= Em_list[j2]; // at this stage j2 is behind vBuffer_index by v_step steps. Subtract final Em list entries from Em counter as they transition from lag state
        
        Im_birth_death = rbinom1(Im, prob_EmIm_birth_death); // mosquito birth or death in Im state
        Im_birth = rbinom1(Im_birth_death, prob_EmIm_birth); // mosquito birth by members of Im state
        Im_death = Im_birth_death - Im_birth; // mosquito death in Im state
        
        Im += Em_list[j2] - Im_death; // update Im
        Sm += Im_birth; // add new births to Sm
        M += Im_birth - Im_death; // update M with births and deaths
        
        // store values
        Sh_vec[i] = Sh;
        Eh_vec[i] = Eh;
        Ih_vec[i] = Ih;
        Sm_vec[i] = Sm;
        Em_vec[i] = Em;
        Im_vec[i] = Im;
        M_vec[i] = M;
        K_vec[i] = K;
        
    }
    
    // return values
    return Rcpp::List::create(Rcpp::Named("Sh")=Sh_vec,
                              Rcpp::Named("Eh")=Eh_vec,
                              Rcpp::Named("Ih")=Ih_vec,
                              Rcpp::Named("Sm")=Sm_vec,
                              Rcpp::Named("Em")=Em_vec,
                              Rcpp::Named("Im")=Im_vec,
                              Rcpp::Named("M")=M_vec,
                              Rcpp::Named("K")=K_vec
                              );
    
}
