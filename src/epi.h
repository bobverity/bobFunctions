
#ifndef __bobFunctions__epi__
#define __bobFunctions__epi__

//#include <string.h>

//------------------------------------------------
// Draw from asynchronous stochastic SIS model. Return infectives at all time points at which new infection or recovery occurs. Stop when maxIterations is reached.
Rcpp::List SIS_stochastic_async_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from stochastic SIS model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
Rcpp::List SIS_stochastic_hybrid_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from synchronous stochastic SIS model. Return infectives at defined time points.
Rcpp::List SIS_stochastic_sync_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from asynchronous stochastic SIR model. Return state of the system at all time points at which an event occurs. Stop when maxIterations is reached. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
Rcpp::List SIR_stochastic_async_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from stochastic SIR model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
Rcpp::List SIR_stochastic_hybrid_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from synchronous stochastic SIR model. Return state of the system at all time points at which an event occurs. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
Rcpp::List SIR_stochastic_sync_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from asynchronous stochastic SLIR model, where L is an incubation (lag) stage of defined length. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
Rcpp::List SLIR_stochastic_async_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from stochastic SLIR model, where L is an incubation (lag) stage of defined length, using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
Rcpp::List SLIR_stochastic_hybrid_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from asynchronous stochastic Ross-Macdonald model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
Rcpp::List RM1_stochastic_async_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from stochastic Ross-Macdonald model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
Rcpp::List RM1_stochastic_hybrid_cpp(Rcpp::List args);

//------------------------------------------------
// Draw from synchronous stochastic Ross-Macdonald model. Return state of the system at known time points. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
Rcpp::List RM1_stochastic_sync_cpp(Rcpp::List args);

#endif