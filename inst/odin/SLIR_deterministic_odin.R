
# delay states
S_delay <- delay(S, dur_lag, 0)
I_delay <- delay(I, dur_lag, 0)

# derivatives
deriv(S) <- -beta*S*I/N
deriv(L) <- beta*S*I/N - beta*S_delay*I_delay/N
deriv(I) <- beta*S_delay*I_delay/N - r*I
deriv(R) <- r*I

# initial conditions
initial(S) <- N - I_init - R_init
initial(L) <- 0
initial(I) <- I_init
initial(R) <- R_init

# parameters
beta <- user()
dur_lag <- user()
r <- user()
I_init <- user()
R_init <- user()
N <- user()