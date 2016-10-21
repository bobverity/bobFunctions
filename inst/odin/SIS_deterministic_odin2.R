
# derivatives
deriv(S) <- -beta*S*I/N + r*I
deriv(I) <- beta*S*I/N - r*I

# initial conditions
initial(S) <- N - I_init
initial(I) <- I_init

# parameters
beta <- user()
r <- user()
I_init <- user()
N <- user()
