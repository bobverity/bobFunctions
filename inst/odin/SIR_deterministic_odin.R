
# derivatives
deriv(S) <- -beta*S*I/N + mu*(S+I+R) - mu*S
deriv(I) <- beta*S*I/N - r*I - mu*I
deriv(R) <- r*I - mu*R

# initial conditions
initial(S) <- N - I_init - R_init
initial(I) <- I_init
initial(R) <- R_init

# parameters
beta <- user()
r <- user()
mu <- user()
I_init <- user()
R_init <- user()
N <- user()