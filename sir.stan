data {
  int<lower=0> N;         // Total population size
  int<lower=0> T;         // Number of time steps
  array[T] int<lower=0, upper=N> I;  // Infectious individuals at each time step
}

parameters {
  real<lower=0> beta;     // Transmission rate
  real<lower=0> gamma;    // Recovery rate
}

model {
  // Prior distributions for parameters
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  
  I[1] = 0;
  
  for (t in 2:T) {
    // Infectious individuals follow a Poisson distribution
    I[t] ~ poisson((beta * (N - I[t-1]) * I[t-1] / N) - gamma * I[t-1]);
  }
}