generated quantities {
  array[T] int<lower=0, upper=N> S;  // Simulated susceptible individuals at each time step
  array[T] int<lower=0, upper=N> R;  // Simulated recovered individuals at each time step
  
  // Initial values for S and R
  S[1] = N - I[1];
  R[1] = 0;
  
  for (t in 2:T) {
    // Simulate susceptible individuals
    S[t] = S[t-1] - poisson_rng(beta * S[t-1] * I[t-1] / N);
    
    // Simulate recovered individuals
    R[t] = R[t-1] + poisson_rng(gamma * I[t-1]);
  }
}
