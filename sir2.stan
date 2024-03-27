functions {
  array[] real sir_ode(real t_unused, array[] real state, array[] real theta, 
             array[] real sir_reals, array[] int sir_ints) {
    
    real S = state[1];
    real I = state[2];
    real R = state[3];
    
    real N = sir_ints[1];
    
    real beta = theta[1];
    real gamma = theta[2];
    
    real dS_dt = -beta * I * S ;
    real dI_dt =  beta * I * S  - gamma * I;
    real dR_dt =  gamma * I;
    
    return {dS_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> N;         // Total population size
  int<lower=1> T;         // Number of time steps
  array[T] int<lower=0, upper=N> I;  // Infectious individuals at each time step
}
transformed data {
  int T_0 = 1; // T[1] 
  array[T-1] int T_ev = linspaced_int_array(T-1,2,T);
  array[0] real sir_reals;
  array[1] int sir_ints = { N };
  array[T-1] int I_ev = I[2:T];
}
parameters {
  real<lower=0.001> gamma;
  real<lower=0.001> beta;
  real<lower=0.001> phi_inv;
}
transformed parameters{
  array[T-1, 3] real states;
  
  real<lower=0> S_0 = N-I[1];
  real<lower=0> I_0 = I[1];
  real<lower=0> R_0 = 0;
  
  real phi = 1. / phi_inv;
  {
    array[2] real theta;
    theta[1] = beta;
    theta[2] = gamma;
    
    array[3] real state_0 = {S_0,I_0,R_0};
    
    states = integrate_ode_rk45(sir_ode, state_0, T_0, T_ev, theta, sir_reals, sir_ints);
  }
}
model {
  //priors
  beta ~ normal(0.2, 1);
  gamma ~ normal(0.05, 0.5);
  phi_inv ~ exponential(5);
  S_0 ~ normal(N-1,0.5);
  I_0 ~ normal(1,0.5);
  R_0 ~ normal(0,0.1);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  I_ev ~ neg_binomial_2(col(to_matrix(states), 2), phi);
}
generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  array[T-1] real pred_cases;
  pred_cases = neg_binomial_2_rng(col(to_matrix(states), 2), phi);
}
