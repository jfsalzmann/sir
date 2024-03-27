functions {
  vector sir_ode(real t_unused, vector state, real beta, real gamma, real N) {
    real S = state[1];
    real I = state[2];
    real R = state[3];
    
    vector[3] d_dt;
    d_dt[1] = -beta * I * S / N;                 // dS_dt
    d_dt[2] =  beta * I * S / N - gamma * I;     // dI_dt
    d_dt[3] = -(d_dt[1] + d_dt[2]);              // dR_dt = gamma * I
    return d_dt;
  }
}
data {
  int<lower=1> N;                    // Total population size
  int<lower=1> T;                    // Number of time steps
  array[T] int<lower=0, upper=N> I;  // Infectious individuals at each time step
}
transformed data {
  int T_0 = 1;
  int T_head = 4;
  array[T-1] int T_ev = linspaced_int_array(T-1,2,T);
  array[T-1] int I_ev = I[2:T];
  array[T_head] real<lower=0> I_head = I[1:T_head];
}
parameters {
  real<lower=0.001> I_0;
  real<lower=0.001> beta;
  real<lower=0.001> gamma;
  real<lower=0.001> phi_inv;
  //real<lower=0.001> sigma;
}
transformed parameters{
  real<lower=0> S_0 = N-I_0;
  real<lower=0> R_0 = 0.001;
  vector[3] state_0 = [S_0,I_0,R_0]';
  
  array[T-1] vector[3] states_unbound = ode_rk45(sir_ode, state_0, T_0, T_ev, beta, gamma, N);
  
  array[T-1] vector<lower=0>[3] states;
  for(t in 1:T-1){
    states[t] = [fmax(states_unbound[t,1],0.001), fmax(states_unbound[t,2],0.001), fmax(states_unbound[t,3],0.001)]';
  }
}
model {
  // priors
  beta ~ gamma(2, 0.3);
  gamma ~ gamma(2, 0.1);
  I_head ~ normal(I_0,1);
  I_0 ~ normal(1,0.3);
  
  phi_inv ~ exponential(5);
  //sigma ~ lognormal(0,1);

  // likelihood
  I_ev ~ neg_binomial_2(states[:,2], 1./phi_inv+0.001);
  //I_ev ~ normal(states[:,2], sigma);
}
generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  array[T-1] real pred_cases = neg_binomial_2_rng(states[:,2], 1./phi_inv);
  //array[T-1] real pred_cases = normal_rng(states[:,2], sigma);
}
