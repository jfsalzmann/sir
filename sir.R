################################################################################
#####
##### Bayesian Hierarchical & Latent Variable Modeling
##### SIR Model
##### Johann-Friedrich Salzmann
#####
################################################################################

# Initialisation

rm(list = ls(envir = globalenv()),envir = globalenv())

# Adjust relative path so that required R & Stan files will be found
path = "sir/"

path = ifelse(endsWith(path,"/"),path,paste0(path,"/"))
source(paste0(path,"init.R"))

################################################################################

# Step 1: Stan model

basename = "sir3"

if(file.exists(path %.% basename %.% ((Sys.info()["sysname"]=="Windows") %?% ".exe" %:% ""))) {
  model = cmdstan_model(exe_file = path %.% basename %.% ((Sys.info()["sysname"]=="Windows") %?% ".exe" %:% ""))
} else {
  model = cmdstan_model(stan_file = path %.% basename %.% ".stan")
}

# Step 2: Simulation data

set.seed(123)

# Assumptions - stable: 0.2/0.05
true_beta = 0.2
true_gamma = 0.1
t_sim = 150
n_sim = 500

# Variables
S = numeric(t_sim)
I = numeric(t_sim)
R = numeric(t_sim)

# Initial conditions
S[1] = n_sim - 1
I[1] = 1
R[1] = 0

'for (t in 2:t_sim) {
  # Compute new values of S, I, R
  dS = -true_beta * S[t-1] * I[t-1] / n_sim
  dI = true_beta * S[t-1] * I[t-1] / n_sim - true_gamma * I[t-1]
  dR = true_gamma * I[t-1]
  
  # Update S, I, R
  S[t] = S[t-1] + dS
  I[t] = I[t-1] + dI
  R[t] = R[t-1] + dR
}
'

#####################################

# Alternative approach: ODE RK45

sir = function(t, y_0, param) {
  with(as.list(c(y_0, param)), {
    dS = -beta * S * I / N
    dI = beta * S * I / N - gamma * I
    dR = -(dS + dI)
    list(c(dS, dI, dR))
  })
}

y_0 = c(S=S[1],I=I[1],R=R[1]+0.001)
param = c(beta=true_beta,gamma=true_gamma,N=n_sim)
ts=1:t_sim

bias = 0
msd = 1
n_tail = 0

data = ode(y_0, ts, sir, param, method = "ode45",rtol = 1e-6, atol = 1e-6, maxsteps = 1e6) %>% as_tibble() %>% rename(t=time) %>%
  mutate(across(everything(),~as.numeric(.x))) %>%
  mutate(across(!t,~.x+rnorm(.x, mean=0+bias, sd=msd))) %>%
  mutate(across(!t),round(.)) %>% 
  mutate(across(!t,~pmax(.x,1))) %>%
  slice_tail(n = t_sim-n_tail)

#data[1,] = list(1,S[1],I[1],R[1])


#####################################

# Measurement noise
'bias = 0
msd = 1

data = tibble(t=1:t_sim,S,I,R) %>%
#  mutate(across(!t,~.x+rnorm(.x, mean=0+bias, sd=msd))) %>%
#  mutate(across(!t),round(.)) %>%
  mutate(across(!t,~pmax(.x,1)))'

data %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = I, color = "Infectious")) +
  geom_line(aes(y = R, color = "Recovered")) +
  labs(x = "Time", y = "Count", color = "Category") +
  scale_color_manual(values = c("Susceptible" = "blue", "Infectious" = "red", "Recovered" = "green")) +
  theme_minimal()



'observations = tibble(y = rnorm(n_sim, mean = true_alpha + bias, sd = 0.01))
obs_mean = mean(observations$y)

observations %>% ggplot() +
  geom_histogram(aes(x = y), bins = 20, fill = "skyblue", color = "black") +
  geom_vline(xintercept=obs_mean, linetype="dashed", color = "darkblue", linewidth=0.5) +
  labs(title = "Simulated Data Distribution", x = "Observation", y = "Frequency")'


# Step 3: Draw from posterior distribution

meta_data = list(
  N = n_sim,
  T = t_sim-n_tail # t_sim
)

vector_data = data %>% select(I) %>% as.list

stan_data = c(meta_data,vector_data)

fitted = model$sample(
  data = stan_data,
  #init = list(init_data),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 200,
  refresh = 100,
  output_dir = getwd() %.% "/" %.% path,
  output_basename = basename,
  #  save_warmup = T,
  #  sig_figs = 18,
  seed = 1457L
)

posterior = fitted$draws(format = "df")

extracted = fitted %>% recover_types(stan_data) %>% spread_draws(states[T,sir] | sir) %>% rename(S=`1`,I=`2`,R=`3`) %>% mutate(T = T+1+n_tail)
add = data %>% rename(T=t) %>% mutate(.iteration=0,.chain=0,.draw=0)

extracted %>% filter(.chain == 1) %>% bind_rows(add) %>% mutate(obs = as.character(.iteration==0)) %>%
  ggplot(aes(x = T,group=.iteration,size=obs,alpha=obs)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = I, color = "Infectious")) +
  geom_line(aes(y = R, color = "Recovered")) +
  labs(x = "Time", y = "Count", color = "Category") +
  scale_color_manual(values = c("Susceptible" = "blue", "Infectious" = "red", "Recovered" = "green")) +
  scale_alpha_manual(values=c("TRUE"=1,"FALSE"=.01)) +
  scale_size_manual(values=c("TRUE"=1.2,"FALSE"=1)) +
  theme_minimal()




test = fitted %>% recover_types(stan_data) %>% spread_draws(gamma, beta, pred_cases[T], states_unbound[T,sir] | sir) %>% rename(S=`1`,I=`2`,R=`3`) %>%
  mutate(T = T+1) %>%
  filter(gamma < true_gamma+0.00005, gamma > true_gamma-0.00005, beta < true_beta+0.0001, beta > true_beta-0.0001)

add = data %>% rename(T=t) %>% mutate(.iteration=0,.chain=0,.draw=0,pred_cases=0)
test %>% bind_rows(add) %>% mutate(obs = as.character(.iteration==0)) %>%
  ggplot(aes(x = T,group=.iteration,linewidth=obs,alpha=obs,linetype=obs)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = I, color = "Infectious")) +
  geom_line(aes(y = R, color = "Recovered")) +
  geom_line(aes(y = pred_cases, color = "Predicted")) +
  labs(x = "Time", y = "Count", color = "Category") +
  scale_color_manual(values = c("Susceptible" = "blue", "Infectious" = "red", "Recovered" = "green","Predicted"="grey")) +
  scale_alpha_manual(values=c("TRUE"=1,"FALSE"=1)) +
  scale_linewidth_manual(values=c("TRUE"=1.2,"FALSE"=1)) +
  scale_linetype_manual(values=c("TRUE"="longdash","FALSE"="solid")) +
  theme_minimal()






ggplot() +
  geom_step(aes(x=t, y = S, color = "Susceptible"),data=data) +
  geom_step(aes(x=t, y = S, color = "Susceptible"),data=data,direction = "vh") +
  geom_step(aes(x=t, y = I, color = "Infectious"),data=data) +
  geom_step(aes(x=t, y = I, color = "Infectious"),data=data,direction = "vh") +
  geom_step(aes(x=t, y = R, color = "Recovered"),data=data) +
  geom_step(aes(x=t, y = R, color = "Recovered"),data=data,direction = "vh") +
  geom_line(aes(x=T, y = S, color = "Susceptible"),data=test) +
  geom_line(aes(x=T, y = I, color = "Infectious"),data=test) +
  geom_line(aes(x=T, y = R, color = "Recovered"),data=test) +
  geom_point(aes(x=t, y=I, color = "Simulated"),data=data_ode) +
  geom_line(aes(x=T, y = pred_cases, color = "Predicted"),data=test) +
  scale_color_manual(values = c("Susceptible" = "blue", "Infectious" = "red", "Recovered" = "green","Predicted"="grey", "Simulated"="violet")) +
  theme_minimal()






x=colnames(posterior) %>% as_tibble %>% filter(!startsWith(value,"states")) %>% filter(!startsWith(value,"pred"))

posterior %>% select(beta,gamma,S_0,I_0,R_0) %>% slice(1:200) %>% ggpairs(diag = list(continuous=wrap("barDiag",bins=40)))


# Step 4: Analyse results

fitted
fitted$cmdstan_summary()

summary = posterior %>% summarise_draws(prob=mean, MCSE=mcse_mean)
summary_chains = tibble(chain = 1:nchains) %>%
  mutate(summary = 
           map(chain, ~ posterior %>%
                subset_draws(chain = .x) %>%
                summarise_draws(prob = mean, MCSE = mcse_mean)
           ),
         pmean = map_dbl(summary, ~ .x %>% filter(variable == "alpha") %$% prob)
  )
posterior_mean = mean(summary_chains$pmean)


# Step 5: Visualise


posterior %>% ggplot(aes(x = .iteration, y = alpha, color = -lp__)) +
  geom_line() +
  geom_point() +
  labs(x = "Iteration", y = "Alpha", color = "Log probability") +
  facet_wrap(~ factor("Chain " %.% .chain), scales = "free", ncol=1) +
  scale_color_gradient(low = "grey", high = "red",trans = "log") +
  geom_hline(yintercept=posterior_mean, linetype="dotted", color = "black", linewidth=0.5) +
  theme_minimal()

posterior_wu %>% ggplot(aes(x = .iteration, y = alpha, color = -lp__)) +
  geom_line() +
  geom_point() +
  labs(x = "Iteration", y = "Alpha", color = "Log probability") +
  facet_wrap(~ factor("Chain " %.% .chain), scales = "free", ncol=1) +
  scale_color_gradient(low = "grey", high = "red",trans = "log") +
  geom_hline(yintercept=obs_mean, linetype="dashed", color = "darkblue", linewidth=0.5) +
  geom_hline(yintercept=posterior_mean, linetype="dotted", color = "black", linewidth=0.5) +
  theme_minimal()

posterior %>% ggplot(aes(x = lp__, y = alpha, color = factor(.chain))) +
  geom_line() +
  labs(x = "Log probability", y = "Alpha", color = "Chain") +
  facet_wrap(~ .chain, scales = "free", ncol=1) +
  geom_hline(yintercept=posterior_mean, linetype="dotted", color = "black", linewidth=0.5) +
  theme_minimal()

posterior_wu %>% ggplot(aes(x = lp__, y = alpha, color = factor(.chain))) +
  geom_line() +
  labs(x = "Log probability", y = "Alpha", color = "Chain") +
  facet_wrap(~ .chain, scales = "free", ncol=1) +
  geom_hline(yintercept=obs_mean, linetype="dashed", color = "darkblue", linewidth=0.5) +
  geom_hline(yintercept=posterior_mean, linetype="dotted", color = "black", linewidth=0.5) +
  theme_minimal()

posterior %>% ggplot(aes(x = alpha, y = lp__, color = factor(.chain))) +
  geom_line() +
  labs(x = "Alpha", y = "Log probability", color = "Chain") +
  facet_wrap(~ .chain, scales = "free", ncol=1) +
  geom_vline(xintercept=posterior_mean, linetype="dotted", color = "black", linewidth=0.5) +
  theme_minimal()

posterior_wu %>% ggplot(aes(x = alpha, y = lp__, color = factor(.chain))) +
  geom_line() +
  labs(x = "Alpha", y = "Log probability", color = "Chain") +
  facet_wrap(~ .chain, scales = "free", ncol=1) +
  geom_vline(xintercept=obs_mean, linetype="dashed", color = "darkblue", linewidth=0.5) +
  geom_vline(xintercept=posterior_mean, linetype="dotted", color = "black", linewidth=0.5) +
  theme_minimal()

mcmc_areas(fitted$draws("alpha"))
mcmc_hist(fitted$draws("alpha"))


#######

# Checks

c1 = posterior %$% mean(alpha)
c2 = posterior %$% weighted.mean(alpha,w=lp__)
c3 = posterior %$% weighted.mean(alpha,w=exp(lp__))
c4 = posterior %$% weighted.mean(alpha,w=exp(lp__)-max(lp__))
