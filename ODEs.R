### Model 1: susceptible-colonized model (MS equation 1)
ODEs_model1 <- function(t, states, params_pre){
  
  ## update parameter values:
  # (1/3): extract parameter values for compound parameters
  alpha = params_pre['alpha'];
  gamma = params_pre['gamma']; c = params_pre['c'];
  a = params_pre['a']; r_R = params_pre['r_R']; theta_C = params_pre['theta_C'];
  
  # (2/3): calculate compound parameter values
  alpha_R = as.numeric(alpha*(1-a*(1-r_R))); # endogenous acquisition only when not exposed to effective antibiotic therapy
  gamma_R = as.numeric(gamma*(1+c)); # higher clearance rate in C_R due to cost of resistance
  sigma_R = as.numeric(a*(1-r_R)*theta_C); # antibiotic pathogen clearance
  
  # (3/3): update parameter list
  params <- params_pre;
  params['alpha_R'] <- alpha_R;
  params['gamma_R'] <- gamma_R;
  params['sigma_R'] <- sigma_R;
  
  ### write derivatives
  with(as.list(c(states, params)),{
    
    dS <- (1-(f_C*f_R))*mu - (S*mu) - S*(beta*C_R + alpha_R) + C_R*(gamma_R + sigma_R)
    
    dC_R <- (f_C*f_R)*mu - (C_R*mu) + S*(beta*C_R + alpha_R) - C_R*(gamma_R + sigma_R)
    
    # incidence
    
    dC_S_trans <- 0
    
    dC_S_acq <- 0
    
    dC_R_trans <- beta*(C_R)*S
    
    dC_R_acq <- alpha_R*(S)
    
    dC_R_hgt <- 0
    
    # return derivatives
    der <- c(dS, dC_R, dC_S_trans, dC_S_acq, dC_R_trans, dC_R_acq, dC_R_hgt)
    list(der)
  }) 
}


### Model 2: strain competition model (MS equation 3)
ODEs_model2 <- function(t, states, params_pre){
  
  ## update parameter values:
  # (1/3): extract parameter values needed for compound parameters
  alpha = params_pre['alpha'];
  gamma = params_pre['gamma']; c = params_pre['c'];
  a = params_pre['a']; r_S = params_pre['r_S'];  r_R = params_pre['r_R']; theta_C = params_pre['theta_C'];
  
  # (2/3): calculate compound parameter values
  alpha_S = as.numeric(alpha*(1-a*(1-r_S)));# endogenous acquisition only when not exposed to effective antibiotic therapy
  alpha_R = as.numeric(alpha*(1-a*(1-r_R)));
  gamma_S = as.numeric(gamma); 
  gamma_R = as.numeric(gamma*(1+c)); # higher clearance rate in C_R due to cost of resistance
  sigma_S = as.numeric(a*(1-r_S)*theta_C); # antibiotic pathogen clearance
  sigma_R = as.numeric(a*(1-r_R)*theta_C);
  
  # (3/3): update parameter list
  params <- params_pre
  params['alpha_S'] <- alpha_S; params['alpha_R'] <- alpha_R;
  params['gamma_S'] <- gamma_S; params['gamma_R'] <- gamma_R;
  params['sigma_S'] <- sigma_S; params['sigma_R'] <- sigma_R;
  
  ### write derivatives
  with(as.list(c(states, params)),{
    
    dS <- (1-f_C)*mu - (S*mu) - S*(beta*(C_S+C_R) + alpha_S + alpha_R) + C_S*(gamma_S + sigma_S) + C_R*(gamma_R + sigma_R)
    
    dC_S <- f_C*(1-f_R)*mu - (C_S*mu) + S*(beta*C_S + alpha_S) - C_S*(gamma_S + sigma_S)
    
    dC_R <- f_C*f_R*mu - (C_R*mu) + S*(beta*C_R + alpha_R) - C_R*(gamma_R + sigma_R)
    
    # incidence
    dC_S_trans <- beta*(C_S)*S
    
    dC_S_acq <- alpha_S*(S)
    
    dC_R_trans <- beta*(C_R)*S
    
    dC_R_acq <- alpha_R*(S)
    
    dC_R_hgt <- 0
    
    # return derivatives
    der <- c(dS, dC_S, dC_R,
             dC_S_trans, dC_S_acq, dC_R_trans, dC_R_acq, dC_R_hgt)
    list(der)
  }) 
}


### Model 3: microbiome competition model (MS equation 4)
ODEs_model3 <- function(t, states, params_pre){
  
  ## update parameter values:
  # (1/3): extract parameter values needed for compound parameters
  beta = params_pre['beta']; epsilon = params_pre['epsilon']
  alpha = params_pre['alpha']; phi = params_pre['phi'];
  gamma = params_pre['gamma']; c = params_pre['c']; eta = params_pre['eta'];
  a = params_pre['a']; r_R = params_pre['r_R']; 
  theta_C = params_pre['theta_C']; theta_m = params_pre['theta_m']; delta = params_pre['delta'];
  
  # (2/3): calculate compound parameter values
  beta_epsilon = as.numeric(beta*(1-epsilon)) # transmission rate reduced by colonization resistance
  alpha_R = as.numeric(alpha*(1-a*(1-r_R))); # endogenous acquisition only when not exposed to effective antibiotic therapy
  alpha_R_phi = as.numeric(alpha_R*phi); # endogenous acquisition augmented by ecological release
  gamma_R = as.numeric(gamma*(1+c)); # higher clearance rate in C_R due to cost of resistance
  gamma_R_eta = as.numeric(gamma_R*(1-eta)) # clearance rate lowered during dysbiosis (less resource competition)
  sigma_R = as.numeric(a*(1-r_R)*theta_C); # antibiotic pathogen clearance
  sigma_m = as.numeric(a*theta_m)  # antibiotic microbiome dysbiosis
  delta = as.numeric(delta*(1-a)) # microbiome recovery
  
  # (3/3): update parameter list
  params <- params_pre
  params['beta_epsilon'] <- beta_epsilon;
  params['alpha_R'] <- alpha_R; params['alpha_R_phi'] <- alpha_R_phi;
  params['gamma_R'] <- gamma_R; params['gamma_R_eta'] <- gamma_R_eta;
  params['sigma_R'] <- sigma_R; params['sigma_m'] <- sigma_m;
  params['delta'] <- delta;
  
  ### write derivatives
  with(as.list(c(states, params)),{
    
    dS_e <- (1-(f_C*f_R))*(1-f_d)*mu - (S_e*mu) - S_e*(beta_epsilon*(C_R_e+C_R_d) + alpha_R + sigma_m) + S_d*delta + C_R_e*(gamma_R + sigma_R)
    
    dS_d <- (1-(f_C*f_R))*f_d*mu - (S_d*mu) + S_e*sigma_m - S_d*(beta*(C_R_e+C_R_d) + alpha_R_phi + delta) + C_R_d*(gamma_R_eta + sigma_R)
    
    dC_R_e <- f_C*f_R*(1-f_d)*mu - (C_R_e*mu) + S_e*(beta_epsilon*(C_R_e+C_R_d) + alpha_R) - C_R_e*(gamma_R + sigma_m + sigma_R) + (C_R_d*delta)
    
    dC_R_d <- f_C*f_R*f_d*mu - (C_R_d*mu) + S_d*(beta*(C_R_e + C_R_d) + alpha_R_phi) + C_R_e*sigma_m - C_R_d*(gamma_R_eta + delta + sigma_R)
    
    # incidence
    
    dC_S_trans <- 0
    
    dC_S_acq <- 0
    
    dC_R_trans <-       
      beta_epsilon*(C_R_e + C_R_d)*S_e + 
      beta*(C_R_e + C_R_d)*S_d
    
    dC_R_acq <- 
      alpha_R*(S_e)+
      alpha_R_phi*(S_d)
    
    dC_R_hgt <- 0
    
    # return derivatives
    der <- c(dS_e, dS_d, dC_R_e, dC_R_d,
             dC_S_trans, dC_S_acq, dC_R_trans, dC_R_acq, dC_R_hgt)
    list(der)
  }) 
}


### Model 4: two-strain microbiome competition model (MS equation 6)
ODEs_model4 <- function(t, states, params_pre){
  
  ## update parameter values:
  # (1/3): extract parameter values needed for compound parameters
  beta = params_pre['beta']; epsilon = params_pre['epsilon']
  alpha = params_pre['alpha']; phi = params_pre['phi'];
  gamma = params_pre['gamma']; c = params_pre['c']; eta = params_pre['eta'];
  a = params_pre['a']; r_S = params_pre['r_S']; r_R = params_pre['r_R']; 
  theta_C = params_pre['theta_C']; theta_m = params_pre['theta_m']; delta = params_pre['delta']
  
  # (2/3): calculate compound parameter values
  beta_epsilon = as.numeric(beta*(1-epsilon)) # transmission rate reduced by colonization resistance
  alpha_S = as.numeric(alpha*(1-a*(1-r_S))); # endogenous acquisition only when not exposed to effective antibiotic therapy
  alpha_R = as.numeric(alpha*(1-a*(1-r_R)));
  alpha_S_phi = as.numeric(alpha_S*phi); # endogenous acquisition augmented by ecological release
  alpha_R_phi = as.numeric(alpha_R*phi);
  gamma_S = as.numeric(gamma); 
  gamma_R = as.numeric(gamma*(1+c)); # higher clearance rate in C_R due to cost of resistance
  gamma_S_eta = as.numeric(gamma_S*(1-eta));  # clearance rate lowered during dysbiosis (less resource competition)
  gamma_R_eta = as.numeric(gamma_R*(1-eta));
  sigma_S = as.numeric(a*(1-r_S)*theta_C) # antibiotic pathogen clearance
  sigma_R = as.numeric(a*(1-r_R)*theta_C);
  sigma_m = as.numeric(a*theta_m)  # antibiotic microbiome dysbiosis
  delta = as.numeric(delta*(1-a)); # microbiome recovery
  
  # (3/3): update parameter list
  params <- params_pre
  params['beta_epsilon'] <- beta_epsilon;
  params['alpha_S'] <- alpha_S; params['alpha_R'] <- alpha_R;
  params['alpha_S_phi'] <- alpha_S_phi; params['alpha_R_phi'] <- alpha_R_phi;
  params['gamma_S'] <- gamma_S; params['gamma_R'] <- gamma_R;
  params['gamma_S_eta'] <- gamma_S_eta; params['gamma_R_eta'] <- gamma_R_eta;
  params['sigma_S'] <- sigma_S; params['sigma_R'] <- sigma_R; params['sigma_m'] <- sigma_m;
  params['delta'] <- delta;
  
  ### calculate derivatives
  with(as.list(c(states, params)),{
    
    dS_e <- (1-(f_C))*(1-f_d)*mu - (S_e*mu) - S_e*(beta_epsilon*(C_S_e+C_S_d+C_R_e+C_R_d) + alpha_S + alpha_R + sigma_m) + S_d*delta + C_S_e*(gamma_S + sigma_S) + C_R_e*(gamma_R + sigma_R)
    
    dS_d <- (1-(f_C))*f_d*mu - (S_d*mu) + S_e*sigma_m - S_d*(beta*(C_S_e+C_S_d+C_R_e+C_R_d) + alpha_S_phi + alpha_R_phi + delta) + C_S_d*(gamma_S_eta + sigma_S) + C_R_d*(gamma_R_eta + sigma_R)
    
    dC_S_e <- f_C*(1-f_R)*(1-f_d)*mu - (C_S_e*mu) + S_e*(beta_epsilon*(C_S_e+C_S_d) + alpha_S) - C_S_e*(gamma_S + sigma_m + sigma_S) + (C_S_d*delta)
    
    dC_S_d <- f_C*(1-f_R)*f_d*mu - (C_S_d*mu) + S_d*(beta*(C_S_e + C_S_d) + alpha_S_phi) + C_S_e*sigma_m - C_S_d*(gamma_S_eta + delta + sigma_S)
    
    dC_R_e <- f_C*f_R*(1-f_d)*mu - (C_R_e*mu) + S_e*(beta_epsilon*(C_R_e+C_R_d) + alpha_R) - C_R_e*(gamma_R + sigma_m + sigma_R) + (C_R_d*delta)
    
    dC_R_d <- f_C*f_R*f_d*mu - (C_R_d*mu) + S_d*(beta*(C_R_e + C_R_d) + alpha_R_phi) + C_R_e*sigma_m - C_R_d*(gamma_R_eta + delta + sigma_R)
    
    # incidence
    dC_S_trans <- 
      beta_epsilon*(C_S_e + C_S_d)*S_e + 
      beta*(C_S_e + C_S_d)*S_d
    
    dC_S_acq <- 
      alpha_S*(S_e) + 
      alpha_S_phi*(S_d)
    
    dC_R_trans <-       
      beta_epsilon*(C_R_e + C_R_d)*S_e + 
      beta*(C_R_e + C_R_d)*S_d
    
    dC_R_acq <- 
      alpha_R*(S_e)+
      alpha_R_phi*(S_d)
    
    dC_R_hgt <- 0
    
    # return derivatives
    der <- c(dS_e, dS_d, dC_S_e, dC_S_d, dC_R_e, dC_R_d,
             dC_S_trans, dC_S_acq, dC_R_trans, dC_R_acq, dC_R_hgt)
    list(der)
  }) 
}



### Model 5: two-strain microbiome competition model with HGT (MS equation XXXXXX)
ODEs_model5 <- function(t, states, params_pre){
  
  ## update parameter values:
  # (1/3): extract parameter values needed for compound parameters
  beta = params_pre['beta']; epsilon = params_pre['epsilon']
  alpha = params_pre['alpha']; phi = params_pre['phi'];
  gamma = params_pre['gamma']; c = params_pre['c']; eta = params_pre['eta'];
  a = params_pre['a']; r_S = params_pre['r_S']; r_R = params_pre['r_R']; 
  theta_C = params_pre['theta_C']; theta_m = params_pre['theta_m']; delta = params_pre['delta']
  
  # (2/3): calculate compound parameter values
  beta_epsilon = as.numeric(beta*(1-epsilon)) # transmission rate reduced by colonization resistance
  alpha_S = as.numeric(alpha*(1-a*(1-r_S))); # endogenous acquisition only when not exposed to effective antibiotic therapy
  alpha_R = as.numeric(alpha*(1-a*(1-r_R))); 
  alpha_S_phi = as.numeric(alpha_S*phi); # endogenous acquisition augmented by ecological release
  alpha_R_phi = as.numeric(alpha_R*phi); 
  gamma_S = as.numeric(gamma);
  gamma_R = as.numeric(gamma*(1+c)); # higher clearance rate in C_R due to cost of resistance
  gamma_S_eta = as.numeric(gamma_S*(1-eta)); # clearance rate lowered during dysbiosis (less resource competition)
  gamma_R_eta = as.numeric(gamma_R*(1-eta));
  sigma_S = as.numeric(a*(1-r_S)*theta_C) # antibiotic pathogen clearance
  sigma_R = as.numeric(a*(1-r_R)*theta_C);
  sigma_m = as.numeric(a*theta_m) # antibiotic microbiome dysbiosis
  delta = as.numeric(delta*(1-a)); # microbiome recovery
  
  # (3/3): update parameter list
  params <- params_pre
  params['beta_epsilon'] <- beta_epsilon;
  params['alpha_S'] <- alpha_S; params['alpha_R'] <- alpha_R;
  params['alpha_S_phi'] <- alpha_S_phi; params['alpha_R_phi'] <- alpha_R_phi;
  params['gamma_S'] <- gamma_S; params['gamma_R'] <- gamma_R;
  params['gamma_S_eta'] <- gamma_S_eta; params['gamma_R_eta'] <- gamma_R_eta;
  params['sigma_S'] <- sigma_S; params['sigma_R'] <- sigma_R; params['sigma_m'] <- sigma_m;
  params['delta'] <- delta;
  
  ### derivatives
  with(as.list(c(states, params)),{
    
    dS_e_s <- (1-(f_C))*(1-f_d)*(1-f_w)*mu - (S_e_s*mu) - S_e_s*(beta_epsilon*(C_S_e_s+C_S_e_r+C_S_d_s+C_S_d_r+C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_S + alpha_R + sigma_m) + S_d_s*delta + C_S_e_s*(gamma_S + sigma_S) + C_R_e_s*(gamma_R + sigma_R)
    
    dS_e_r <- (1-(f_C))*(1-f_d)*f_w*mu - (S_e_r*mu) - S_e_r*(beta_epsilon*(C_S_e_s+C_S_e_r+C_S_d_s+C_S_d_r+C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_S + alpha_R + sigma_m) + S_d_r*delta + C_S_e_r*(gamma_S + sigma_S) + C_R_e_r*(gamma_R + sigma_R)
    
    dS_d_s <- (1-(f_C))*f_d*(1-f_w)*mu - (S_d_s*mu) + S_e_s*((1-omega)*sigma_m) - S_d_s*(beta*(C_S_e_s+C_S_e_r+C_S_d_s+C_S_d_r+C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_S_phi + alpha_R_phi + delta) + C_S_d_s*(gamma_S_eta + sigma_S) + C_R_d_s*(gamma_R_eta + sigma_R)
    
    dS_d_r <- (1-(f_C))*f_d*f_w*mu - (S_d_r*mu) + S_e_s*(omega*sigma_m) + S_e_r*sigma_m - S_d_r*(beta*(C_S_e_s+C_S_e_r+C_S_d_s+C_S_d_r+C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_S_phi + alpha_R_phi + delta) + C_S_d_r*(gamma_S_eta + sigma_S) + C_R_d_r*(gamma_R_eta + sigma_R)
    
    dC_S_e_s <- f_C*(1-f_R)*(1-f_d)*(1-f_w)*mu - (C_S_e_s*mu) + S_e_s*(beta_epsilon*(C_S_e_s+C_S_e_r+C_S_d_s+C_S_d_r) + alpha_S) - C_S_e_s*(gamma_S + sigma_m + sigma_S) + C_S_d_s*delta
    
    dC_S_e_r <- f_C*(1-f_R)*(1-f_d)*f_w*mu - (C_S_e_r*mu) + S_e_r*(beta_epsilon*(C_S_e_s+C_S_e_r+C_S_d_s+C_S_d_r) + alpha_S) - C_S_e_r*(gamma_S + sigma_m + sigma_S + chi_e) + C_S_d_r*delta
    
    dC_S_d_s <- f_C*(1-f_R)*f_d*(1-f_w)*mu - (C_S_d_s*mu) + S_d_s*(beta*(C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r) + alpha_S_phi) + C_S_e_s*((1-omega)*sigma_m) - C_S_d_s*(gamma_S_eta + delta + sigma_S)
    
    dC_S_d_r <- f_C*(1-f_R)*f_d*f_w*mu - (C_S_d_r*mu) + S_d_r*(beta*(C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r) + alpha_S_phi) + C_S_e_s*(omega*sigma_m) + C_S_e_r*sigma_m - C_S_d_r*(gamma_S_eta + delta + sigma_S + chi_d)
    
    dC_R_e_s <- f_C*f_R*(1-f_d)*(1-f_w)*mu - (C_R_e_s*mu) + S_e_s*(beta_epsilon*(C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_R) - C_R_e_s*(gamma_R + sigma_m + sigma_R + chi_e) + (C_R_d_s*delta)
    
    dC_R_e_r <- f_C*f_R*(1-f_d)*f_w*mu - (C_R_e_r*mu) + S_e_r*(beta_epsilon*(C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_R) + C_S_e_r*chi_e + C_R_e_s*chi_e - C_R_e_r*(gamma_R + sigma_m + sigma_R) + C_R_d_r*delta
    
    dC_R_d_s <- f_C*f_R*f_d*(1-f_w)*mu - (C_R_d_s*mu) + S_d_s*(beta*(C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_R_phi) + C_R_e_s*((1-omega)*sigma_m) - C_R_d_s*(gamma_R_eta + delta + sigma_R + chi_d)
    
    dC_R_d_r <- f_C*f_R*f_d*f_w*mu - (C_R_d_r*mu) + S_d_r*(beta*(C_R_e_s+C_R_e_r+C_R_d_s+C_R_d_r) + alpha_R_phi) + C_S_d_r*chi_d + C_R_e_s*(omega*sigma_m) + C_R_e_r*sigma_m + C_R_d_s*chi_d - C_R_d_r*(gamma_R_eta + delta + sigma_R)
    
    # incidence
    dC_S_trans <- 
      beta_epsilon*(C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r)*S_e_s + 
      beta_epsilon*(C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r)*S_e_r + 
      beta*(C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r)*S_d_s + 
      beta*(C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r)*S_d_r
    
    dC_S_acq <- 
      alpha_S*(S_e_s + S_e_r) + 
      alpha_S_phi*(S_d_s + S_d_r)
    
    dC_R_trans <-       
      beta_epsilon*(C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r)*S_e_s + 
      beta_epsilon*(C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r)*S_e_r + 
      beta*(C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r)*S_d_s + 
      beta*(C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r)*S_d_r
    
    dC_R_acq <- 
      alpha_R*(S_e_s + S_e_r)+
      alpha_R_phi*(S_d_s + S_d_r)
    
    dC_R_hgt <-
      (chi_e)*(C_S_e_r) +
      (chi_d)*(C_S_d_r)
    
    
    # return derivatives
    der <- c(dS_e_s, dS_e_r, dS_d_s, dS_d_r, dC_S_e_s, dC_S_e_r, dC_S_d_s, dC_S_d_r, dC_R_e_s, dC_R_e_r, dC_R_d_s, dC_R_d_r,
             dC_S_trans, dC_S_acq, dC_R_trans, dC_R_acq, dC_R_hgt)
    list(der)
  }) 
}
