##############
### STATES ###
##############

### DEFAULT STATES
### NB: model-specific
states_model1 <- c(S = 0.9, C_R = 0.1, 
                   dC_S_trans = 0, dC_S_acq = 0, dC_R_trans = 0, dC_R_acq = 0, dC_R_hgt = 0)
states_model2 <- c(S = 0.8, C_S = 0.1, C_R = 0.1, 
                   dC_S_trans = 0, dC_S_acq = 0, dC_R_trans = 0, dC_R_acq = 0, dC_R_hgt = 0)
states_model3 <- c(S_e = 0.7, S_d = 0.1, C_R_e = 0.1, C_R_d = 0.1,
                   dC_S_trans = 0, dC_S_acq = 0, dC_R_trans = 0, dC_R_acq = 0, dC_R_hgt = 0)
states_model4 <- c(S_e = 0.7, S_d = 0.1, C_S_e = 0.05, C_S_d = 0.05, C_R_e = 0.05, C_R_d = 0.05,
                   dC_S_trans = 0, dC_S_acq = 0, dC_R_trans = 0, dC_R_acq = 0, dC_R_hgt = 0)
states_model5 <- c(S_e_s = 0.3, S_e_r = 0.03, S_d_s = 0.2, S_d_r = 0.02, C_S_e_s = 0.1, C_S_e_r = 0.05, C_S_d_s = 0.05, C_S_d_r = 0.05, C_R_e_s = 0.05, C_R_e_r = 0.05, C_R_d_s = 0.05, C_R_d_r = 0.05,
                   dC_S_trans = 0, dC_S_acq = 0, dC_R_trans = 0, dC_R_acq = 0, dC_R_hgt = 0)
 
##################
### PARAMETERS ###
##################
### NB: not model-specific

### DEFAULT PARAMETERS
params_default <- c(
  # pathogen colonization
  beta = 0.2,
  alpha = 0.01,
  gamma = 0.03,
  c = 1,
  # patient demography
  mu = 0.1,
  f_d = 0,
  f_C = 0.1,#0.02,
  f_R = 0.5,
  f_w = 0,
  # antibiotics
  a = 0.2,
  r_S = 0,
  r_R = 0.8,
  theta_C = 0.2,
  theta_m = 1,
  # microbiome ecology
  epsilon = 0.5,
  eta = 0.5,
  phi = 5,
  chi_e = 0.01,
  chi_d = 0.1,
  delta = 1/7,
  omega = 0.01
)

##############################
### ALTERNATIVE PARAMETERS ###
##############################

### FIGURE 1
params_perfectR = params_default
params_perfectR['r_R'] <- 1

### variable epsilon (colonization resistance)
# defaults
params_default_epsilon_low = params_default
params_default_epsilon_med = params_default
params_default_epsilon_high = params_default
params_default_epsilon_low['epsilon'] <- 0.2; params_default_epsilon_low['eta'] <- 0 ; params_default_epsilon_low['phi'] <- 1
params_default_epsilon_med['epsilon'] <- 0.5; params_default_epsilon_med['eta'] <- 0 ; params_default_epsilon_med['phi'] <- 1
params_default_epsilon_high['epsilon'] <- 0.8; params_default_epsilon_high['eta'] <- 0 ; params_default_epsilon_high['phi'] <- 1

# perfect r_R
params_perfectR_epsilon_low = params_perfectR
params_perfectR_epsilon_med = params_perfectR
params_perfectR_epsilon_high = params_perfectR
params_perfectR_epsilon_low['epsilon'] <- 0.2; params_perfectR_epsilon_low['eta'] <- 0 ; params_perfectR_epsilon_low['phi'] <- 1
params_perfectR_epsilon_med['epsilon'] <- 0.5; params_perfectR_epsilon_med['eta'] <- 0 ; params_perfectR_epsilon_med['phi'] <- 1
params_perfectR_epsilon_high['epsilon'] <- 0.8; params_perfectR_epsilon_high['eta'] <- 0 ; params_perfectR_epsilon_high['phi'] <- 1

### variable eta (resource competition)
# defaults
params_default_eta_low = params_default
params_default_eta_med = params_default
params_default_eta_high = params_default
params_default_eta_low['epsilon'] <- 0; params_default_eta_low['eta'] <- 0.2 ; params_default_eta_low['phi'] <- 1
params_default_eta_med['epsilon'] <- 0; params_default_eta_med['eta'] <- 0.5 ; params_default_eta_med['phi'] <- 1
params_default_eta_high['epsilon'] <- 0; params_default_eta_high['eta'] <- 0.8 ; params_default_eta_high['phi'] <- 1

# perfect r_R
params_perfectR_eta_low = params_perfectR
params_perfectR_eta_med = params_perfectR
params_perfectR_eta_high = params_perfectR
params_perfectR_eta_low['epsilon'] <- 0; params_perfectR_eta_low['eta'] <- 0.2 ; params_perfectR_eta_low['phi'] <- 1
params_perfectR_eta_med['epsilon'] <- 0; params_perfectR_eta_med['eta'] <- 0.5 ; params_perfectR_eta_med['phi'] <- 1
params_perfectR_eta_high['epsilon'] <- 0; params_perfectR_eta_high['eta'] <- 0.8 ; params_perfectR_eta_high['phi'] <- 1

### variable phi (ecological release)
# defaults
params_default_phi_low = params_default
params_default_phi_med = params_default
params_default_phi_high = params_default
params_default_phi_low['epsilon'] <- 0; params_default_phi_low['eta'] <- 0 ; params_default_phi_low['phi'] <- 2
params_default_phi_med['epsilon'] <- 0; params_default_phi_med['eta'] <- 0 ; params_default_phi_med['phi'] <- 5
params_default_phi_high['epsilon'] <- 0; params_default_phi_high['eta'] <- 0 ; params_default_phi_high['phi'] <- 8

# perfect r_R
params_perfectR_phi_low = params_perfectR
params_perfectR_phi_med = params_perfectR
params_perfectR_phi_high = params_perfectR
params_perfectR_phi_low['epsilon'] <- 0; params_perfectR_phi_low['eta'] <- 0 ; params_perfectR_phi_low['phi'] <- 2
params_perfectR_phi_med['epsilon'] <- 0; params_perfectR_phi_med['eta'] <- 0 ; params_perfectR_phi_med['phi'] <- 5
params_perfectR_phi_high['epsilon'] <- 0; params_perfectR_phi_high['eta'] <- 0 ; params_perfectR_phi_high['phi'] <- 8


### FIGURE 3
params_lowInt_lowR = params_default; params_lowInt_medR = params_default; params_lowInt_highR = params_default; 
params_medInt_lowR = params_default; params_medInt_medR = params_default; params_medInt_highR = params_default; 
params_highInt_lowR = params_default; params_highInt_medR = params_default; params_highInt_highR = params_default; 

params_lowInt_lowR['epsilon'] = 0.2; params_lowInt_lowR['eta'] = 0.2; params_lowInt_lowR['phi'] = 2; params_lowInt_lowR['r_R'] = 0.2
params_lowInt_medR['epsilon'] = 0.2; params_lowInt_medR['eta'] = 0.2; params_lowInt_medR['phi'] = 2; params_lowInt_medR['r_R'] = 0.5
params_lowInt_highR['epsilon'] = 0.2; params_lowInt_highR['eta'] = 0.2; params_lowInt_highR['phi'] = 2; params_lowInt_highR['r_R'] = 0.8
params_medInt_lowR['epsilon'] = 0.5; params_medInt_lowR['eta'] = 0.5; params_medInt_lowR['phi'] = 5; params_medInt_lowR['r_R'] = 0.2
params_medInt_medR['epsilon'] = 0.5; params_medInt_medR['eta'] = 0.5; params_medInt_medR['phi'] = 5; params_medInt_medR['r_R'] = 0.5
params_medInt_highR['epsilon'] = 0.5; params_medInt_highR['eta'] = 0.5; params_medInt_highR['phi'] = 5; params_medInt_highR['r_R'] = 0.8
params_highInt_lowR['epsilon'] = 0.8; params_highInt_lowR['eta'] = 0.8; params_highInt_lowR['phi'] = 8; params_highInt_lowR['r_R'] = 0.2
params_highInt_medR['epsilon'] = 0.8; params_highInt_medR['eta'] = 0.8; params_highInt_medR['phi'] = 8; params_highInt_medR['r_R'] = 0.5
params_highInt_highR['epsilon'] = 0.8; params_highInt_highR['eta'] = 0.8; params_highInt_highR['phi'] = 8; params_highInt_highR['r_R'] = 0.8

### FIGURE 4
# Default pars
params_noInt_noHGT = params_default; params_noInt_lowHGT = params_default; params_noInt_highHGT = params_default; 
params_withInt_noHGT = params_default; params_withInt_lowHGT = params_default; params_withInt_highHGT = params_default

params_noInt_noHGT['epsilon'] = 0; params_noInt_noHGT['eta'] = 0; params_noInt_noHGT['phi'] = 1; params_noInt_noHGT['chi_e'] = 0; params_noInt_noHGT['chi_d'] = 0
params_noInt_lowHGT['epsilon'] = 0; params_noInt_lowHGT['eta'] = 0; params_noInt_lowHGT['phi'] = 1; params_noInt_lowHGT['chi_e'] = 0.01; params_noInt_lowHGT['chi_d'] = 0.1
params_noInt_highHGT['epsilon'] = 0; params_noInt_highHGT['eta'] = 0; params_noInt_highHGT['phi'] = 1; params_noInt_highHGT['chi_e'] = 0.1; params_noInt_highHGT['chi_d'] = 1
params_withInt_noHGT['epsilon'] = 0.5; params_withInt_noHGT['eta'] = 0.5; params_withInt_noHGT['phi'] = 5; params_withInt_noHGT['chi_e'] = 0; params_withInt_noHGT['chi_d'] = 0
params_withInt_lowHGT['epsilon'] = 0.5; params_withInt_lowHGT['eta'] = 0.5; params_withInt_lowHGT['phi'] = 5; params_withInt_lowHGT['chi_e'] = 0.01; params_withInt_lowHGT['chi_d'] = 0.1
params_withInt_highHGT['epsilon'] = 0.5; params_withInt_highHGT['eta'] = 0.5; params_withInt_highHGT['phi'] = 5; params_withInt_highHGT['chi_e'] = 0.1; params_withInt_highHGT['chi_d'] = 1

# Medium resistance (r_R = 0.2)
params_noInt_noHGT_medR = params_noInt_noHGT; params_noInt_lowHGT_medR = params_noInt_lowHGT; params_noInt_highHGT_medR = params_noInt_highHGT; 
params_withInt_noHGT_medR = params_withInt_noHGT; params_withInt_lowHGT_medR = params_withInt_lowHGT; params_withInt_highHGT_medR = params_withInt_highHGT

params_noInt_noHGT_medR['r_R'] = 0.5; 
params_noInt_lowHGT_medR['r_R'] = 0.5; 
params_noInt_highHGT_medR['r_R'] = 0.5; 
params_withInt_noHGT_medR['r_R'] = 0.5;
params_withInt_lowHGT_medR['r_R'] = 0.5;
params_withInt_highHGT_medR['r_R'] = 0.5

# Low resistance (r_R = 0.2)
params_noInt_noHGT_lowR = params_noInt_noHGT; params_noInt_lowHGT_lowR = params_noInt_lowHGT; params_noInt_highHGT_lowR = params_noInt_highHGT; 
params_withInt_noHGT_lowR = params_withInt_noHGT; params_withInt_lowHGT_lowR = params_withInt_lowHGT; params_withInt_highHGT_lowR = params_withInt_highHGT

params_noInt_noHGT_lowR['r_R'] = 0.2; 
params_noInt_lowHGT_lowR['r_R'] = 0.2; 
params_noInt_highHGT_lowR['r_R'] = 0.2; 
params_withInt_noHGT_lowR['r_R'] = 0.2;
params_withInt_lowHGT_lowR['r_R'] = 0.2;
params_withInt_highHGT_lowR['r_R'] = 0.2


### HGT supplementary Figure S7
# a
params_HGTsupp_noHGT = params_default; params_HGTsupp_noHGT['chi_e'] = 0; params_HGTsupp_noHGT['chi_d'] = 0; 
params_HGTsupp_lowHGT = params_default; params_HGTsupp_lowHGT['chi_e'] = 0.01; params_HGTsupp_lowHGT['chi_d'] = 0.1; 
params_HGTsupp_highHGT = params_default; params_HGTsupp_highHGT['chi_e'] = 0.1; params_HGTsupp_highHGT['chi_d'] = 1
# b
params_HGTsupp_varyHGT1 = params_default; params_HGTsupp_varyHGT1['epsilon'] = 0; params_HGTsupp_varyHGT1['eta'] = 0; params_HGTsupp_varyHGT1['phi'] = 1; 
params_HGTsupp_varyHGT1['chi_e'] = 0.05; params_HGTsupp_varyHGT1['chi_d'] = 0.05; 
params_HGTsupp_varyHGT2 = params_HGTsupp_varyHGT1; params_HGTsupp_varyHGT2['chi_d'] = 0.05*2
params_HGTsupp_varyHGT3 = params_HGTsupp_varyHGT1; params_HGTsupp_varyHGT3['chi_d'] = 0.05*4
params_HGTsupp_varyHGT4 = params_HGTsupp_varyHGT1; params_HGTsupp_varyHGT4['chi_d'] = 0.05*8
params_HGTsupp_varyHGT5 = params_HGTsupp_varyHGT1; params_HGTsupp_varyHGT5['chi_d'] = 0.05*16
# c
params_HGTsupp_vary_c = params_default; params_HGTsupp_vary_c['epsilon'] = 0; params_HGTsupp_vary_c['eta'] = 0; params_HGTsupp_vary_c['phi'] = 1; 
params_HGTsupp_vary_c0 = params_HGTsupp_vary_c; params_HGTsupp_vary_c0['chi_e'] = 0; params_HGTsupp_vary_c0['chi_d'] = 0; 
params_HGTsupp_vary_c1 = params_HGTsupp_vary_c; params_HGTsupp_vary_c1['chi_e'] = 0.01; params_HGTsupp_vary_c1['chi_d'] = 0.05; 

params_HGTsupp_vary_c0_1 = params_HGTsupp_vary_c0; params_HGTsupp_vary_c0_1['c'] = -0.5
params_HGTsupp_vary_c1_1 = params_HGTsupp_vary_c1; params_HGTsupp_vary_c1_1['c'] = -0.5

params_HGTsupp_vary_c0_2 = params_HGTsupp_vary_c0; params_HGTsupp_vary_c0_2['c'] = 0
params_HGTsupp_vary_c1_2 = params_HGTsupp_vary_c1; params_HGTsupp_vary_c1_2['c'] = 0

params_HGTsupp_vary_c0_3 = params_HGTsupp_vary_c0; params_HGTsupp_vary_c0_3['c'] = 1
params_HGTsupp_vary_c1_3 = params_HGTsupp_vary_c1; params_HGTsupp_vary_c1_3['c'] = 1

params_HGTsupp_vary_c0_4 = params_HGTsupp_vary_c0; params_HGTsupp_vary_c0_4['c'] = 2
params_HGTsupp_vary_c1_4 = params_HGTsupp_vary_c1; params_HGTsupp_vary_c1_4['c'] = 2

params_HGTsupp_vary_c0_5 = params_HGTsupp_vary_c0; params_HGTsupp_vary_c0_5['c'] = 4
params_HGTsupp_vary_c1_5 = params_HGTsupp_vary_c1; params_HGTsupp_vary_c1_5['c'] = 4


### Dynamic responses to public health interventions

states_invasion <- c(S_e_s = 0.649, S_e_r = 0, S_d_s = 0.2, S_d_r = 0, C_S_e_s = 0.1, C_S_e_r = 0, C_S_d_s = 0.05, C_S_d_r = 0, C_R_e_s = 0.001, C_R_e_r = 0, C_R_d_s = 0, C_R_d_r = 0,
                   dC_S_trans = 0, dC_S_acq = 0, dC_R_trans = 0, dC_R_acq = 0, dC_R_hgt = 0)

params_invasion = params_default; 

par_intervention1 = 'r_R'
par_intervention2 = 'theta_m'
par_intervention3 = 'a'

val_intervention1 = 0.4
val_intervention2 = 0.5
val_intervention3 = 0.1







