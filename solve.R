library(tidyr)
library(dplyr)
library(magrittr)
library(deSolve)

source('ODEs.R')
source('functions.R')
source('parameters.R')


############
### PLAY ###
############
# examples of outputs from model and functions
# adjust parameters in param_play to evaluate impact on model behaviour

### DYNAMICS
# vary any of these parameters to evaluate impacts on model outputs:
params_play = params_default
params_play['r_R'] = 0.8
params_play['beta'] = 0.2
params_play['alpha'] = 0.01
params_play['gamma'] = 0.03

# run dynamics until when?
time_out = 1000

# MODEL 1
# integrate ODEs 
output_model1 = ode(y = states_model1, times = c(0:time_out), func = ODEs_model1, parms = params_play); 
# print initial outputs and final outputs
head(output_model1); tail(output_model1)
# model compartments should sum to 1
sum(output_model1[time_out,2:3])

# MODEL 2
output_model2 = ode(y = states_model2, times = c(0:time_out), func = ODEs_model2, parms = params_play)
head(output_model2); tail(output_model2)
sum(output_model2[time_out,2:4])

# MODEL 3
output_model3 = ode(y = states_model3, times = c(0:time_out), func = ODEs_model3, parms = params_play)
head(output_model3); tail(output_model3)
sum(output_model3[time_out,2:5])

# MODEL 4
output_model4 = ode(y = states_model4, times = c(0:time_out), func = ODEs_model4, parms = params_play)
head(output_model4); tail(output_model4)
sum(output_model4[time_out,2:7])

# MODEL 5
output_model5 = ode(y = states_model5, times = c(0:time_out), func = ODEs_model5, parms = params_play)
head(output_model5); tail(output_model5)
sum(output_model5[time_out,2:13])

### EQUILIBRIA: univar for each model
univar_par = 'a'
univar_par_range = seq(0,1,0.1)

out_model1 = f_eqbm_univar("model1", ODEs_model1, states_model1, params_play, univar_par, univar_par_range)
out_model2 = f_eqbm_univar("model2", ODEs_model2, states_model2, params_play, univar_par, univar_par_range)
out_model3 = f_eqbm_univar("model3", ODEs_model3, states_model3, params_play, univar_par, univar_par_range)
out_model4 = f_eqbm_univar("model4", ODEs_model4, states_model4, params_play, univar_par, univar_par_range)
out_model5 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_play, univar_par, univar_par_range)


### EQUILIBRIA: bivar for each model
bivar_par1 = 'theta_C'
bivar_par1_range = seq(0,1,0.2)

bivar_par2 = 'theta_m'
bivar_par2_range = seq(0,1,0.2)

out_model1_thetas = f_eqbm_bivar("model1", ODEs_model1, states_model1, params_play, bivar_par1, bivar_par1_range, bivar_par2, bivar_par2_range)
out_model2_thetas = f_eqbm_bivar("model2", ODEs_model2, states_model2, params_play, bivar_par1, bivar_par1_range, bivar_par2, bivar_par2_range)
out_model3_thetas = f_eqbm_bivar("model3", ODEs_model3, states_model3, params_play, bivar_par1, bivar_par1_range, bivar_par2, bivar_par2_range)
out_model4_thetas = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_play, bivar_par1, bivar_par1_range, bivar_par2, bivar_par2_range)
out_model5_thetas = f_eqbm_bivar("model5", ODEs_model5, states_model5, params_play, bivar_par1, bivar_par1_range, bivar_par2, bivar_par2_range)


################
### FIGURE 1 ###
################

range_univar = seq(0,1,0.02)

### find numerical equilibrium (solve ODEs) as a function of "a" (antibiotic exposure prevalence )
### model 1 :
# (a) r_R = 0.8
fig1_model1a = f_eqbm_univar("model1", ODEs_model1, states_model1, params_default, 'a', range_univar)%>%
  dplyr::select(par, par_val, prevalence_C_S, prevalence_C_R, R_rate, incidence_C_S_daily, incidence_C_R_daily)
# (b) r_R = 1
fig1_model1b = f_eqbm_univar("model1", ODEs_model1, states_model1, params_perfectR, 'a', range_univar)%>%
  dplyr::select(par, par_val, prevalence_C_S, prevalence_C_R, R_rate, incidence_C_S_daily, incidence_C_R_daily)


### model 2:
# (a) r_R = 0.8
fig1_model2a = f_eqbm_univar("model2", ODEs_model2, states_model2, params_default, 'a', range_univar)%>%
  dplyr::select(par, par_val, prevalence_C_S, prevalence_C_R, R_rate, incidence_C_S_daily, incidence_C_R_daily)

fig1_model2a_strains = fig1_model2a%>%
  dplyr::select(par, par_val, prevalence_C_S, prevalence_C_R)%>%
  pivot_longer(-c(par, par_val))%>%
  mutate(strain = ifelse(name == "prevalence_C_S", 'Sensitive', 'Resistant'))%>%
  dplyr::select(-name)

# (b) r_R = 1
fig1_model2b = f_eqbm_univar("model2", ODEs_model2, states_model2, params_perfectR, 'a', range_univar)%>%
  dplyr::select(par, par_val, prevalence_C_S, prevalence_C_R, R_rate, incidence_C_S_daily, incidence_C_R_daily)

fig1_model2b_strains = fig1_model2b%>%
  dplyr::select(par, par_val, prevalence_C_S, prevalence_C_R)%>%
  pivot_longer(-c(par, par_val))%>%
  mutate(strain = ifelse(name == "prevalence_C_S", 'Sensitive', 'Resistant'))%>%
  dplyr::select(-name)

### model 3:
## (a) r_R = 0.8
# (i) epsilon
fig1_model3a_epsilon_low = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_epsilon_low, 'a', range_univar)%>%
  mutate(interaction = 'epsilon', value = 'low')
fig1_model3a_epsilon_med = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_epsilon_med, 'a', range_univar)%>%
  mutate(interaction = 'epsilon', value = 'medium')
fig1_model3a_epsilon_high = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_epsilon_high, 'a', range_univar)%>%
  mutate(interaction = 'epsilon', value = 'high')
# (ii) eta
fig1_model3a_eta_low = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_eta_low, 'a', range_univar)%>%
  mutate(interaction = 'eta', value = 'low')
fig1_model3a_eta_med = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_eta_med, 'a', range_univar)%>%
  mutate(interaction = 'eta', value = 'medium')
fig1_model3a_eta_high = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_eta_high, 'a', range_univar)%>%
  mutate(interaction = 'eta', value = 'high')
# (iii) phi
fig1_model3a_phi_low = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_phi_low, 'a', range_univar)%>%
  mutate(interaction = 'phi', value = 'low')
fig1_model3a_phi_med = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_phi_med, 'a', range_univar)%>%
  mutate(interaction = 'phi', value = 'medium')
fig1_model3a_phi_high = f_eqbm_univar("model3", ODEs_model3, states_model3, params_default_phi_high, 'a', range_univar)%>%
  mutate(interaction = 'phi', value = 'high')
# (iv) combined
fig1_model3a_interactions = rbind(fig1_model3a_epsilon_low,
                                  fig1_model3a_epsilon_med,
                                  fig1_model3a_epsilon_high,
                                  fig1_model3a_eta_low,
                                  fig1_model3a_eta_med,
                                  fig1_model3a_eta_high,
                                  fig1_model3a_phi_low,
                                  fig1_model3a_phi_med,
                                  fig1_model3a_phi_high)%>%
  dplyr::select(par, par_val, prevalence_C_R, interaction, value)%>%
  pivot_wider(names_from = value, values_from = prevalence_C_R)

## (b) r_R = 1
# (i) epsilon
fig1_model3b_epsilon_low = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_epsilon_low, 'a', range_univar)%>%
  mutate(interaction = 'epsilon', value = 'low')
fig1_model3b_epsilon_med = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_epsilon_med, 'a', range_univar)%>%
  mutate(interaction = 'epsilon', value = 'medium')
fig1_model3b_epsilon_high = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_epsilon_high, 'a', range_univar)%>%
  mutate(interaction = 'epsilon', value = 'high')
# (ii) eta
fig1_model3b_eta_low = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_eta_low, 'a', range_univar)%>%
  mutate(interaction = 'eta', value = 'low')
fig1_model3b_eta_med = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_eta_med, 'a', range_univar)%>%
  mutate(interaction = 'eta', value = 'medium')
fig1_model3b_eta_high = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_eta_high, 'a', range_univar)%>%
  mutate(interaction = 'eta', value = 'high')
# (iii) phi
fig1_model3b_phi_low = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_phi_low, 'a', range_univar)%>%
  mutate(interaction = 'phi', value = 'low')
fig1_model3b_phi_med = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_phi_med, 'a', range_univar)%>%
  mutate(interaction = 'phi', value = 'medium')
fig1_model3b_phi_high = f_eqbm_univar("model3", ODEs_model3, states_model3, params_perfectR_phi_high, 'a', range_univar)%>%
  mutate(interaction = 'phi', value = 'high')
# (iv) combined
fig1_model3b_interactions = rbind(fig1_model3b_epsilon_low,
                                  fig1_model3b_epsilon_med,
                                  fig1_model3b_epsilon_high,
                                  fig1_model3b_eta_low,
                                  fig1_model3b_eta_med,
                                  fig1_model3b_eta_high,
                                  fig1_model3b_phi_low,
                                  fig1_model3b_phi_med,
                                  fig1_model3b_phi_high)%>%
  dplyr::select(par, par_val, prevalence_C_R, interaction, value)%>%
  pivot_wider(names_from = value, values_from = prevalence_C_R)


################
### FIGURE 2 ###
################

range_bivar = seq(0,0.5,0.025)

fig2_default = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_default, 'theta_C', range_bivar, 'theta_m', range_bivar)

fig2_lowInt_lowR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_lowInt_lowR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_lowInt_medR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_lowInt_medR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_lowInt_highR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_lowInt_highR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_medInt_lowR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_medInt_lowR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_medInt_medR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_medInt_medR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_medInt_highR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_medInt_highR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_highInt_lowR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_highInt_lowR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_highInt_medR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_highInt_medR, 'theta_C', range_bivar, 'theta_m', range_bivar)
fig2_highInt_highR = f_eqbm_bivar("model4", ODEs_model4, states_model4, params_highInt_highR, 'theta_C', range_bivar, 'theta_m', range_bivar)

fig2_mixed = rbind(fig2_lowInt_lowR%>%mutate(label = 1),
                   fig2_lowInt_medR%>%mutate(label = 2),
                   fig2_lowInt_highR%>%mutate(label = 3),
                   fig2_medInt_lowR%>%mutate(label = 4),
                   fig2_medInt_medR%>%mutate(label = 5),
                   fig2_medInt_highR%>%mutate(label = 6),
                   fig2_highInt_lowR%>%mutate(label = 7),
                   fig2_highInt_medR%>%mutate(label = 8),
                   fig2_highInt_highR%>%mutate(label = 9))

fig2_mixed_labels = c("Low interaction strengths\nLow resistance level (r = 0.2)",
                      "Low interaction strengths\nMedium resistance level (r = 0.5)",
                      "Low interaction strengths\nHigh resistance level (r = 0.8)",
                      "Medium interaction strengths\nLow resistance level (r = 0.2)",
                      "Medium interaction strengths\nMedium resistance level (r = 0.5)",
                      "Medium interaction strengths\nHigh resistance level (r = 0.8)",
                      "High interaction strengths\nLow resistance level (r = 0.2)",
                      "High interaction strengths\nMedium resistance level (r = 0.5)",
                      "High interaction strengths\nHigh resistance level (r = 0.8)")

fig2_mixed$label = factor(fig2_mixed$label, levels = 1:9, labels = fig2_mixed_labels)

################
### FIGURE 3 ###
################

# Default parameters
fig3_noInt_noHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_noHGT, 'a', range_univar)
fig3_noInt_lowHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_lowHGT, 'a', range_univar)
fig3_noInt_highHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_highHGT, 'a', range_univar)
fig3_withInt_noHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_noHGT, 'a', range_univar)
fig3_withInt_lowHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_lowHGT, 'a', range_univar)
fig3_withInt_highHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_highHGT, 'a', range_univar)

fig3 = rbind(fig3_noInt_noHGT%>%mutate(HGT = 'null', Interactions = 'none'),
             fig3_noInt_lowHGT%>%mutate(HGT = 'low', Interactions = 'none'),
             fig3_noInt_highHGT%>%mutate(HGT = 'high', Interactions = 'none'),
             fig3_withInt_noHGT%>%mutate(HGT = 'null', Interactions = 'yes'),
             fig3_withInt_lowHGT%>%mutate(HGT = 'low', Interactions = 'yes'),
             fig3_withInt_highHGT%>%mutate(HGT = 'high', Interactions = 'yes'))

# Medium resistance level (r_R = 0.5)
fig3_noInt_noHGT_medR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_noHGT_medR, 'a', range_univar)
fig3_noInt_lowHGT_medR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_lowHGT_medR, 'a', range_univar)
fig3_noInt_highHGT_medR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_highHGT_medR, 'a', range_univar)
fig3_withInt_noHGT_medR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_noHGT_medR, 'a', range_univar)
fig3_withInt_lowHGT_medR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_lowHGT_medR, 'a', range_univar)
fig3_withInt_highHGT_medR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_highHGT_medR, 'a', range_univar)

fig3_medR = rbind(fig3_noInt_noHGT_medR%>%mutate(HGT = 'null', Interactions = 'none'),
                  fig3_noInt_lowHGT_medR%>%mutate(HGT = 'low', Interactions = 'none'),
                  fig3_noInt_highHGT_medR%>%mutate(HGT = 'high', Interactions = 'none'),
                  fig3_withInt_noHGT_medR%>%mutate(HGT = 'null', Interactions = 'yes'),
                  fig3_withInt_lowHGT_medR%>%mutate(HGT = 'low', Interactions = 'yes'),
                  fig3_withInt_highHGT_medR%>%mutate(HGT = 'high', Interactions = 'yes'))

# Low resistance level (r_R = 0.2)
fig3_noInt_noHGT_lowR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_noHGT_lowR, 'a', range_univar)
fig3_noInt_lowHGT_lowR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_lowHGT_lowR, 'a', range_univar)
fig3_noInt_highHGT_lowR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_noInt_highHGT_lowR, 'a', range_univar)
fig3_withInt_noHGT_lowR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_noHGT_lowR, 'a', range_univar)
fig3_withInt_lowHGT_lowR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_lowHGT_lowR, 'a', range_univar)
fig3_withInt_highHGT_lowR = f_eqbm_univar("model5", ODEs_model5, states_model5, params_withInt_highHGT_lowR, 'a', range_univar)

fig3_lowR = rbind(fig3_noInt_noHGT_lowR%>%mutate(HGT = 'null', Interactions = 'none'),
                  fig3_noInt_lowHGT_lowR%>%mutate(HGT = 'low', Interactions = 'none'),
                  fig3_noInt_highHGT_lowR%>%mutate(HGT = 'high', Interactions = 'none'),
                  fig3_withInt_noHGT_lowR%>%mutate(HGT = 'null', Interactions = 'yes'),
                  fig3_withInt_lowHGT_lowR%>%mutate(HGT = 'low', Interactions = 'yes'),
                  fig3_withInt_highHGT_lowR%>%mutate(HGT = 'high', Interactions = 'yes'))


#################
### Figure S5 ###
#################

### dynamics

### no interactions
# initial conditions
params_invasion_none = params_invasion; params_invasion_none['epsilon'] = 0; params_invasion_none['eta'] = 0; params_invasion_none['phi'] = 1; params_invasion_none['chi_e'] = 0; params_invasion_none['chi_d'] = 0; params_invasion_none['omega'] = 0; 
dynamics_none = ode(y = states_invasion, times = c(0:1000), func = ODEs_model5, parms = params_invasion_none, method = 'lsoda')
states_none_preintervention = dynamics_none[nrow(dynamics_none),2:ncol(dynamics_none)]
# from beginning
dynamics_none = ode(y = states_none_preintervention, times = c(0:90), func = ODEs_model5, parms = params_invasion_none, method = 'lsoda')
# introduce intervention
params_invasion_none_intervention = params_invasion_none; params_invasion_none_intervention[par_intervention1] = val_intervention1;
states_none_intervention = dynamics_none[nrow(dynamics_none),2:ncol(dynamics_none)]
dynamics_none_2 = ode(y = states_none_intervention, times = c(91:180), func = ODEs_model5, parms = params_invasion_none_intervention, method = 'lsoda')
# introduce intervention 2
params_invasion_none_intervention2 = params_invasion_none_intervention; params_invasion_none_intervention2[par_intervention2] = val_intervention2;
states_none_intervention2 = dynamics_none_2[nrow(dynamics_none_2),2:ncol(dynamics_none_2)]
dynamics_none_3 = ode(y = states_none_intervention2, times = c(181:270), func = ODEs_model5, parms = params_invasion_none_intervention2, method = 'lsoda')
# introduce intervention 3
params_invasion_none_intervention3 = params_invasion_none_intervention2; params_invasion_none_intervention3[par_intervention3] = val_intervention3;
states_none_intervention3 = dynamics_none_3[nrow(dynamics_none_3),2:ncol(dynamics_none_3)]
dynamics_none_4 = ode(y = states_none_intervention3, times = c(271:364), func = ODEs_model5, parms = params_invasion_none_intervention3, method = 'lsoda')
# combine
dynamics_none_intervention = rbind(dynamics_none, dynamics_none_2, dynamics_none_3, dynamics_none_4)%>%as.data.frame()

### colonization resistance
# initial conditions
params_invasion_epsilon = params_invasion_none; params_invasion_epsilon['epsilon'] = 0.5; 
dynamics_epsilon = ode(y = states_invasion, times = c(0:1000), func = ODEs_model5, parms = params_invasion_epsilon, method = 'lsoda')
states_epsilon_preintervention = dynamics_epsilon[nrow(dynamics_epsilon),2:ncol(dynamics_epsilon)]
# from beginning
dynamics_epsilon = ode(y = states_epsilon_preintervention, times = c(0:90), func = ODEs_model5, parms = params_invasion_epsilon, method = 'lsoda')
# introduce intervention
params_invasion_epsilon_intervention = params_invasion_epsilon; params_invasion_epsilon_intervention[par_intervention1] = val_intervention1
states_epsilon_intervention = dynamics_epsilon[nrow(dynamics_epsilon),2:ncol(dynamics_epsilon)]
dynamics_epsilon_2 = ode(y = states_epsilon_intervention, times = c(91:180), func = ODEs_model5, parms = params_invasion_epsilon_intervention, method = 'lsoda')
# introduce intervention 2
params_invasion_epsilon_intervention2 = params_invasion_epsilon_intervention; params_invasion_epsilon_intervention2[par_intervention2] = val_intervention2;
states_epsilon_intervention2 = dynamics_epsilon_2[nrow(dynamics_epsilon_2),2:ncol(dynamics_epsilon_2)]
dynamics_epsilon_3 = ode(y = states_epsilon_intervention2, times = c(181:270), func = ODEs_model5, parms = params_invasion_epsilon_intervention2, method = 'lsoda')
# introduce intervention 3
params_invasion_epsilon_intervention3 = params_invasion_epsilon_intervention2; params_invasion_epsilon_intervention3[par_intervention3] = val_intervention3;
states_epsilon_intervention3 = dynamics_epsilon_3[nrow(dynamics_epsilon_3),2:ncol(dynamics_epsilon_3)]
dynamics_epsilon_4 = ode(y = states_epsilon_intervention3, times = c(271:364), func = ODEs_model5, parms = params_invasion_epsilon_intervention3, method = 'lsoda')
# combine
dynamics_epsilon_intervention = rbind(dynamics_epsilon, dynamics_epsilon_2, dynamics_epsilon_3, dynamics_epsilon_4)%>%as.data.frame()

### resource competition
# initial conditions
params_invasion_eta = params_invasion_none; params_invasion_eta['eta'] = 0.5; 
dynamics_eta = ode(y = states_invasion, times = c(0:1000), func = ODEs_model5, parms = params_invasion_eta, method = 'lsoda')
states_eta_preintervention = dynamics_eta[nrow(dynamics_eta),2:ncol(dynamics_eta)]
# from beginning
dynamics_eta = ode(y = states_eta_preintervention, times = c(0:90), func = ODEs_model5, parms = params_invasion_eta, method = 'lsoda')
# introduce intervention
params_invasion_eta_intervention = params_invasion_eta; params_invasion_eta_intervention[par_intervention1] = val_intervention1
states_eta_intervention = dynamics_eta[nrow(dynamics_eta),2:ncol(dynamics_eta)]
dynamics_eta_2 = ode(y = states_eta_intervention, times = c(91:180), func = ODEs_model5, parms = params_invasion_eta_intervention, method = 'lsoda')
# introduce intervention 2
params_invasion_eta_intervention2 = params_invasion_eta_intervention; params_invasion_eta_intervention2[par_intervention2] = val_intervention2;
states_eta_intervention2 = dynamics_eta_2[nrow(dynamics_eta_2),2:ncol(dynamics_eta_2)]
dynamics_eta_3 = ode(y = states_eta_intervention2, times = c(181:270), func = ODEs_model5, parms = params_invasion_eta_intervention2, method = 'lsoda')
# introduce intervention 3
params_invasion_eta_intervention3 = params_invasion_eta_intervention2; params_invasion_eta_intervention3[par_intervention3] = val_intervention3;
states_eta_intervention3 = dynamics_eta_3[nrow(dynamics_eta_3),2:ncol(dynamics_eta_3)]
dynamics_eta_4 = ode(y = states_eta_intervention3, times = c(271:364), func = ODEs_model5, parms = params_invasion_eta_intervention3, method = 'lsoda')
# combine
dynamics_eta_intervention = rbind(dynamics_eta, dynamics_eta_2, dynamics_eta_3, dynamics_eta_4)%>%as.data.frame()

### ecological release
# initial conditions
params_invasion_phi = params_invasion_none; params_invasion_phi['phi'] = 5; 
dynamics_phi = ode(y = states_invasion, times = c(0:1000), func = ODEs_model5, parms = params_invasion_phi, method = 'lsoda')
states_phi_preintervention = dynamics_phi[nrow(dynamics_phi),2:ncol(dynamics_phi)]
# from beginning
dynamics_phi = ode(y = states_phi_preintervention, times = c(0:90), func = ODEs_model5, parms = params_invasion_phi, method = 'lsoda')
# introduce intervention
params_invasion_phi_intervention = params_invasion_phi; params_invasion_phi_intervention[par_intervention1] = val_intervention1
states_phi_intervention = dynamics_phi[nrow(dynamics_phi),2:ncol(dynamics_phi)]
dynamics_phi_2 = ode(y = states_phi_intervention, times = c(91:180), func = ODEs_model5, parms = params_invasion_phi_intervention, method = 'lsoda')
# introduce intervention 2
params_invasion_phi_intervention2 = params_invasion_phi_intervention; params_invasion_phi_intervention2[par_intervention2] = val_intervention2;
states_phi_intervention2 = dynamics_phi_2[nrow(dynamics_phi_2),2:ncol(dynamics_phi_2)]
dynamics_phi_3 = ode(y = states_phi_intervention2, times = c(181:270), func = ODEs_model5, parms = params_invasion_phi_intervention2, method = 'lsoda')
# introduce intervention 3
params_invasion_phi_intervention3 = params_invasion_phi_intervention2; params_invasion_phi_intervention3[par_intervention3] = val_intervention3;
states_phi_intervention3 = dynamics_phi_3[nrow(dynamics_phi_3),2:ncol(dynamics_phi_3)]
dynamics_phi_4 = ode(y = states_phi_intervention3, times = c(271:364), func = ODEs_model5, parms = params_invasion_phi_intervention3, method = 'lsoda')
# combine
dynamics_phi_intervention = rbind(dynamics_phi, dynamics_phi_2, dynamics_phi_3, dynamics_phi_4)%>%as.data.frame()



### all
# initial conditions
params_invasion_all = params_invasion_none; params_invasion_all['epsilon'] = 0.5; params_invasion_all['eta'] = 0.5; params_invasion_all['phi'] = 5; #params_invasion_all['chi_e'] = 0.01; params_invasion_all['chi_d'] = 0.05;
dynamics_all = ode(y = states_invasion, times = c(0:1000), func = ODEs_model5, parms = params_invasion_all, method = 'lsoda')
states_all_preintervention = dynamics_all[nrow(dynamics_all),2:ncol(dynamics_all)]
# from beginning
dynamics_all = ode(y = states_all_preintervention, times = c(0:90), func = ODEs_model5, parms = params_invasion_all, method = 'lsoda')
# introduce intervention
params_invasion_all_intervention = params_invasion_all; params_invasion_all_intervention[par_intervention1] = val_intervention1
states_all_intervention = dynamics_all[nrow(dynamics_all),2:ncol(dynamics_all)]
dynamics_all_2 = ode(y = states_all_intervention, times = c(91:180), func = ODEs_model5, parms = params_invasion_all_intervention, method = 'lsoda')
# introduce intervention 2
params_invasion_all_intervention2 = params_invasion_all_intervention; params_invasion_all_intervention2[par_intervention2] = val_intervention2;
states_all_intervention2 = dynamics_all_2[nrow(dynamics_all_2),2:ncol(dynamics_all_2)]
dynamics_all_3 = ode(y = states_all_intervention2, times = c(181:270), func = ODEs_model5, parms = params_invasion_all_intervention2, method = 'lsoda')
# introduce intervention 3
params_invasion_all_intervention3 = params_invasion_all_intervention2; params_invasion_all_intervention3[par_intervention3] = val_intervention3;
states_all_intervention3 = dynamics_all_3[nrow(dynamics_all_3),2:ncol(dynamics_all_3)]
dynamics_all_4 = ode(y = states_all_intervention3, times = c(271:364), func = ODEs_model5, parms = params_invasion_all_intervention3, method = 'lsoda')
# combine
dynamics_all_intervention = rbind(dynamics_all, dynamics_all_2, dynamics_all_3, dynamics_all_4)%>%as.data.frame()

# combine all
dynamics_combined = rbind(dynamics_none_intervention%>%mutate(Interaction = 'none'),
                          dynamics_epsilon_intervention%>%mutate(Interaction = 'colonization resistance'),
                          dynamics_eta_intervention%>%mutate(Interaction = 'resource competition'),
                          dynamics_phi_intervention%>%mutate(Interaction = 'ecological release'),
                          #dynamics_chi_intervention%>%mutate(Interaction = 'HGT'),
                          dynamics_all_intervention%>%mutate(Interaction = 'all'))%>%
  mutate(prevalence_C_R = C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r,
         R_rate = (C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r)/(C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r + C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r))

dynamics_combined$Interaction = factor(dynamics_combined$Interaction,
                                       levels = c('none', 'colonization resistance', 'resource competition', 'ecological release', 'all'))


#################
### Figure S7 ###
#################


# panel A: effect of HGT varying over parameter ranges
pars_varied = c('a', 'r_R', 'c', 'epsilon', 'eta', 'phi', 'beta', 'gamma', 'alpha', 'delta', 'omega', 'f_w')
pars_varied_labels = c('antibiotic exposure prevalence', 'resistance level', 'cost of resistance',
                       'colonization resistance', 'resource competition', 'ecological release',
                       'transmission rate', 'clearance rate', 'endogenous acquisition rate',
                       'dysbiosis recovery rate', 'plasmid acquisition rate', 'admission fraction (plasmid)')

figS7_a_final = data.frame()

for(par_i in pars_varied){
  
  print(par_i)
  
  if(par_i %in% c('a', 'r_R', 'epsilon', 'eta')){iter_min = 0; iter_max = 1; iter_interval = 0.02}
  if(par_i %in% c('c')){iter_min = 0; iter_max = 10; iter_interval = 0.2}
  if(par_i %in% c('phi')){iter_min = 1; iter_max = 10; iter_interval = 0.2}
  if(par_i %in% c('delta')){iter_min = 0; iter_max = 0.5; iter_interval = 0.002*5}
  if(par_i %in% c('beta')){iter_min = 0; iter_max = 0.2; iter_interval = 0.005}
  if(par_i %in% c('alpha', 'gamma', 'f_w', 'omega')){iter_min = 0; iter_max = 0.1; iter_interval = 0.002}
  
  range_univar_loop = seq(iter_min, iter_max, iter_interval)
  
  figS7_a_noHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_noHGT, par_i, range_univar_loop)
  figS7_a_lowHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_lowHGT, par_i, range_univar_loop)
  figS7_a_highHGT = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_highHGT, par_i, range_univar_loop)
  
  figS7_a = rbind(figS7_a_noHGT%>%mutate(rate = 'none')%>%
                    dplyr::select(par, par_val, prevalence_C_R, rate),
                  figS7_a_lowHGT%>%mutate(rate = 'low')%>%
                    dplyr::select(par, par_val, prevalence_C_R, rate), 
                  figS7_a_highHGT%>%mutate(rate = 'high')%>%
                    dplyr::select(par, par_val, prevalence_C_R, rate))%>%
    pivot_wider(values_from = prevalence_C_R, names_from = rate)%>%
    mutate(diff_none = none - none,
           diff_low = low - none,
           diff_high = high - none)
  
  figS7_a_final = rbind(figS7_a_final, figS7_a)
  
}

figS7_a_final$par = factor(figS7_a_final$par, levels = pars_varied, labels = pars_varied_labels)

# panel B: varying chi_d/chi_e
range_univar = seq(0,1,0.02)
figS7_b_varyHGT1 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_varyHGT1, 'a', range_univar)
figS7_b_varyHGT2 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_varyHGT2, 'a', range_univar)
figS7_b_varyHGT3 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_varyHGT3, 'a', range_univar)
figS7_b_varyHGT4 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_varyHGT4, 'a', range_univar)
figS7_b_varyHGT5 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_varyHGT5, 'a', range_univar)

figS7_varyHGT = rbind(figS7_b_varyHGT1%>%mutate(ratio = 1),
                      figS7_b_varyHGT2%>%mutate(ratio = 2),
                      figS7_b_varyHGT3%>%mutate(ratio = 4),
                      figS7_b_varyHGT4%>%mutate(ratio = 8),
                      figS7_b_varyHGT5%>%mutate(ratio = 16))

# panel C: varying c over a

# cost 1: c=-0.5 (actually a fitness benefit)
figS7_c_vary_c0_1 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c0_1, 'a', range_univar)%>%
  mutate(prevalence0 = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence0)
figS7_c_vary_c1_1 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c1_1, 'a', range_univar)%>%
  mutate(prevalence = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence)
figS7_c_vary_1 = left_join(figS7_c_vary_c0_1, figS7_c_vary_c1_1)%>%mutate(prev_difference = prevalence - prevalence0)%>%
  mutate(cost = -0.5)

# cost 2: c=0 (no cost)
figS7_c_vary_c0_2 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c0_2, 'a', range_univar)%>%
  mutate(prevalence0 = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence0)
figS7_c_vary_c1_2 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c1_2, 'a', range_univar)%>%
  mutate(prevalence = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence)
figS7_c_vary_2 = left_join(figS7_c_vary_c0_2, figS7_c_vary_c1_2)%>%mutate(prev_difference = prevalence - prevalence0)%>%
  mutate(cost = 0)

# cost 3: c=1 (baseline cost)
figS7_c_vary_c0_3 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c0_3, 'a', range_univar)%>%
  mutate(prevalence0 = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence0)
figS7_c_vary_c1_3 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c1_3, 'a', range_univar)%>%
  mutate(prevalence = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence)
figS7_c_vary_3 = left_join(figS7_c_vary_c0_3, figS7_c_vary_c1_3)%>%mutate(prev_difference = prevalence - prevalence0)%>%
  mutate(cost = 1)

# cost 4: c=2 (higher cost)
figS7_c_vary_c0_4 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c0_4, 'a', range_univar)%>%
  mutate(prevalence0 = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence0)
figS7_c_vary_c1_4 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c1_4, 'a', range_univar)%>%
  mutate(prevalence = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence)
figS7_c_vary_4 = left_join(figS7_c_vary_c0_4, figS7_c_vary_c1_4)%>%mutate(prev_difference = prevalence - prevalence0)%>%
  mutate(cost = 2)

# cost 5: c=4 (highest cost)
figS7_c_vary_c0_5 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c0_5, 'a', range_univar)%>%
  mutate(prevalence0 = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence0)
figS7_c_vary_c1_5 = f_eqbm_univar("model5", ODEs_model5, states_model5, params_HGTsupp_vary_c1_5, 'a', range_univar)%>%
  mutate(prevalence = prevalence_C_R + prevalence_C_S)%>%
  dplyr::select(par, par_val, prevalence)
figS7_c_vary_5 = left_join(figS7_c_vary_c0_5, figS7_c_vary_c1_5)%>%mutate(prev_difference = prevalence - prevalence0)%>%
  mutate(cost = 4)
  
# combine 
figS7_c_vary = rbind(figS7_c_vary_1, figS7_c_vary_2, figS7_c_vary_3, figS7_c_vary_4, figS7_c_vary_5)




