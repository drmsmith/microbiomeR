### over single set of ODE output:
f_test_dynamics = function(out){
  out_df = data.frame(out)%>%dplyr::select(-time, -dC_S_trans, -dC_S_acq, -dC_R_trans, -dC_R_acq, -dC_R_hgt)
  
  # test constant N
  ones = out_df%>%(rowSums)
  if(all.equal(length(ones),sum(ones))){}else{stop("Non-constant N")}
  
  # test equilibrated by final t
  for(i in 1:ncol(out_df)){
    if(nth(out_df[,i],-1)-nth(out_df[,i],-3)<0.0001){
      if(nth(out_df[,i],-1)-nth(out_df[,i],-2)<0.0001){}else{stop("Periodicity?")}}
    else{stop("Not at equilibrium")}
  }
}


### f_eqbm_univar
# Function to produce equilibrium values while varying one parameter over a range
f_eqbm_univar = function(model, ODEs, states_input, par_input,var, range){
  if(!model %in% c('model1', 'model2', 'model3', 'model4', 'model5')){warning('select a model'); stop()}
  max_time = 10000
  time_step = 1000
  output = data.frame()
  for(i in range){
    par_input[var] = i
    out <- ode(y = states_input, times = seq(0, max_time, by = time_step), func = ODEs, parms = par_input, method = 'lsoda')
    
    ### determine system has found equilibrium (error function; will stop if not)
    f_test_dynamics(out)
    
    ### model 1
    if(model == 'model1'){
      out_df = data.frame(out)%>%
        filter(time %in% c(max_time, max_time - time_step))%>%
        cbind(par = var)%>%
        cbind(par_val = i)%>%
        mutate(prevalence_C_S = 0,
               prevalence_C_R = C_R,
               R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
               incidence_C_S = dC_S_trans + dC_S_acq,
               incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
    }
    ### model 2
    if(model == 'model2'){
      out_df = data.frame(out)%>%
        filter(time %in% c(max_time, max_time - time_step))%>%
        cbind(par = var)%>%
        cbind(par_val = i)%>%
        mutate(prevalence_C_S = C_S,
               prevalence_C_R = C_R,
               R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
               incidence_C_S = dC_S_trans + dC_S_acq,
               incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
    }
    
    ### model 3
    if(model == 'model3'){
      out_df = data.frame(out)%>%
        filter(time %in% c(max_time, max_time - time_step))%>%
        cbind(par = var)%>%
        cbind(par_val = i)%>%
        mutate(prevalence_C_S = 0,
               prevalence_C_R = C_R_e + C_R_d,
               R_rate = 1,
               incidence_C_S = dC_S_trans + dC_S_acq,
               incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
    }
    
    ### model 4
    if(model == 'model4'){
      out_df = data.frame(out)%>%
        filter(time %in% c(max_time, max_time - time_step))%>%
        cbind(par = var)%>%
        cbind(par_val = i)%>%
        mutate(prevalence_C_S = C_S_e + C_S_d,
               prevalence_C_R = C_R_e + C_R_d,
               R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
               incidence_C_S = dC_S_trans + dC_S_acq,
               incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
    }
    
    ### model 5
    if(model == 'model5'){
      out_df = data.frame(out)%>%
        filter(time %in% c(max_time, max_time - time_step))%>%
        cbind(par = var)%>%
        cbind(par_val = i)%>%
        mutate(prevalence_C_S = C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r,
               prevalence_C_R = C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r,
               R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
               incidence_C_S = dC_S_trans + dC_S_acq,
               incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
    }
    
    # update data: calculate daily incidence from cumulative incidence
    out_df_dailyIncidence = out_df%>%
      dplyr::select(time, incidence_C_S, incidence_C_R)%>%
      summarise(time = diff(time), incidence_C_S = diff(incidence_C_S), incidence_C_R = diff(incidence_C_R))%>%
      mutate(incidence_C_S_daily = incidence_C_S/time,
             incidence_C_R_daily = incidence_C_R/time)%>%
      dplyr::select(incidence_C_S_daily, incidence_C_R_daily)%>%
      cbind(out_df%>%filter(time == max_time),.)
    
    
    
    # ones = out_df%>%dplyr::select(-time, -par, -par_val, 
    #                               -dC_S_trans, -dC_S_acq, -dC_R_trans, -dC_R_acq, -dC_R_hgt, 
    #                               -incidence_C_S, -incidence_C_R)%>%(rowSums)
    # if(ones<0.999 | ones > 1.001){warning('N != 1'); stop("N != 1")}
    
    output = rbind(output,out_df_dailyIncidence)
  }
  return(output)
}

### f_eqbm_bivar
# Function to produce equilibrium values while varying two parameters, each over a range
f_eqbm_bivar = function(model, ODEs, states_input, par_input,var1, range1, var2, range2){
  if(!model %in% c('model1', 'model2', 'model3', 'model4', 'model5')){warning('select a model'); stop()}
  max_time = 10000
  time_step = 1000
  output = data.frame()
  for(i in range1){
    par_input[var1] = i
    for(j in range2){
      par_input[var2] = j
      out <- ode(y = states_input, times = seq(0, max_time, by = time_step), func = ODEs, parms = par_input, method = 'lsoda')
      
      ### determine system has found equilibrium (error function; will stop if not)
      f_test_dynamics(out)
      
      ### model 1
      if(model == 'model1'){
        out_df = data.frame(out)%>%
          filter(time %in% c(max_time, max_time - time_step))%>%
          cbind(par1 = var1)%>%
          cbind(par1_val = i)%>%
          cbind(par2 = var2)%>%
          cbind(par2_val = j)%>%
          mutate(prevalence_C_S = 0,
                 prevalence_C_R = C_R,
                 R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
                 incidence_C_S = dC_S_trans + dC_S_acq,
                 incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
      }
      
      ### model 2
      if(model == 'model2'){
        out_df = data.frame(out)%>%
          filter(time %in% c(max_time, max_time - time_step))%>%
          cbind(par1 = var1)%>%
          cbind(par1_val = i)%>%
          cbind(par2 = var2)%>%
          cbind(par2_val = j)%>%
          mutate(prevalence_C_S = C_S,
                 prevalence_C_R = C_R,
                 R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
                 incidence_C_S = dC_S_trans + dC_S_acq,
                 incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
      }
      
      ### model 3
      if(model == 'model3'){
        out_df = data.frame(out)%>%
          filter(time %in% c(max_time, max_time - time_step))%>%
          cbind(par1 = var1)%>%
          cbind(par1_val = i)%>%
          cbind(par2 = var2)%>%
          cbind(par2_val = j)%>%
          mutate(prevalence_C_S = 0,
                 prevalence_C_R = C_R_e + C_R_d,
                 R_rate = 1,
                 incidence_C_S = dC_S_trans + dC_S_acq,
                 incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
          
      }
      
      ### model 4
      if(model == 'model4'){
        out_df = data.frame(out)%>%
          filter(time %in% c(max_time, max_time - time_step))%>%
          cbind(par1 = var1)%>%
          cbind(par1_val = i)%>%
          cbind(par2 = var2)%>%
          cbind(par2_val = j)%>%
          mutate(prevalence_C_S = C_S_e + C_S_d,
                 prevalence_C_R = C_R_e + C_R_d,
                 R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
                 incidence_C_S = dC_S_trans + dC_S_acq,
                 incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
      }
      
      ### model 5
      if(model == 'model5'){
        out_df = data.frame(out)%>%
          filter(time %in% c(max_time, max_time - time_step))%>%
          cbind(par1 = var1)%>%
          cbind(par1_val = i)%>%
          cbind(par2 = var2)%>%
          cbind(par2_val = j)%>%
          mutate(prevalence_C_S = C_S_e_s + C_S_e_r + C_S_d_s + C_S_d_r,
                 prevalence_C_R = C_R_e_s + C_R_e_r + C_R_d_s + C_R_d_r,
                 R_rate = prevalence_C_R/(prevalence_C_S+prevalence_C_R),
                 incidence_C_S = dC_S_trans + dC_S_acq,
                 incidence_C_R = dC_R_trans + dC_R_acq + dC_R_hgt)
      }
      
      # update data: calculate daily incidence from cumulative incidence
      out_df_dailyIncidence = out_df%>%
        dplyr::select(time, incidence_C_S, incidence_C_R)%>%
        summarise(time = diff(time), incidence_C_S = diff(incidence_C_S), incidence_C_R = diff(incidence_C_R))%>%
        mutate(incidence_C_S_daily = incidence_C_S/time,
               incidence_C_R_daily = incidence_C_R/time)%>%
        dplyr::select(incidence_C_S_daily, incidence_C_R_daily)%>%
        cbind(out_df%>%filter(time == max_time),.)
      
      output = rbind(output,out_df_dailyIncidence)

    }
  }
  return(output)
}
