source("code/simulation.R")
source("code/delay.R")

#-------------------------------------------------------------------------------
#Code from Gostic et al. (2020) https://github.com/cobeylab/Rt_estimation
#-------------------------------------------------------------------------------

simulationData <- function(
  weekend = TRUE, #simulate weekend delay patterns
  p = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1),
  N = 2e6, #total population size
  E_init = 0,
  I_init = 60,
  t_E = 4, # mean time in E (latent period)
  t_I = 4, # mean time in I (duration of infectiousness)
  n_t = 300, # total timesteps
  pre_intervention_R0 = 2.0, # Initial R0 before interventions
  intervention_R0 = 0.8, # Final R0 after interventions
  partially_lifeted_R0 = 1.15,
  intervention_time_1 = 60, # Timepoint at which intervention starts (at which underlying transmission rate begins to fall)
  intervention_time_2 = 60+30,
  days_intervention_to_min = c(7), # Days from intervention start until transmission rate hits min_R0
  days_to_Rt_rise = 7,
  model_types = c('seir'), # Can also choose sir
  methods = c('ode', 'stochastic') # could also choose ode
)
{
  sim_list <- simulation(N,E_init,I_init,t_E,t_I,n_t,pre_intervention_R0,intervention_R0,
                         partially_lifeted_R0,intervention_time_1,intervention_time_2,
                         days_intervention_to_min,days_to_Rt_rise,model_types,methods)
  trueI <- round(sim_list$sim_df$incidence)
  trueI <- trueI[-c(1,length(trueI))]
  N <- length(trueI)
  K <- 25
  d <- delayDistr(5.3,3.2,4.5,4.9,100000,K)
  C <- get_tObs_from_tInf(trueI,1:N,simulated_delay_dist)
  N <- length(trueI)
  C <- C[1:N]
  trueRt <- sim_list$sim_df$true_rt[-c(1,length(sim_list$sim_df$true_rt))]
  GI_mean <- (t_E+t_I)
  GI_var <- 2*(GI_mean/2)^2
  w_shape <- GI_mean^2/GI_var
  w_rate <- GI_mean/GI_var
  w <- dgamma((1:20),shape = w_shape, rate = w_rate)
  S <- length(w)
  
  if(weekend){
    wkday <- sapply(0:(N-1), function(x) x%%7 + 1) # start on 1 = Monday
    dC <- rep(0,N)
    for(i in 1:(N-2)){
      if(C[i]>0){
        if(wkday[i] < 6) #Normal week days
        {
          delayed <- rbinom(1,C[i],p[wkday[i]])
          dC[i+1] <- dC[i+1] + delayed
          dC[i] <- dC[i] + (C[i] - delayed)
        }
        else if(wkday[i] == 6) #Saturday
        {
          delayed <- rbinom(1,C[i],p[6] + p[7])
          delayedMon <- rbinom(1,delayed, p[6]/(p[6] + p[7]))
          delayedSun <- delayed-delayedMon
          dC[i] <- dC[i] + (C[i] - delayed)
          dC[i+1] <- dC[i+1] + delayedSun
          dC[i+2] <- dC[i+2] + delayedMon
        }
        else #Sunday
        {
          delayed <- rbinom(1,C[i],p[8])
          dC[i+1] <- dC[i+1] + delayed
          dC[i] <- dC[i] + (C[i] - delayed)
        }  
      }

    } #for
    C <- dC
  }# if(weekend)
  else{
    wkday <- NULL
    p <- NULL
  }
  return(list(C = C,N = N,d = d,K = K,trueI = trueI, trueRt = trueRt,
              GI_mean = GI_mean, GI_std = sqrt(GI_var), wkday = wkday))
}

