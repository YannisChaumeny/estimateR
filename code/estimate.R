library(lubridate)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source('code/delay.R')

estimate <- function(cases,
                     date = NULL,
                     SI_mean = 4.8, #Linton et al-
                     SI_std = 2.3,
                     mean_incub = 5.3,
                     std_incub = 3.2,
                     mean_test = 4.5,
                     std_test = 4.9,
                     sigma = 0.05,
                     simulation = FALSE,
                     delay_type = 1,
                     N_pad = 10,
                     warmup = 400,
                     iter = 1000,
                     chains = 4
                     ){
  N <- length(cases)
  if(simulation){
    wkday <- rep(1:7,ceiling(N/7))[1:N] #start fake weekdays at Monday = 1
  }
  else{
    if(!(length(cases) == length(date))) stop("Cases and dates have different length")
    wkday <- wday(as.Date(date),week_start=1)
  }
  
  I_upper <- round(10*max(cases))
  d <- delayDistr(mean_incub,std_incub,mean_test,std_test,1e5,30) #remove hardcoded 30
  C_pad <- rep(round(sum(cases[(N-6):N])/7),N_pad)
  data <- list(
    C = cases,
    N_obs = N,
    wkday = wkday,
    I_upper = I_upper,
    sigmaR = sigma,
    SI_mean = SI_mean,
    SI_std = SI_std,
    d = d,
    K = length(d),
    S = 20,    #move w outside of stan and determine when w goes below a tolerance?
    delay = delay_type,
    N_pad = N_pad,
    C_pad = C_pad
  )
  
  print(format(Sys.time(), "%a %b %d %X %Y"))
  stan_fit <- stan(file = 'estimateR.stan',data = data, warmup = warmup, iter = iter, chains = chains)
  
  return(list(
    R_summary = summary(stan_fit,pars="R")$summary,
    I_summary = summary(stan_fit,pars="I")$summary,
    pr_summary = summary(stan_fit,pars="pr")$summary
  ))
  
}