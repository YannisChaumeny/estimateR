delayDistr <-
  function(mean_incub,
           sd_incub,
           mean_test,
           sd_test,
           n_iter,
           days_max) {
    shape_incub <- mean_incub ^ 2 / sd_incub ^ 2
    rate_incub  <- mean_incub / sd_incub ^ 2
    shape_test <- mean_test ^ 2 / sd_test ^ 2
    rate_test  <- mean_test / sd_test ^ 2
    sample_incub <- rgamma(n_iter, shape_incub, rate_incub)
    sample_test <- rgamma(n_iter, shape_test, rate_test)
    return(diff(ecdf(sample_incub + sample_test)(1:(days_max + 1)-0.5)))
  }


fullDelay <- function(N, dates, gammaDelay) {
  d <- matrix(0, N, N)
  weekdays <- format(dates, format = "%a")
  for (t in 1:(N - 1)) {
    for (s in t:N) {
      if (weekdays[s] == "Sat" | weekdays[s] == "Sun") {
        d[t, s] = 0.6 * gammaDelay[s - t + 1]
      }
      else{
        d[t, s] = gammaDelay[s - t + 1]
      }
    }
  }
  
  for (t in 1:(N - 1)) {
    d[t, 1:N] = d[t, 1:N] / sum(d[t, 1:N])
  }
  return(d)
}

belgianDelay <- function(N, dates, gammaDelay) {
  d <- matrix(0, N, N)
  weekdays <- format(dates, format = "%a")
  for (t in 1:(N - 1)) {
    for (s in t:N) {
      if (weekdays[s] == "Sun" | weekdays[s] == "Mon") {
        d[t, s] = 0.6 * gammaDelay[s - t + 1]
      }
      else{
        d[t, s] = gammaDelay[s - t + 1]
      }
    }
  }
  
  for (t in 1:(N - 1)) {
    d[t, 1:N] = d[t, 1:N] / sum(d[t, 1:N])
  }
  return(d)
}

simulated_delay_dist <- function(nn) {
  r_inc_dist <-
    function(n) {
      rgamma(n, shape = (5.3 / 3.2) ^ 2, rate = 5.3 / (3.2 ^ 2))
    } # Incubation period (infection -> symptoms)
  r_sym_to_obs_dist <-
    function(n) {
      rgamma(n, shape = (4.5 / 4.9) ^ 2, rate = 4.5 / (4.9 ^ 2))
    } # Additional delay from symptoms -> observation
  r_inc_dist(nn) + r_sym_to_obs_dist(nn)
}

get_tObs_from_tInf <- function(n_dS,
                               times,
                               r_delay_dist,
                               return_times = FALSE) {
  stopifnot(length(n_dS) == length(times))
  n_dS <- ifelse(is.na(n_dS), 0, n_dS)
  stopifnot(n_dS == round(n_dS))
  
  
  # This function draws times of observation for each of the ndS cases incident at a given time point
  get_obs_times_for_one_timestep <- function(ndS, tt) {
    (tt -
       runif(ndS) +  ## Subtract a uniform between 0 and 1 to get exact time of day of onset on the previous day
       r_delay_dist(ndS))  %>%
      ceiling()  # Round up: infections are recorded at the end of the day in progress
  }
  
  ## Get a vector of imputed observation times for each infection
  obs_time_vec <-
    mapply(FUN = get_obs_times_for_one_timestep, ndS = n_dS, tt = times) %>%
    unlist()
  stopifnot(length(obs_time_vec) == sum(n_dS))
  
  ## Reformat: count the number of observed infections at each time, and output
  data.frame(time = obs_time_vec) %>%
    group_by(time) %>%
    summarise(n = n()) %>%
    # Pad with 0s at times with no observed cases
    complete(time = times, fill = list(n = 0)) -> out
  
  if (return_times)
    out
  else
    out$n
}



# time from infection to symptom onset ~ lognormal
# references Lauer et al., Liu et al., Linton et al., Backer, Pellis (see assumptions nowcasting Christel)
delayIncubation <- function(nsim = 100000,
                            mu = 1.63,
                            sd = 0.28) {
  u <- rlnorm(nsim, meanlog = mu, sdlog = sd)
  delay.time <- u
  delay.prob <- table(round(delay.time)) / length(delay.time)
  m <- min(as.numeric(names(delay.prob)))
  if (m == 0) {
    pds <- delay.prob[-1]
  }
  if (m > 0) {
    pds <- c(rep(0, (m - 1)), delay.prob)
  }
  return(pds)
}


# delay time between symptom onset and confirmed test ~ weibull
# based on Belgian data (see updated results nowcasting Christel)
BELdelayCase <-
  function(nsim = 100000,
           mu = 1.63,
           sd = 0.28,
           shape = 1.15164,
           scale = 6.362063) {
    u <- rlnorm(nsim, meanlog = mu, sdlog = sd)
    r <- rweibull(nsim, shape = shape, scale = scale)
    delay.time <- u + r
    delay.prob <- table(round(delay.time)) / length(delay.time)
    m <- min(as.numeric(names(delay.prob)))
    if (m == 0) {
      pd <- delay.prob[-1]
    }
    if (m > 0) {
      pd <- c(rep(0, (m - 1)), delay.prob)
    }
    return(pd)
  }

# delay time between symptom onset and hospitalization ~ weibull
# based on Belgian data (see updated results nowcasting Christel)
BELdelayHosp <-
  function(nsim = 100000,
           mu = 1.63,
           sd = 0.28,
           shape=1.112,
           scale=5.970) {
    u <- rlnorm(nsim, meanlog = mu, sdlog = sd)
    r <- rweibull(nsim, shape = shape, scale = scale)
    delay.time <- u + r
    delay.prob <- table(round(delay.time)) / length(delay.time)
    m <- min(as.numeric(names(delay.prob)))
    if (m == 0) {
      pd <- delay.prob[-1]
    }
    if (m > 0) {
      pd <- c(rep(0, (m - 1)), delay.prob)
    }
    return(pd)
  }
