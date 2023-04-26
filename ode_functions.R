# Model function
COVID.base = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Calculate the total population size
    Population = S + E + IA + IS + C
    
    # Calculate intervention efficacy
    beta = (1-lambda)*beta
    
    # Calculate the average force of infection imposed on each susceptible individual
    force_of_infection = beta*(psi*IA + IS)/Population
    
    # Calculate the net (instantaneous) change in each state variable
    dSdt = A - force_of_infection*S - mu*S
    dEdt = force_of_infection*S - E*(alphaA + alphaS + mu) 
    dIAdt = E*alphaA - IA*(omega + deltaA + theta + mu)
    dISdt = E*alphaS + IA*omega - IS*(deltaS + epsilonI + mu)
    dCdt = IA*deltaA + IS*(deltaS)
    
    # Return net changes as list
    return(list(
      c(
        dSdt,
        dEdt,
        dIAdt,
        dISdt,
        dCdt
      )
    ))
  })
}

solve.base.model = function(y_ = state,
                            times_ = times,
                            func. = COVID.base,
                            parms_ = parameters) {
  
  out = ode(
    y = y_,
    times = as.numeric(times_ - times_[1]),
    func = func.,
    parms = parms_
  )
  
  # Calculate the prevalence, incidence and cumulative incidence (for comparison with data)
  out = as.data.frame(out) %>%
    mutate(
      # Prevalence = E + IA + IS + C,
      Incidence = E*parms_["alphaA"] + E*parms_["alphaS"],
      Cumulative_incidence = cumsum(Incidence) + Incidence[1],
      # Cumulative_C = cumsum(C),
      Population = S + E + IA + IS + C,
      Date = times_)
  
  return(out)
}

negative_log_likelihood = function(transformed_parameters,
                                   data = uncontrolled_period$Cases[-1],
                                   state_base = state,
                                   times_ = times,
                                   func. = COVID.base,
                                   parms_base = parameters
) {
  
  # Untransform parameters
  E_0 = exp(transformed_parameters[1])
  IA_0 = exp(transformed_parameters[2])
  IS_0 = exp(transformed_parameters[3])
  lambda = exp(transformed_parameters[4])
  deltaA = exp(transformed_parameters[5])
  deltaS = exp(transformed_parameters[6])
  
  # Overwrite baseline parameters with proposed parameters
  state_base["E"] = E_0
  state_base["IA"] = IA_0
  state_base["IS"] = IS_0
  parms_base["lambda"] = lambda
  parms_base["deltaA"] = deltaA
  parms_base["deltaS"] = deltaS
  
  # Solve model with updated parameters
  out = solve.base.model(state_base,
                         times_,
                         func.,
                         parms_base)
  
  return(-sum(dpois(x=data,
                    lambda=out$Incidence[-1],
                    log=TRUE)))
}
