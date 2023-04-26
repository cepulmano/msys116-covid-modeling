##########################################
## Modeling COVID-19 in the Philippines ##
##########################################

# Import the required packages
required_packages <- (c("dplyr","deSolve","ggplot2"))
for (pkg in required_packages) {
  if (!pkg %in% rownames(installed.packages())) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}  

source('./ode_functions.R')

# Import dataset from the Philippines COVID-19 tracker: https://doh.gov.ph/covid19tracker
df = read.csv("covidtracker.csv")
head(df)

# Transform the dataset to get the daily new cases
df = select(df,Report_Date)
df$Report_Date <- as.Date(df$Report_Date,format ="%Y-%m-%d")

first_wave = df %>% group_by(Report_Date) %>% summarise(n = n())
colnames(first_wave) <- c('Date','Cases')
head(first_wave)

# Set inclusive period
start_date = as.Date("2020-04-14")
end_date = start_date + 14 - 1
end_modeling_date = end_date + 14

# Initial conditions
N = 109035343  # Population of Philippines
S_0 = 0.75*N
E_0 = 1 # for fitting
IA_0 = 2 # for fitting
IS_0 = 2 # for fitting
C_0 = 5132

# Initial model parameters taken from
# Literature
# Actual data
parameters = c(
  mu = 4.05e-5,
  A =  N*0.020177/365,
  beta = 0.4343,
  lambda = 0.1, # for fitting
  psi = 1,
  alphaA = 0.2,
  alphaS = 0.2,
  omega = 0.33,
  epsilonI = 0.0031,
  epsilonT = 0.0031,
  deltaA = 0.2, #for fitting
  deltaS = 0.2, #for fitting
  theta = 0.714
)

# State variables
state = c(S = S_0,
          E = E_0,
          IA = IA_0,
          IS = IS_0,
          C = C_0)

## Define model as per the previous session
# Time window
times = seq(start_date, end_date, by=1)

# transform the parameters
# The search / optimization algorithm we employ searches over
# the range (-Inf, Inf) for each variable. In our application
# we are only interested in solutions in the range (0, Inf)
# Therefore, it would be a good idea to apply a transformation to our
# search space so that we are not wasting time exploring infeasible regions of parameter space.
initial_transformed_parameters = log(c(1,2,2,0.1,0.2,0.2))
initial_transformed_parameters

# Filter data to capture period prior to interventions. This has already been done in the last session.
uncontrolled_period = first_wave %>%
  filter(Date >= start_date, Date <= end_date)

# Solve the model using the initial set of parameters
out_init = solve.base.model(y_ = state,
                            times_ = times,
                            func. = COVID.base,
                            parms_ = parameters) 

# Plot the filtered data and the first model "guess"
ggplot(uncontrolled_period) +
  geom_col(aes(x=Date, y=Cases), width=1, fill="dodgerblue2", colour="black") +
  geom_line(data=out_init, aes(x=Date, y=Incidence), size=2, colour="darkgreen") +
  ylab("Daily cases") +
  xlab("") +
  ggtitle("Philippines' Daily Confirmed Cases - Guess parameters") +
  theme_bw()

# Run the fitting procedure to find optimal parameters
optim_out <- optim(initial_transformed_parameters,
                   fn = negative_log_likelihood,
                   lower = c(log(1), log(1), log(1), log(0.01), log(0.01), log(0.01)),
                   upper = c(log(1000), log(1000), log(1000), log(1), log(1), log(1)),
                   method = "L-BFGS-B",
                   control = list(maxit=500),
                   hessian = TRUE
)

# Check for convergence
optim_out$convergence # 0 - converged; 1 - failed to converge

# Inspect solution
optim_out$par

# Back-transform parameters
optimum_parameters <- exp(optim_out$par)
optimum_parameters

# Inspect the optimal negative log-likelihood
optimum_NLL <- optim_out$value

#### Optimal model solution ####

# # Create the optimal parameter and state vectors
optimal_initial_state = state
optimal_initial_state["E"] = optimum_parameters[1]
optimal_initial_state["IA"] = optimum_parameters[2]
optimal_initial_state["IS"] = optimum_parameters[3]

optimal_parameters = parameters
optimal_parameters["lambda"] = optimum_parameters[4]
optimal_parameters["deltaA"] = optimum_parameters[5]
optimal_parameters["deltaS"] = optimum_parameters[6]

# Solve the model given the optimal parameters and initial conditions
optimal_solution = solve.base.model(
  y_ = optimal_initial_state,
  times_ = times,
  func. = COVID.base,
  parms = optimal_parameters
)

# Plot the optimal solution
ggplot(uncontrolled_period) +
  geom_col(aes(x=Date, y=Cases), width=1, fill="dodgerblue2", colour="black") +
  geom_line(data=optimal_solution, aes(x=Date, y=Incidence), size=2, colour="darkgreen") +
  ylab("Daily cases") +
  xlab("") +
  ggtitle("Philippines' Daily Confirmed Cases - Fitted") +
  theme_bw()

# Projection using optimal parameters based on the optimization algorithm

# Solve the model given the optimal parameters and initial conditions
times_extended = seq(start_date, end_modeling_date, by=1)
optimal_solution = solve.base.model(
  y_ = optimal_initial_state,
  times_ = times_extended,
  func. = COVID.base,
  parms = optimal_parameters
)

actual_cases_training = first_wave %>%
  filter(Date >= start_date, Date <= end_date)
actual_cases_projections = first_wave %>%
  filter(Date > end_date, Date <= end_modeling_date)
forecast_base = first_wave %>%
  filter(Date > end_date, Date <= end_modeling_date)
trained = optimal_solution %>%
  filter(Date <= end_date)
forecasted = optimal_solution %>%
  filter(Date > end_date, Date <= end_modeling_date)

poisson_log_likelihood_projected = -sum(dpois(forecast_base$Cases,forecasted$Incidence, log = TRUE))

# negative log likelihood for the projected data
poisson_log_likelihood_projected

# Plot the model optimal model projections
ggplot(actual_cases_training) +
  geom_col(aes(x=Date, y=Cases), width=1, fill="dodgerblue2", colour="black") +
  geom_col(data=actual_cases_projections, aes(x=Date, y=Cases), width=1, fill="firebrick", colour="black") +
  geom_line(data=optimal_solution, aes(x=Date, y=Incidence), size=2, colour="darkgreen") +
  ylab("Daily cases") +
  xlab("") +
  ggtitle("Projections using L-BFGS-B") +
  theme_bw()