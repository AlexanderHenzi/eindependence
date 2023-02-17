#-------------------------------------------------------------------------------
# functions and packages
source("simulation_functions.R")

#-------------------------------------------------------------------------------
# parameters
n <- 2000
p <- 4 # or p <- 8 for higher dimension
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
betap <- seq(0, 1, 0.1)
betap <- betap[id %% length(betap) + 1]
n0 <- 5 * p

#-------------------------------------------------------------------------------
# data generation
set.seed(id)
Mu <- rep(0, p)
Sigma <- toeplitz(1 / seq_len(p) * c(1, -1)) # change to c(1, -1) for negative correlations
z <- MASS::mvrnorm(n = n, mu = Mu, Sigma = Sigma)
x <- z[, p]
z <- z[, -p, drop = FALSE]
intercept <- runif(1, -1, 1)
pars <- runif(p - 1, -1, 1)
probs <- expit(intercept + z %*% pars + x * betap)
y <- rbinom(n, 1, probs)
sim_fun <- function(z) sim_normal_covariate(z, Mu = Mu, Sigma = Sigma, nsim = 500)
n0 <- 5 * p

#-------------------------------------------------------------------------------
# different testing methods

## e-values with conditional randomization and different truncation thresholds
# ecrt0 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0)
# ecrt0.01 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0.01)
ecrt0.05 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0.05)
# ecrt0.1 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0.1)

## oracle conditional randomization
oracle <- e_crt_logistic_oracle(x = x, y = y, z = z, sim_fun = sim_fun, par = c(intercept, pars, betap))

## running mle
rmle <- running_mle_logistic(x = x, y = y, z = z, n0 = n0)

## penalized running mle
rmlep <- running_mle_logistic_penalized(x = x, y = y, z = z, n0 = n0, not_penalize = NULL) # not_penalize = 1 for not penalizing the "x"-coefficient

## p-values
pvals <- crt_lrt_logistic(x = x, y = y, z = z, stops = seq(25, n, 25), sim_fun = sim_fun)

#-------------------------------------------------------------------------------
# export results
out <- list(
  # ecrt0 = ecrt0,
  # ecrt0.01 = ecrt0.01,
  ecrt0.05 = ecrt0.05,
  # ecrt0.1 = ecrt0.1,
  oracle = oracle,
  rmle = rmle,
  rmlep = rmlep,
  pvals = pvals,
  p = p,
  n = n,
  betap = betap,
  id = id
)

save(list = "out", file = paste0("simulation_1_", id, ".rda"))
