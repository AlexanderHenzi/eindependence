#-------------------------------------------------------------------------------
# functions and packages
source("simulation_functions.R")
require(MASS)

#-------------------------------------------------------------------------------
# parameters for different settings

## dimension; the following values were tested
##     q = 4
##     q = 8
q <- 4

## correlation (positive/negative correlations); set to
##    "pos" for all correlations positive
##     "negative_cor" for positive and negative correlations
correlation <- "pos"

## not_penalize_rmle; do (or not) penalize the parameter of interest in RMLE
##     FALSE penalize all parameters (default setting in glmnet)
##     TRUE do not penalize the parameter of interest (the first one, for x)
not_penalize_rmle <- TRUE

# check for validity of parameter values
if (!is.numeric(q) || length(q) != 1)
  stop("invalid p")
if (!(identical(correlation, "pos") || identical(correlation, "negative_cor")))
  stop("invalid correlation")
if (!(isTRUE(not_penalize_rmle) || isFALSE(not_penalize_rmle)))
  stop("invalid not_penalize_rmle")

# ------------------------------------------------------------------------------
## parse to parameter names and values used below
p <- q
if (not_penalize_rmle) {
  not_penalize_rmle <- 1
} else {
  not_penalize_rmle <- NULL
}
if (correlation == "pos") {
  correlation <- c(1, 1)
} else {
  correlation <- c(1, -1)
}

## other fixed or derived parameters, do not change
n <- 2000
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
betap <- seq(0, 1, 0.1)
betap <- betap[id %% length(betap) + 1]
n0 <- 5 * p

#-------------------------------------------------------------------------------
# data generation
set.seed(id)
Mu <- rep(0, p)
Sigma <- toeplitz(1 / seq_len(p) * correlation)
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
ecrt0 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0)
ecrt0.05 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0.05)
ecrt0.01 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0.01)
ecrt0.1 <- e_crt_logistic(x = x, y = y, z = z, sim_fun = sim_fun, n0 = n0, eps = 0.1)

## oracle conditional randomization
oracle <- e_crt_logistic_oracle(x = x, y = y, z = z, sim_fun = sim_fun, par = c(intercept, pars, betap))

## running mle
rmle <- running_mle_logistic(x = x, y = y, z = z, n0 = n0)

## penalized running mle
rmlep <- running_mle_logistic_penalized(x = x, y = y, z = z, n0 = n0, not_penalize = not_penalize_rmle)

## p-values
pvals <- crt_lrt_logistic(x = x, y = y, z = z, stops = seq(25, n, 25), sim_fun = sim_fun)

#-------------------------------------------------------------------------------
# export results
out <- c(
  ecrt0 = ecrt0,
  ecrt0.01 = ecrt0.01,
  ecrt0.05 = ecrt0.05,
  ecrt0.1 = ecrt0.1,
  oracle = oracle,
  rmle = rmle,
  rmlep = rmlep,
  pvals = pvals,
  p = p,
  n = n,
  betap = betap,
  id = id
)

# also export parameters
save(
  list = c(
    "out",
    "p",
    "correlation",
    "not_penalize_rmle"
  ),
  file = paste0("simulation_1_", id, ".rda")
)
