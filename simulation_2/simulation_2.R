#-------------------------------------------------------------------------------
# functions and packages
source("simulation_functions.R")

#-------------------------------------------------------------------------------
# parameters for different settings

## dimension; the following values were tested
##     q = 4
##     q = 8
q <- 4

## correlation (positive/negative correlations); set to
##    "pos" for all correlations positive
##    "negative_cor" for positive and negative correlations
correlation <- "pos"

# check for validity of parameter values
if (!is.numeric(q) || length(q) != 1)
  stop("invalid p")
if (!(identical(correlation, "pos") || identical(correlation, "negative_cor")))
  stop("invalid correlation")

#-------------------------------------------------------------------------------
# other fixed or derived parameters, do not change
n <- 2000
p <- q
if (correlation == "pos") {
  correlation <- c(1, 1)
} else {
  correlation <- c(1, -1)
}
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
specs1 <- expand.grid(
  theta = seq(0, 2, 0.05),
  misspec = c("quadratic", "cubic", "tanh"),
  n_sample = Inf,
  update_reuse = FALSE,
  stringsAsFactors = FALSE
)
specs2 <- expand.grid(
  theta = 0,
  misspec = "none",
  n_sample = seq(50, 500, 50),
  update_reuse = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
specs <- rbind(specs1, specs2)
specs$theta[specs$misspec == "tanh"] <- specs$theta[specs$misspec == "tanh"] * 3
specs$theta[specs$misspec == "quadratic"] <- specs$theta[specs$misspec == "quadratic"] * 0.6
n0 <- 5 * p
betap <- 0
nsim <- 500

#-------------------------------------------------------------------------------
# loop over all parameter specifications
out <- vector("list", length = nrow(specs))
for (k in seq_len(nrow(specs))) {
  n_sample <- specs$n_sample[k]
  misspec <- specs$misspec[k]
  theta <- specs$theta[k]
  misspec_fun <- switch(
    misspec,
    none = identity,
    quadratic = function(x) x + theta * x^2,
    cubic = function(x) x - theta * x^3,
    tanh = if (theta == 0) identity else function(x) tanh(x * theta) / theta
  )
  update_reuse <- specs$update_reuse[k]
  
  #-----------------------------------------------------------------------------
  # data generation and computation of e-values (quadratic error)
  set.seed(k * 100 + id)

  Mu <- rep(0, p)
  Sigma <- toeplitz(1 / seq_len(p) * correlation)
  z <- MASS::mvrnorm(n = n, mu = Mu, Sigma = Sigma)
  x <- z[, p]
  z <- z[, -p]
  coefs <- runif(p - 1, -1, 1)
  intercept <- runif(1, -1, 1)
  probs <- c(expit(intercept + x * betap + z %*% coefs))
  y <- rbinom(n, 1, probs)
  
  if (is.infinite(n_sample)) {
    xz_sample <- NULL
    coef_vec <- c(
      solve(Sigma[-p, -p], Sigma[p, -p]),
      1 - Sigma[p, -p] %*% solve(Sigma[-p, -p], Sigma[p, -p])
    )
  } else {
    coef_vec <- NULL
    xz_sample <- MASS::mvrnorm(n = n_sample, mu = Mu, Sigma = Sigma)[, c(p, seq_len(p - 1))]
  }
  
  ecrt <- e_crt_logistic_gaussian(
    x = x,
    y = y,
    z = z,
    n0 = n0,
    eps = 0.05,
    coef = coef_vec,
    xz_sample = xz_sample,
    update_sample = update_reuse,
    bias_fun = misspec_fun,
    nsim = nsim
  )
  pvals <- crt_logistic_gaussian(
    x = x,
    y = y,
    z = z,
    stops = seq(200, n, 200),
    coef = coef_vec,
    xz_sample = xz_sample,
    reuse_data = update_reuse,
    bias_fun = misspec_fun,
    nsim = nsim
  )

  out[[k]] <- list(
    ecrt = ecrt,
    pvals = pvals,
    theta = theta,
    betap = betap,
    misspec = misspec,
    n_sample = n_sample,
    update_reuse = update_reuse,
    sim = id
  )
  cat(k, "\n")
}

#-------------------------------------------------------------------------------
# collect and export results
save(
  list = c("out", "p", "correlation"),
  file = paste0("simulation_2_", id, ".rda")
)

