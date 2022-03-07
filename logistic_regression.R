# functions
expit <- function(x) 1/(1 + exp(-x))

# setup simulation parameters
b0 <- -1
b1 <- 0.5
b2 <- seq(0, 3, 0.5)
alpha <- 0.05
ee <- pp <- numeric(length(b2))
alpha <- 0.05
eps <- 0.01 # truncate predictions at [eps, 1-eps] to avoid zero e-values
n <- 50
m <- 500 # number of simulations
start_at <- 10 # start parameter estimation at the "start_at"th observation

# run simulations for all values of beta_2

for (k in seq_along(b2)) {
  e <- numeric(m)
  p <- numeric(m)
  
  # run simulations
  pb <- txtProgressBar(max = m)
  for (i in seq_len(m)) {
    setTxtProgressBar(pb, i)
    Z <- rnorm(n)
    X <- rbinom(n, 1, pnorm(Z))
    Y <- rbinom(n, 1, expit(b0 + b1*Z + b2[k]*X))
    e_tmp <- rep(1, n)
    for (j in start_at:n) {
      y <- Y[seq_len(j - 1)]
      x <- X[seq_len(j - 1)]
      z <- Z[seq_len(j - 1)]
      model <- glm(y ~ z + x, family = binomial())
      pred_probs <- predict(
        model,
        data.frame(x = c(0, 1), z = Z[j]), # first entry: X = 0, second entry: X = 1
        type = "response"
      )
      pred_probs <- pmin(pmax(pred_probs, eps), 1 - eps)
      alternative <- if (X[j] == 0) pred_probs[1] else pred_probs[2]
      err <- runif(1, -0.1, 0.1) # set to 0 for exact distribution
      randomize_prob <- pmax(0, pmin(1, pnorm(Z[j] + err))) 
      denom <- (1 - randomize_prob) * (pred_probs[1]) + randomize_prob * pred_probs[2]
      e_tmp[j] <- if (Y[j] == 1) alternative / denom else 
        (1 - alternative) / (1 - denom)
    }
    e[i] <- max(cumprod((e_tmp)))
    model <- glm(Y ~ Z + X, family = binomial())
    p[i] <- summary(model)$coefficients[3, 4]
  }
  close(pb)
  ee[k] <- mean(e > 1/alpha)
  pp[k] <- mean(p < alpha)
}

# plot results
pdf("logistic_regression2.pdf", width = 8, height = 6)
plot(
  b2,
  ee,
  type = "b",
  pch = 16,
  xlab = expression(beta[2]),
  ylab = "Frequency of rejection",
  main = "Rejection rates",
  ylim = c(0, 1)
)
points(b2, pp, type = "b", col = 2, pch = 16)
abline(h = c(alpha, 1), col = "darkgray", lty = 5)
legend("topleft", col = c(1, 2), lty = c(1, 1), legend = c("p-value", "e-value"))
dev.off()
