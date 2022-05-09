# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# ggplot theme settings
theme_set(theme_bw())

# Load generated data
load("Simulation1.rda")

# compute p-values and e-values
pvals <- numeric(length(cs))
for (i in seq_len(length(cs))) {
  ps <- numeric(nrep)
  for (j in seq_len(nrep)) {
    ps[j] <- (1 + sum(cpt_cor[j, i, ] >= true_cor[i, j])) / (M + 1)
  }
  pvals[i] <- mean(ps <= alpha)
}

get_evalue <- function(t_c, cpt_c, alpha, f = NULL) {
  nc <- nrow(t_c)
  nr <- ncol(t_c)
  out <- matrix(nrow = nc, ncol = 3)
  if (is.null(f)) {
    f <- function(obs, sample) obs / mean(sample)
  }
  for (i in seq_len(nc)) {
    e_mean <- numeric(nr)
    for (j in seq_len(nr)) {
      e_mean[j] <- f(t_c[i, j], cpt_c[j, i, ])
    }
    out[i, 1] <- mean(e_mean)
    out[i, 2] <- mean(e_mean >= 1/alpha)
    out[i, 3] <- mean(log(e_mean))
  }
  out
}
g <- function(obs, sample) {
  M <- length(sample)
  p <- (1 + sum(sample >= obs)) / (length(sample) + 1)
  1/sqrt(p) - 1
}
# Approach 1: Naive e-value T(X,Y,Z) / E_{X ~ Q(.|Z)}[T(X,Y,Z)] and p-to-E-value
evals <- get_evalue(true_dens, cpt_dens, alpha)
# evals_cod <- get_evalue(exp(true_codec), exp(cpt_codec), alpha)
evals_pe <- get_evalue(true_cor, cpt_cor, alpha, g)

a1_rr <- tibble(c = cs, CPT = pvals,  `LR E-value (1/E)` = evals[, 2], 
                                      # `FOCI E-value (1/E)` = evals_cod[, 2],  
                                      `1/sqrt(p[CPT]) - 1` = evals_pe[, 2]) %>%
  gather(key = "method", value = "rr", -c) %>%
  mutate(
    method = factor(
      method,
      levels = c("LR E-value (1/E)", 
                 # "FOCI E-value (1/E)", 
                 "1/sqrt(p[CPT]) - 1", 
                 "CPT"),
      labels = c("LR E-value (1/E)", 
                 #"FOCI E-value (1/E)", 
                 "1/sqrt(p[CPT]) - 1", 
                 "CPT"),
      ordered = TRUE
    )
  ) %>%
  ggplot() +
  geom_point(aes(x = c, y = rr, shape = method, group = method)) +
  geom_line(aes(x = c, y = rr, group = method, linetype = method)) +
  geom_hline(yintercept = alpha, lty = 3) +
  theme(legend.position = c(0.2, 0.6), legend.title = element_blank()) +
  labs(
    x = expression("Parameter"~ italic("c")),
    y = "Rejection rate"
  ) +
  ggtitle("(a) Rejection rates of CPT and density-based e-value")


a1_mean <- tibble(
  c = cs,
  `LR E-value (1/E)` = evals[, 1],
  #`FOCI E-value (1/E)` = evals_cod[, 1],
  `1/sqrt(p[CPT]) - 1` = evals_pe[, 1]
) %>%
  gather(key = "method", value = "e", -c) %>%
  mutate(
    method = factor(
      method,
      levels = c("LR E-value (1/E)", 
                 #"FOCI E-value (1/E)", 
                 "1/sqrt(p[CPT]) - 1"),
      labels = c("LR E-value (1/E)", 
                 #"FOCI E-value (1/E)", 
                 "1/sqrt(p[CPT]) - 1"),
      ordered = TRUE
    )
  ) %>%
  filter(is.finite(e)) %>%
  ggplot() +
  geom_point(aes(x = c, y = e, group = method, shape = method)) +
  geom_line(aes(x = c, y = e, group = method, linetype = method)) +
  labs(
    x = expression("Parameter"~ italic("c")),
    y = "E-value (mean)"
  ) +
  geom_hline(yintercept = 1, lty = 3) +
  ggtitle("(b) Mean of transformed e-value") +
  theme(legend.position = c(0.1, 1e11), legend.title = element_blank()) +
  scale_y_log10()

a1_log_mean <- tibble(
  c = cs,
  `LR E-value (1/E)` = evals[, 3],
  #`FOCI E-value (1/E)` = evals_cod[, 3],
  `1/sqrt(p[CPT]) - 1` = evals_pe[, 3]
) %>%
  gather(key = "method", value = "e", -c) %>%
  mutate(
    method = factor(
      method,
      levels = c("LR E-value (1/E)", 
                 #"FOCI E-value (1/E)", 
                 "1/sqrt(p[CPT]) - 1"),
      labels = c("LR E-value (1/E)", 
                 #"FOCI E-value (1/E)", 
                 "1/sqrt(p[CPT]) - 1"),
      ordered = TRUE
    )
  )%>%
  filter(is.finite(e)) %>%
  ggplot() +
  geom_point(aes(x = c, y = e, group = method, shape = method)) +
  geom_line(aes(x = c, y = e, group = method, linetype = method)) +
  labs(
    x = expression("Parameter"~ italic("c")),
    y = "E-value (mean log)"
  ) +
  ggtitle("(c) Mean of log of transformed e-value")  +
  theme(legend.position = c(0.1, 4), legend.title = element_blank())

a1_plot <- ggarrange(
  a1_rr,
  ggarrange(a1_mean, a1_log_mean, ncol = 2),
  nrow = 2
)

pdf(width = 8, height = 6, file = "plots.pdf")
print(a1_plot)
dev.off()
