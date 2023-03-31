#-------------------------------------------------------------------------------
# packages and global ggplot settings
library(tidyverse)
library(ggthemes)
library(ldbounds)
library(ggpubr)

#-------------------------------------------------------------------------------
# parameters
# argument(s) to this script represent the following:
#   1) path; the path to the collected output of simulation_2;
# parameter must be passed, e.g.:
# Rscript simulation_2_plots.R "simulation_2.rda"

path <- toString(args[1])

#-------------------------------------------------------------------------------
## load functions and data
source("simulation_functions.R")
load(path)

## parse parameter values for output file name
if (p == 8) {
  dimension <- "_q8"
} else {
  dimension <- ""
}
if (all.equal(correlation, c(1, -1))) {
  correlation <- "_negative_cor"
} else {
  correlation <- ""
}

nulls <- sapply(results_logistic_2, is.null)
sum(nulls)
results_logistic_2 <- results_logistic_2[!nulls]

n_crt <- c(200, 400, 1000, 1400, 2000)
alphas <- c(0.01, 0.05, 0.1)
for (j in seq_along(results_logistic_2)) {
  tmp <- results_logistic_2[[j]]
  for (i in seq_along(tmp)) {
    tmpi <- tmp[[i]]
    tmpi[[1]] <- get_first_rejection(
      cumprod(do.call(get_e, c(as.list(tmpi[[1]]), M = 500, unbiased = TRUE))),
      alphas = alphas,
      tmax = Inf
    )
    for (k in seq_along(n_crt)) {
      tmpi[[1]][paste0("crt_", n_crt[k])] <- 
        (tmpi[[2]]$p_crt[tmpi[[2]]$stops == n_crt[k]] <= alphas)
      tmpi[[1]][paste0("crt_cor_", n_crt[k])] <- 
        (tmpi[[2]]$p_crt_cor[tmpi[[2]]$stops == n_crt[k]] <= alphas)
    }
    tmpi <- tmpi[-2]
    tmpi <- unnest(as_tibble(tmpi), ecrt)
    tmp[[i]] <- tmpi
  }
  results_logistic_2[[j]] <- do.call(rbind, tmp)
}
results_logistic_2 <- do.call(rbind, results_logistic_2)
results_logistic_2 <- results_logistic_2 %>%
  rename(evalues = rejected) %>%
  gather(key = "method", value = "rejected", evalues, starts_with("crt_"))

misspecified_model_data <- results_logistic_2 %>%
  filter(alphas == 0.05 & is.infinite(n_sample)) %>%
  filter(method %in% c(
    "evalues",
    "crt_200",
    "crt_1000",
    "crt_2000",
    "crt_cor_200",
    "crt_cor_2000")
  ) %>%
  mutate(
    method = factor(
      method,
      labels = c("E-value", "CRT (200)", "CRT (1000)", "CRT (2000)", "CRT (200; cor)", "CRT (2000; cor)"),
      levels = c("evalues", "crt_200", "crt_1000", "crt_2000", "crt_cor_200", "crt_cor_2000"),
      ordered = TRUE
    )
  ) %>%
  group_by(alphas, theta, misspec, n_sample, method, update_reuse) %>%
  summarise(rejected = mean(rejected))

theme_set(theme_bw(base_size = 13.5))
misspecified_model <- ggplot() +
  geom_hline(yintercept = 0.05, lty = 3) +
  geom_line(
    data = misspecified_model_data,
    mapping = aes(x = theta, y = rejected, color = method, group = method)
  ) +
  geom_point(
    data = filter(
      group_by(misspecified_model_data, alphas, misspec, method),
      rank(theta) %in% seq(5, 40, 5)
    ),
    aes(x = theta, y = rejected, color = method, shape = method)
  ) +
  facet_grid(cols = vars(misspec), scales = "free_x") +
  coord_cartesian(ylim = c(0, 0.25)) +
  theme(legend.position = "none") +
  labs(
    x = expression(xi),
    y = "Rejection rate",
    color = element_blank(),
    shape = element_blank()
  ) +
  scale_color_colorblind() +
  ggtitle("(a)")

finite_sample <- results_logistic_2 %>%
  filter(alphas == 0.05 & is.finite(n_sample) & !update_reuse) %>%
  filter(method %in% c("evalues", "crt_200", "crt_1000", "crt_2000", "crt_cor_200", "crt_cor_2000")) %>%
  mutate(
    method = factor(
      method,
      labels = c("E-value", "CRT (200)", "CRT (1000)", "CRT (2000)", "CRT (200; cor)", "CRT (2000; cor)"),
      levels = c("evalues", "crt_200", "crt_1000", "crt_2000", "crt_cor_200", "crt_cor_2000"),
      ordered = TRUE
    )
  ) %>%
  group_by(alphas, theta, misspec, n_sample, method, update_reuse) %>%
  summarise(rejected = mean(rejected)) %>%
  ggplot() +
  geom_hline(yintercept = 0.05, lty = 3) +
  geom_line(aes(x = n_sample, y = rejected, color = method, group = method)) +
  geom_point(aes(x = n_sample, y = rejected, color = method, shape = method)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme(legend.position = "none") +
  labs(
    x = "Unlabelled sample size",
    y = "Rejection rate",
    color = element_blank(),
    shape = element_blank()
  ) +
  scale_color_colorblind() +
  ggtitle("(b)")

recycle_data <- results_logistic_2 %>%
  filter(alphas == 0.05 & is.finite(n_sample) & update_reuse) %>%
  filter(method %in% c("evalues", "crt_200", "crt_1000", "crt_2000", "crt_cor_200", "crt_cor_2000")) %>%
  group_by(alphas, theta, misspec, n_sample, method, update_reuse) %>%
  summarise(rejected = mean(rejected)) %>%
  spread(key = "method", value = "rejected") %>%
  ungroup() %>%
  mutate_at(
    .vars = vars(contains("crt_")),
    .funs = function(x) mean(x)
  ) %>%
  gather(key = "method", value = "rejected", evalues, contains("crt_")) %>%
  mutate(
    method = factor(
      method,
      labels = c("E-value", "CRT (200)", "CRT (1000)", "CRT (2000)", "CRT (200; cor)", "CRT (2000; cor)"),
      levels = c("evalues", "crt_200", "crt_1000", "crt_2000", "crt_cor_200", "crt_cor_2000"),
      ordered = TRUE
    )
  ) %>%
  ggplot() +
  geom_hline(yintercept = 0.05, lty = 3) +
  geom_line(aes(x = n_sample, y = rejected, color = method, group = method)) +
  geom_point(aes(x = n_sample, y = rejected, color = method, shape = method)) +
  coord_cartesian(ylim = c(0, 0.16)) +
  theme(legend.position = "right") +
  labs(
    x = "Unlabelled sample size",
    y = "Rejection rate",
    color = element_blank(),
    shape = element_blank()
  ) +
  scale_color_colorblind() +
  ggtitle("(c)")

robustness_plot <- ggarrange(
  misspecified_model,
  ggarrange(finite_sample, recycle_data, ncol = 2, widths = c(1, 1.3)),
  nrow = 2
)

plot_name <- paste0(
  "simulation_robustness_",
  dimension,
  correlation,
  ".pdf"
)
pdf(width = 8, height = 6, file = plot_name)
print(robustness_plot)
dev.off()