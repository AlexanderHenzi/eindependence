library(tidyverse)
source("simulation_functions.R")
theme_set(theme_bw(base_size = 10))

# plots of misspecification functions in setting of CPT paper
set.seed(123)
n <- 5000
b <- rnorm(20)
z <- matrix(nrow = n, ncol = 20, rnorm(n * 20))
x <- rnorm(z %*% b)

mu <- z %*% b
dmu <- density(mu)
lines(dmu$x, dmu$y, col = 2)

mus <- sort(mu)

df_berrett <- data.frame(
  true = mus,
  quadratic_0.05 = mus + 0.05 * mus^2,
  quadratic_0.1 = mus + 0.1 * mus^2,
  quadratic_0.15 = mus + 0.15 * mus^2,
  cubic_0.01 = mus - 0.01 * mus^3,
  cubic_0.02 = mus - 0.02 * mus^3,
  cubic_0.03 = mus - 0.03 * mus^3,
  tanh_0.3 = tanh(0.3 * mus) / 0.3,
  tanh_0.6 = tanh(0.6 * mus) / 0.6,
  tanh_0.9 = tanh(0.9 * mus) / 0.9,
  paper = "Berrett et al."
)
mus_berrett <- mus

# plots of misspecification function in our simulations
p <- 4
Sigma <- toeplitz(1 / seq_len(p) * c(1, 1))
z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
x <- z[, p]
z <- z[, -p]
b <- solve(Sigma[-p, -p], Sigma[p, -p])
sdev <- 1 - Sigma[p, -p] %*% solve(Sigma[-p, -p], Sigma[p, -p])

mus <- sort(z %*% b)
dmu <- density(mus)

df_ours <- data.frame(
  true = mus,
  quadratic_0.4 = mus + 0.4 * mus^2,
  quadratic_0.8 = mus + 0.8 * mus^2,
  quadratic_1.2 = mus + 1.2 * mus^2,
  cubic_0.67 = mus - 2/3 * mus^3,
  cubic_1.33 = mus - 4/3 * mus^3,
  cubic_2 = mus - 2 * mus^3,
  tanh_2 = tanh(2 * mus) / 2,
  tanh_4 = tanh(4 * mus) / 4,
  tanh_6 = tanh(6 * mus) / 6,
  paper = "ours"
)
mus_ours <- mus

df_berrett <- df_berrett %>%
  pivot_longer(
    cols = -c(paper, true),
    names_to = c("misspec", "theta"),
    names_sep = "_",
    values_to = "mu"
  )

df <- df_ours %>%
  pivot_longer(
    cols = -c(paper, true),
    names_to = c("misspec", "theta"),
    names_sep = "_",
    values_to = "mu"
  ) %>%
  bind_rows(df_berrett)

plot_berrett <- df %>%
  filter(paper == "Berrett et al.") %>%
  group_by(misspec) %>%
  mutate(rank = factor(rank(theta, ties.method = "min"))) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, lty = 3) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_density(
    data = data.frame(
      mu = rep(mus_berrett, 3),
      misspec = rep(c("cubic", "quadratic", "tanh"), each = n)
    ),
    mapping = aes(x = mu, y = ..density..*75),
    lty = 5,
    color = "darkgray"
  ) +
  geom_line(aes(x = true, y = mu, color = rank, group = rank)) +
  ggthemes::scale_color_colorblind() +
  facet_grid(cols = vars(misspec)) +
  coord_cartesian(ylim = range(mus_berrett)) +
  theme(legend.position = "none") +
  labs(
    x = element_blank(),
    y = "Misspecified mean"
  )

plot_ours <- df %>%
  filter(paper == "ours") %>%
  group_by(misspec) %>%
  mutate(
    rank = as.integer(factor(theta)),
    rank = factor(
      rank,
      levels = 1:3,
      labels = c(
        "small (0.01/0.05/0.33 and 0.67/0.40/2.00)",
        "medium (0.02/0.10/0.67 and 0.67/0.80/4.00)",
        "large (0.03/0.15/1.00 and 2.00/1.20/6.00)"
      ),
      ordered = TRUE
    )
  ) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, lty = 3) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_density(
    data = data.frame(
      mu = rep(mus_ours, 3),
      misspec = rep(c("cubic", "quadratic", "tanh"), each = n)
    ),
    mapping = aes(x = mu, y = ..density..),
    lty = 5,
    color = "darkgray"
  ) +
  geom_line(aes(x = true, y = mu, color = rank, group = rank)) +
  ggthemes::scale_color_colorblind() +
  facet_grid(cols = vars(misspec)) +
  coord_cartesian(ylim = range(mus_ours)) +
  theme(legend.position = "bottom") +
  labs(
    x = "True mean of X conditional on Z",
    y = "Misspecified mean",
    color = expression(theta)
  )

pdf(file = "misspecification_function.pdf", width = 8, height = 6)
ggpubr::ggarrange(plot_berrett, plot_ours, nrow = 2, heights = c(1, 1.2))
dev.off()
