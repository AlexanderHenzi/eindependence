#-------------------------------------------------------------------------------
# packages
library(tidyverse)
library(lubridate)
theme_set(theme_bw())

#-------------------------------------------------------------------------------
# load data
data <- read.csv(file = "data_application/2011-capitalbikeshare-tripdata.csv")

#-------------------------------------------------------------------------------
# formatting and pre-selection
data <- data %>%
  mutate(
    month = month(Start.date),
    day = day(Start.date),
    minute = hour(Start.date) * 60 + minute(Start.date),
    day_of_week = 
      ((4 + day + cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30))[month]) %% 7),
    holiday = (month == 9 & day == 5) |
      (month == 10 & day == 10) |
      (month == 10 & day == 31) |
      (month == 11 & (abs(day - 24L) <= 1L))
  ) %>%
  filter(day_of_week <= 4 & month >= 9 & month <= 11 & !holiday) %>%
  select(-holiday) %>%
  mutate(Duration = log(Duration))

#-------------------------------------------------------------------------------
# estimate mean and variance of log duration by kernel smoothing (and keep all
# other relevant variables)
duration_mean_by_time <- function(
    duration,
    minute,
    month,
    Member.type,
    End.date
  ) {
  h <- 20
  K <- 
    exp(-(12 * 60 - abs(abs(outer(minute, minute, "-")) - 12 * 60))^2 / 2 / h^2)
  K[, month == 11] <- 0 
  K_weights_total <- rowSums(K)
  Duration_mean <- c(K %*% duration / K_weights_total)
  Duration_var <- c(K %*% (duration^2) / K_weights_total - Duration_mean^2)
  list(
    tibble(
      Duration_mean = Duration_mean,
      Duration_var = Duration_var,
      weights_total = K_weights_total,
      month = month,
      duration = duration,
      Member.type = Member.type,
      End.date = End.date,
      minute = minute
    )
  )
}

df <- data %>%
  group_by(Start.station.number, End.station.number) %>%
  summarise(
    duration_by_time = duration_mean_by_time(
      Duration,
      minute,
      month,
      Member.type,
      End.date
    )
  ) %>% 
  unnest(duration_by_time) %>%
  ungroup()

#-------------------------------------------------------------------------------
# extract validation data and run tests
df1 <- df %>%
  filter(month == 11 & weights_total >= 20) %>%
  arrange(as.POSIXct(End.date))

n0 <- 500
df_test <- data.frame(
  y = as.integer(df1$Member.type == "Member"),
  z = df1$Duration_mean,
  x = df1$duration
)

m <- nrow(df1)
e <- rep(1, m)
pred <- rep(1, m)
denom <- rep(1, m)
for (j in (n0+1):m) {
  df_fit <- df_test[seq_len(j - 1), ]
  df_new <- df_test[j, ]
  sim <- rnorm(500, mean = df1$Duration_mean[j], sd = sqrt(df1$Duration_var[j]))
  sim_data <- data.frame(x = sim, z = df_test$z[j])
  model <- glm(y ~ x + z, family = binomial, data = df_fit)
  pred[j] <- pmin(pmax(0.01, predict(model, df_new, type = "response")), 0.99)
  denom[j] <- mean(
    pmin(pmax(0.01, predict(model, sim_data, type = "response")), 0.99)
  )
  e[j] <- (1 - df_test$y[j] - pred[j]) / (1 - df_test$y[j] - denom[j])
}

which.max(cumprod(e) >= 1e4)
df1$End.date[which.max(cumprod(e) >= 1e4)]
prod(e)

# plot
plot_data <- 
  tibble(e = cumprod(e), date = as.POSIXct(df1$End.date), ind = seq_along(e))

Sys.setlocale("LC_TIME", "en_GB.UTF-8") 
case_study <- ggplot(plot_data) +
  geom_hline(yintercept = 1, lty = 3) +
  geom_step(aes(x = date, y = e)) +
  geom_rug(aes(x = date), size = 0.05) +
  scale_y_log10() +
  labs(x = "Date", y = "E-value")

pdf(width = 8, height = 4, file = "case_study.pdf")
case_study
dev.off()

# conditional randomization test:
cor_true <- abs(cor(df_test$y, df_test$x))
M <- 10000
cor_sim <- numeric(M)
for (j in seq_len(M)) {
  cor_sim[j] <- abs(cor(df_test$y, rnorm(nrow(df_test), 0, sqrt(df1$Duration_var))))
}

(sum(cor_sim >= cor_true) + 1) / (M + 1)
