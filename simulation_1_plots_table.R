#-------------------------------------------------------------------------------
# packages and global ggplot settings
library(tidyverse)
library(ggthemes)
library(ldbounds)
library(ggpubr)

args = commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Please specify an epsilon and path")
}
#-------------------------------------------------------------------------------
# sample size for different type I error and power (with results of simulation_1.R)

# parameters -------------------------------------------------------------------
# arguments to this script represent the following parameters respectively:
#   1) eps; the truncation level of the e-crt;
#      values in {0,0.1,0.05,0.1};
#   2) path; the path to the collected output of simulation_1;
# both arguments must be passed, e.g.:
# Rscript simulation_1_plots_table.R 0.05 "simulation_1.rda"

eps_vals = c(0, 0.01, 0.05, 0.1)
eps <- as.numeric(args[1])
if (!eps %in% eps_vals) {
  stop("Invalid value for parameter epsilon")
}

path <- toString(args[2])

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
if (identical(correlation, c(1, -1))) {
  corr_par <- "_negative_cor"
} else {
  corr_par <- ""
}
if (is.null(not_penalize_rmle)) {
  penalization <- ""
} else {
  penalization <- "_not_penalize"
}

tmax = 2000

## remove null entries (a higher number than 800 is generated, to make sure than
## at least 800 simulations are available when running the computations on the
## cluster and some of the tasks fail due to time/memory limit)
nulls <- sapply(results_logistic_1, is.null)
sum(nulls)
results_logistic_1 <- results_logistic_1[!nulls]

eps_remove <- (1:4)[-match(eps, c(0, 0.01, 0.05, 0.1))]
results_logistic_1 <- lapply(
  results_logistic_1,
  function(x) x[-eps_remove]
)

## the type I level errors to consider; containers for rejection rates of CRT 
## and LRT (CRT and LRT are computed after each new batch of 25 observations,
## with a total of 80 batches).
alphas <- c(0.01, 0.05)
rejections_crt <- rejections_asymptotic <- 
  array(dim = c(tmax/25, length(alphas), 11), 0)

## compute e-values and rejection times, for each simulation
for (k in seq_along(results_logistic_1)) {
  ## extract k-th simulation; store parameters; compute e-values and rejection
  ## times for all methods
  tmp <- results_logistic_1[[k]]
  betap <- tmp$betap
  id <- tmp$id
  rejections_k <- vector("list", 6)
  for (j in 1:2) {
    e <- cumprod(do.call(get_e, c(as.list(tmp[[j]]), M = 500, unbiased = TRUE)))
    rejections <- get_first_rejection(e = e, alphas = alphas, tmax = tmax)
    rejections$method <- names(tmp)[j]
    rejections_k[[j]] <- rejections
  }
  for (j in 3:4) {
    e <- tmp[[j]]$M
    rejections <- get_first_rejection(e = e, alphas = alphas, tmax = tmax)
    rejections$method <- names(tmp)[j]
    rejections_k[[j]] <- rejections
  }
  
  rejections_k <- do.call(rbind, rejections_k[1:4])
  rejections_k$betap <- betap
  rejections_k$id <- id
  results_logistic_1[[k]] <- rejections_k
  
  ## compute rejection times for CRT and LRT
  ind <- which.min(abs(seq(0, 1, 0.1) - betap))
  rejections_asymptotic[, , ind] <- rejections_asymptotic[, , ind] + 
    outer(tmp[[5]]$p_asymptotic, alphas, "<=")
  rejections_crt[, , ind] <- 
    rejections_crt[, , ind] + outer(tmp[[5]]$p_crt, alphas, "<=")
}
results_logistic_1 <- do.call(rbind, results_logistic_1)

## rejections under the null hypothesis (table; all methods except group
## sequential, which are computed below)
table_1 <- results_logistic_1 %>%
  filter(betap == 0) %>%
  group_by(method, alphas) %>%
  summarise(rejected = mean(rejected))

## compute sample sizes required for rejecting at given level; "betas" gives
## the desired power (1 - type II error).

betas <- c(0.8, 0.95)
df <- results_logistic_1 %>%
  group_by(method, alphas, betap) %>%
  nest() %>%
  mutate(
    results = map(
      .x = data,
      .f = function(dat) {
        lapply(
          betas,
          function(b) find_sample_size(dat$size, dat$rejected, b, tmax)
        )
      }
    )
  ) %>%
  unnest(cols = results) %>%
  unnest(cols = results)

## for LRT and CRT
nsamples = length(results_logistic_1) / 11
rejections_asymptotic <- rejections_asymptotic / nsamples
rejections_crt <- rejections_crt / nsamples
df_asymptotic <- vector("list", length(betas))
df_crt <- vector("list", length(betas))
for (k in seq_along(betas)) {
  ## for LRT, CRT: The "minimum" sample size is searched over the sizes
  ## 25, 50, ..., 1975, 2000 which are actually available. In the "worst case",
  ## the actual minimum might be achieved at size 26, 51, ..., so we subtract
  ## 25 from the result to make sure that the minimum sample size is not
  ## over-estimated.
  tmp <- t(apply(
    rejections_asymptotic,
    c(2, 3),
    function(x) pmax(1, min(which(x >= betas[k])) - 1) * 25
  ))
  tmp <- cbind(betap = seq(0, 1, 0.1), tmp)
  colnames(tmp) <- c("betap", alphas)
  tmp <- as_tibble(tmp)
  tmp <- gather(tmp, key = "alphas", value = "worst", -betap)
  tmp$alphas <- as.numeric(tmp$alphas)
  tmp$method <- "asymptotic"
  tmp$data <- vector("list", 1)
  tmp$beta <- betas[k]
  tmp$power_achieved <- is.finite(tmp$worst)
  tmp$worst <- pmin(tmax, tmp$worst)
  tmp$average <- tmp$worst
  df_asymptotic[[k]] <- tmp
  
  ## for CRT
  tmp <- t(apply(
    rejections_crt,
    c(2, 3),
    function(x) pmax(1, min(which(x >= betas[k])) - 1) * 25
  ))
  tmp <- cbind(betap = seq(0, 1, 0.1), tmp)
  colnames(tmp) <- c("betap", alphas)
  tmp <- as_tibble(tmp)
  tmp <- gather(tmp, key = "alphas", value = "worst", -betap)
  tmp$alphas <- as.numeric(tmp$alphas)
  tmp$method <- "crt"
  tmp$data <- vector("list", 1)
  tmp$beta <- betas[k]
  tmp$power_achieved <- is.finite(tmp$worst)
  tmp$worst <- pmin(tmax, tmp$worst)
  tmp$average <- tmp$worst
  df_crt[[k]] <- tmp
}
df_asymptotic <- do.call(rbind, df_asymptotic)
df_crt <- do.call(rbind, df_crt)

## combine all results
df <- rbind(df,  df_crt, df_asymptotic) %>%
  gather(key = "type", value = "size", average, worst) %>%
  ungroup() %>%
  mutate(
    type = factor(
      type,
      levels = c("worst", "average"),
      labels = c("Worst case", "Average")
    ),
    method = factor(
      method,
      levels = c("asymptotic", "crt", "oracle", paste0("ecrt", eps), "rmle", "rmlep"),
      labels = c("LRT", "CRT", "E-CRT-O", "E-CRT", "R-MLE", "R-MLE-P"),
      ordered = TRUE
    )
  )


## generate plot
df <- df %>%
  mutate(
    lty = c("1","5")[1 + (type == "Average" & method %in% c("CRT", "LRT"))],
    type = (as.character(factor(
      type,
      levels = c("Worst case", "Average"),
      labels = c("N(beta, eta)", "N[av](beta,eta)")
    )))
  )

theme_set(theme_bw(base_size = 13.5))
for (a in c(0.01, 0.05)) { # for both alphas (0.01, 0.05)
  simulation_alternative <- ggplot() +
    geom_point(
      data = filter(df, power_achieved & alphas == a),
      aes(
        x = betap,
        y = size,
        color = method,
        group = method,
        shape = method,
      ),
      cex = 2
    ) +
    geom_line(
      data = filter(df, power_achieved & alphas == a),
      aes(
        x = betap,
        y = size,
        color = method,
        group = method,
        linetype = lty
      ),
      lwd = 0.5
    ) +
    ggthemes::scale_color_colorblind() +
    scale_linetype_manual(values = c(1, 5)) +
    facet_grid(rows = vars(type), cols = vars(beta), labeller = label_parsed) +
    labs(
      x = expression(beta),
      y = "Sample size",
      linetype = element_blank(),
      color = element_blank(),
      shape = element_blank()
    ) +
    theme(legend.position = "right", legend.text = element_text(size = 12)) +
    guides(
      color = guide_legend(ncol = 1, byrow = TRUE),
      linetype = "none"
    )
  
  plot_name <- paste0(
    "simulation_alternative_eps",
    sub("0.", "",toString(eps)),
    dimension,
    corr_par,
    penalization,
    "_alpha_",
    a,
    ".pdf"
  )

  pdf(width = 8, height = 4, file = plot_name)
  print(simulation_alternative)
  dev.off()
}


#-------------------------------------------------------------------------------
# comparison with group sequential methods

## load functions
# rm(list = ls())
source("simulation_functions.R")

## reload data
load(path)

## remove null entries (see above)
nulls <- sapply(results_logistic_1, is.null)
sum(nulls)
results_logistic_1 <- results_logistic_1[!nulls]

## select "epsilon" for the E-CRT (as above)
eps_remove <- (1:4)[-match(eps, c(0, 0.01, 0.05, 0.1))]
results_logistic_1 <- lapply(
  results_logistic_1,
  function(x) x[-eps_remove]
)

## fix type I error level, maximum sample size to be considered, and the number
## of stops for group sequential methods
alphas <- c(0.01, 0.05)
nstops <- c(5, 10, 20, 40)
stop_seq <- 
  lapply(nstops, function(x) seq(tmax / 25 /x, tmax / 25, tmax / 25 / x))

## compute the nominal alpha levels for the group sequential methods (Pocock and
## O'Brien & Fleming method)
alphas_pk <- lapply(
  nstops,
  function(k) {
    do.call(
      rbind,
      lapply(
        alphas,
        function(a) commonbounds(looks = k, iuse = "PK", alpha = a)$nom.alpha
      )
    )
  }
)
alphas_of <- lapply(
  nstops,
  function(k) {
    do.call(
      rbind,
      lapply(
        alphas,
        function(a) commonbounds(looks = k, iuse = "OF", alpha = a)$nom.alpha
      )
    )
  }
)
variants <- length(alphas) * length(nstops)

## find first minimum rejection time of e-variables and group sequential
## methods, for all variants (combination of type I error and number of stops)
for (k in seq_along(results_logistic_1)) {
  tmp <- results_logistic_1[[k]]
  betap <- tmp$betap
  id <- tmp$id
  rejections_k <- vector("list", 2 * variants + 1)
  e <- cumprod(do.call(get_e, c(as.list(tmp[[1]][seq_len(tmax), ]), M = 500, unbiased = TRUE)))
  rejections <- get_first_rejection(e = e, alphas = alphas, tmax = tmax)
  rejections$method <- names(tmp)[1]
  rejections$betap <- betap
  rejections_k[[1]] <- rejections
  
  l <- 2
  for (i in seq_along(alphas)) {
    for (j in seq_along(nstops)) {
      rejections <- get_first_rejection(
        e = 1 / tmp$pvals$p_asymptotic[stop_seq[[j]]],
        alphas = alphas_pk[[j]][i, ],
        tmax = nstops[j]
      )
      rejections$size <- (tmax / nstops[j]) * rejections$size
      rejections$method <- paste0("pk_", nstops[j])
      rejections$alphas <- alphas[i]
      rejections$betap <- betap
      rejections_k[[l]] <- rejections
      
      l <- l + 1
      rejections <- get_first_rejection(
        e = 1 / tmp$pvals$p_asymptotic[stop_seq[[j]]],
        alphas = alphas_of[[j]][i, ],
        tmax = nstops[j]
      )
      rejections$size <- (tmax / nstops[j]) * rejections$size
      rejections$method <- paste0("of_", nstops[j])
      rejections$alphas <- alphas[i]
      rejections$betap <- betap
      rejections_k[[l]] <- rejections
      l <- l + 1
    }
  }
  rejections_k <- do.call(rbind, rejections_k)
  results_logistic_1[[k]] <- rejections_k
}
results_logistic_1 <- do.call(rbind, results_logistic_1)

## get rejection rates under null hypothesis (for the table, for group 
## sequential methods)
table_1 <- results_logistic_1 %>%
  filter(betap == 0) %>%
  group_by(method, alphas) %>%
  summarise(rejected = mean(rejected)) %>%
  filter(grepl("20", method)) %>%
  bind_rows(table_1) %>%
  arrange(alphas, method)

## export table of rejection rates under the null
filename <- paste0(
  "table_1_eps_",
  sub("0.", "", toString(eps)),
  dimension,
  corr_par,
  penalization,
  ".txt"
)
write.table(table_1, row.names = FALSE, file = filename)

## plot
theme_set(theme_bw(base_size = 12))
grpseq <- results_logistic_1 %>%
  filter(alphas == 0.05 & betap > 0) %>%
  filter(method %in% c(paste0("ecrt", eps), "of_20", "pk_20")) %>%
  mutate(
    method = factor(
      method,
      levels = c(paste0("ecrt", eps), "pk_20", "of_20"),
      labels = c("E-CRT", "LRT-PK", "LRT-OF"),
      ordered = TRUE
    )
  ) %>%
  mutate(size = size + 1000 * (!rejected)) %>%
  ggplot() +
  geom_hline(yintercept = c(0, 1), color = "darkgray") +
  stat_ecdf(
    aes(x = size, color = method, group = method, linetype = method),
    lwd = 0.7
  ) +
  scale_linetype_manual(values = c("solid", "longdash", "twodash")) +
  facet_wrap(.~betap, nrow = 2) +
  coord_cartesian(xlim = c(0, tmax)) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 8),
    legend.text = element_text(size = 10)
  ) +
  scale_x_continuous(
    breaks = c(0, 1000, tmax)
  ) +
  ggthemes::scale_color_colorblind() +
  labs(
    x = "Sample size at rejection",
    y = "Empirical distribution function",
    color = element_blank(),
    linetype = element_blank()
  )  +
  guides(
    color = guide_legend(ncol = 1, byrow = TRUE)
  )

## export the figure
plot_name <- paste0(
  "simulation_group_sequential_eps",
  sub("0.", "", toString(eps)),
  dimension,
  corr_par,
  ".pdf"
)

pdf(height = 3.5, width = 8, file = plot_name)
print(grpseq)
dev.off()