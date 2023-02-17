m <- 880 * 11
results_logistic_1 <- vector("list", m)
fl <- list.files()
betap_count <- numeric(11)
betap <- seq(0, 1, 0.1)
for (k in seq_len(m)) {
  name <- paste0("simulation_1_", k, ".rda")
  if (!(name %in% fl)) next
  load(name)
  betap_k <- out$betap
  betap_k_pos <- which(betap_k == betap)
  if (betap_count[betap_k_pos] == 800) next
  betap_count[betap_k_pos] <- betap_count[betap_k_pos] + 1
  results_logistic_1[[k]] <- out
}
save(list = "results_logistic_1", file = "results_logistic_1.rda")