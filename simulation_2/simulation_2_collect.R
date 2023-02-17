m <- 880
results_logistic_2 <- vector("list", m)
fl <- list.files()
count <- 0
for (k in seq_len(m)) {
  if (count == 800) next
  name <- paste0(paste0("simulation_2_", k, ".rda"))
  if (!(name %in% fl)) next
  count <- count + 1
  load(name)
  results_logistic_2[[k]] <- out
}
save(list = "results_logistic_2", file = "results_logistic_2.rda")