args = commandArgs(trailingOnly = TRUE)
# --------------------------------------------------------
# argument(s) to this script:
#   1) fileName; where to save the .rda file.
#      DO NOT INCLUDE THE EXTENSION!
# example usage:
# Rscript simulation_1_collect.R "simulation_1"

if(length(args) == 0){
  fileName = "simulation_2"
} else if (length(args) > 1) {
  stop("Invalid number of arguments.")
} else {
  fileName = toString(args[1])
}

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
  unlink(name)
}
save(
  list = c("results_logistic_2", "correlation", "p"),
  file = paste(fileName, "rda", sep=".")
)
