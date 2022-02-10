# load packages
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# parallelization for CPR
cluster <- makeCluster(14)
registerDoParallel(cluster)

# ggplot theme settings
theme_set(theme_bw())


############# Implementation of the conditional permutation test ##################

generate_X_CPT_gaussian = function(nstep,M,X0,mu,sig2){
	# Runs the conditional permutation test using the distribution X | Z=Z[i] ~ N(mu[i],sig2[i])
	log_lik_mat = -(X0^2)%*%t(1/2/sig2) + X0%*%t(mu/sig2)
	# log_lik_mat[i,j] = density at X=X0[i] when Z=Z[j]
	Pi_mat = generate_X_CPT(nstep,M,log_lik_mat)
	X_mat = X0[Pi_mat]
	dim(X_mat) = c(M,length(X0))
	return(t(X_mat))
}

generate_X_CPT = function(nstep,M,log_lik_mat,Pi_init=NULL){
	# log_lik_mat is the n-by-n matrix with entries log(q(X_i|Z_j))
	# this function produces M exchangeable permutations, initialized with permutation Pi_init
	n = dim(log_lik_mat)[1]
	if(length(Pi_init)==0){
		Pi_init = 1:n
	}
	Pi_ = generate_X_CPT_MC(nstep,log_lik_mat,Pi_init)
	Pi_mat = matrix(0,M,n)
	for(m in 1:M){
		Pi_mat[m,] = generate_X_CPT_MC(nstep,log_lik_mat,Pi_)
	}
	return(Pi_mat)
}

generate_X_CPT_MC = function(nstep,log_lik_mat,Pi_){
	# log_lik_mat is the n-by-n matrix with entries log(q(X_i|Z_j))
	# this function runs the MC sampler, initialized with permutation Pi_
	n = length(Pi_)
	npair = floor(n/2)
	for(istep in 1:nstep){
	  perm = sample(n)
		inds_i = perm[1:npair]
		inds_j = perm[(npair+1):(2*npair)]
		# for each k=1,...,npair, deciding whether to swap Pi_[inds_i[k]] with Pi_[inds_j[k]]
    log_odds = (log_lik_mat[cbind(Pi_[inds_i],inds_j)] + log_lik_mat[cbind(Pi_[inds_j],inds_i)]
          - log_lik_mat[cbind(Pi_[inds_i],inds_i)] - log_lik_mat[cbind(Pi_[inds_j],inds_j)])
    swaps = rbinom(npair,1,1/(1+exp(-pmax(-500,log_odds))))
    Pi_[c(inds_i,inds_j)] = Pi_[c(inds_i,inds_j)] + swaps*(Pi_[c(inds_j,inds_i)]-Pi_[c(inds_i,inds_j)])
	}
	return(Pi_)
}

############# Function to generate X,Y,Z data ###################
generate_XYZ = function(type,param){
	a = rnorm(p)/p; b = rnorm(p)
	Z = matrix(rnorm(n*p),n,p)
	if(param==0 | type=='power'){
		X = Z%*%b + rnorm(n)
		Y = Z%*%a + X*param + rnorm(n)
	}else{
		if(type=='quadratic'){
			X = Z%*%b + param * (Z%*%b)^2 + rnorm(n)
		}
		if(type=='cubic'){
			X = Z%*%b - param * (Z%*%b)^3 + rnorm(n)
		}
		if(type=='tanh'){
			X = tanh(param*(Z%*%b))/param + rnorm(n)
		}
		if(type=='t'){
			X = Z%*%b + rt(n,1/param)*sqrt(1-2*param)
		}
		if(type=='skewnormal'){
			X = Z%*%b + (skewnormal(n,param)-sqrt(2/pi)*param/
			  sqrt(1+param^2))/sqrt(1-2/pi*param^2/(1+param^2))
		}
		if(type=='heteroskedastic'){
			X = Z%*%b +  rnorm(dim(Z)[1]) * abs(Z%*%b)^theta / 
             (sum(b^2))^(theta/2) / sqrt((2^theta)*gamma(theta+0.5)/sqrt(pi))
		}
		Y = Z%*%a + rnorm(n)
	}
	X_CPT = generate_X_CPT_gaussian(nstep,M,X,Z%*%b,rep(1,n))
	X_CRT = (Z%*%b)%*%t(rep(1,M)) + matrix(rnorm(n*M),n,M)
	XYZ = list()
	XYZ$X = X; XYZ$Y = Y; XYZ$Z = Z; XYZ$X_CPT = X_CPT; XYZ$X_CRT = X_CRT
	return(XYZ)
}

# rejection sampling from skew-normal distribution ( density 2*phi(t)*Phi(theta*t) )
skewnormal = function(n,theta){
	samples = NULL
	while(length(samples)<n){
		samples_ = rnorm(n)
		samples = c(samples,samples_[which(runif(n)<=pnorm(theta*samples_))])
	}
	return(samples[1:n])
}

############# Simulations under the alternative (testing power) #################


alpha = 0.05; nrep = 1000; n = 50; p = 20; M = 1000; nstep = 50

cs = seq(0, 0.7, 0.1)
cl <- length(cs)
true_cor <- matrix(nrow = cl, ncol = nrep)
cpt_cor <- array(dim = c(nrep, cl, M))
pb <- txtProgressBar(max = length(cs))
for(ic in 1:length(cs)){
  setTxtProgressBar(pb, ic)
	c = cs[ic]
	cors <-
	  foreach(irep = 1:nrep, .inorder = FALSE) %dopar% {
	    XYZ = generate_XYZ('power',c)
	    X = XYZ$X; Y = XYZ$Y; Z = XYZ$Z; X_CPT = XYZ$X_CPT
	    t_c <- c(abs(cor(X,Y)))
	    cpt_c<- c(abs(cor(X_CPT,Y)))
	    list(t_c, cpt_c)
	}
	true_cor[ic, ] <- sapply(cors, function(x) x[[1]])
	cpt_cor[, ic, ] <- do.call(rbind, lapply(cors, function(x) x[[2]]))
}
close(pb)
stopCluster(cluster)

# Save results for later use...
save(
  list = c("true_cor", "cpt_cor", "alpha", "nrep", "n", "p", "cs", "nstep", "M"),
  file = "Simulation1.rda"
)

# compute p-values and e-values
pvals <- numeric(cl)
for (i in seq_len(cl)) {
  p <- numeric(nrep)
  for (j in seq_len(nrep)) {
    p[j] <- (1 + sum(cpt_cor[j, i, ] >= true_cor[i, j])) / (M + 1)
  }
  pvals[i] <- mean(p <= alpha)
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

# Approach 1: Naive e-value T(X,Y,Z) / E_{X ~ Q(.|Z)}[T(X,Y,Z)]
evals <- get_evalue(true_cor, cpt_cor, alpha)

a1_rr <- tibble(c = cs, CPT = pvals, `E-value (1/E)` = evals[, 2]) %>%
  gather(key = "method", value = "rr", -c) %>%
  ggplot() +
  geom_point(aes(x = c, y = rr, shape = method, group = method)) +
  geom_line(aes(x = c, y = rr, group = method, linetype = method)) +
  geom_hline(yintercept = alpha, lty = 3) +
  theme(legend.position = c(0.2, 0.8), legend.title = element_blank()) +
  labs(
    x = expression("Parameter"~ italic("c")),
    y = "Rejection rate"
  ) +
  ggtitle("(a) Rejection rates of CPT and naive e-value")

a1_mean <- tibble(c = cs, e = evals[, 1]) %>%
  ggplot() +
  geom_point(aes(x = c, y = e)) +
  geom_line(aes(x = c, y = e)) +
  labs(
    x = expression("Parameter"~ italic("c")),
    y = "E-value (mean)"
  ) +
  ggtitle("(b) Mean of naive e-value")

a1_log_mean <- tibble(c = cs, e = evals[, 3]) %>%
  ggplot() +
  geom_point(aes(x = c, y = e)) +
  geom_line(aes(x = c, y = e)) +
  labs(
    x = expression("Parameter"~ italic("c")),
    y = "E-value (mean log)"
  ) +
  ggtitle("(c) Mean of log of naive e-value")

a1_plot <- ggarrange(
  a1_rr,
  ggarrange(a1_mean, a1_log_mean, ncol = 2),
  nrow = 2
)

pdf(width = 8, height = 6, file = "cpt_approach_1.pdf")
a1_plot
dev.off()

# Illustrate distribution of correlation
k <- 250
kc <- 6
tc <- true_cor[kc, k]
sc <- cpt_cor[k, kc, ]
title <- paste0(
  "Estimated density of cor(X,Y) in one simulation (c = ",
  cs[kc],
  ")"
)

a1_problem_plot <- ggplot() +
  geom_density(data = tibble(x = sc), aes(x = x)) +
  geom_vline(
    data = tibble(x = c(tc, mean(sc)), type = c("observed", "mean")),
    aes(xintercept = x, linetype = type)
  ) +
  scale_linetype_manual(values = c(5, 3)) +
  theme(legend.title = element_blank()) +
  labs(x = "cor(X,Y)", y = "Density") +
  ggtitle(title)

pdf(width = 8, height = 4, file = "cpt_approach_1_problem.pdf")
a1_problem_plot
dev.off()

# Approach 2: Take 1/f(cor(x,y)) / E[1/f(cor(X,y))] as e-value, where f is an
# estimate of the density of cor(X, y).
f <- function(obs, sample) {
  from <- min(sample, obs)
  to <- max(sample, obs)
  n <- 1024
  dens <- density(
    x = sample,
    from = from,
    to = to,
    n = n
  )
  dens_x <- dens$x
  dens_y <- dens$y
  obs_e <- 1/dens_y[which.min(abs(obs - dens_x))]^3
  n <- length(sample)
  e <- numeric(n)
  for (i in seq_len(n)) {
    e[i] <- 1/dens_y[which.min(abs(sample[i] - dens_x))]^3
  }
  obs_e / mean(e)
}
g <- function(obs, sample) {
  M <- length(sample)
  p <- (1 + sum(sample >= obs)) / (length(sample) + 1)
  1/sqrt(p) - 1
}
evals_2 <- get_evalue(true_cor, cpt_cor, alpha, f)
evals_pe <- get_evalue(true_cor, cpt_cor, alpha, g)

a2_rr <- tibble(
    c = cs,
    CPT = pvals,
    `E-value (1/E)` = evals_2[, 2],
    `1/sqrt(p[CPT]) - 1` = evals_pe[, 2]
  ) %>%
  gather(key = "method", value = "rr", -c) %>%
  mutate(
    method = factor(
      method,
      levels = c("E-value (1/E)", "1/sqrt(p[CPT]) - 1", "CPT"),
      labels = c("E-value (1/E)", "1/sqrt(p[CPT]) - 1", "CPT"),
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
  ggtitle("(a) Rejection rates of CPT and transformed e-value")

a2_mean <- tibble(
    c = cs,
    `E-value (1/E)` = evals_2[, 1],
    `1/sqrt(p[CPT]) - 1` = evals_pe[, 1]
  ) %>%
  gather(key = "method", value = "e", -c) %>%
  mutate(
    method = factor(
      method,
      levels = c("E-value (1/E)", "1/sqrt(p[CPT]) - 1"),
      labels = c("E-value (1/E)", "1/sqrt(p[CPT]) - 1"),
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

a2_log_mean <- tibble(
    c = cs,
    `E-value (1/E)` = evals_2[, 3],
    `1/sqrt(p[CPT]) - 1` = evals_pe[, 3]
  ) %>%
  gather(key = "method", value = "e", -c) %>%
  mutate(
    method = factor(
      method,
      levels = c("E-value (1/E)", "1/sqrt(p[CPT]) - 1"),
      labels = c("E-value (1/E)", "1/sqrt(p[CPT]) - 1"),
      ordered = TRUE
    )
  ) %>%
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

a2_plot <- ggarrange(
  a2_rr,
  ggarrange(a2_mean, a2_log_mean, ncol = 2),
  nrow = 2
)

pdf(width = 8, height = 6, file = "cpt_approach_2.pdf")
a2_plot
dev.off()
