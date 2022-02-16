# load packages
library(doParallel)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(tidyr)
library(ggpubr)
library(FOCI)

# parallelization for CPR
nCores <- detectCores() - 1 # system("nproc", intern=TRUE)
cluster <- makeCluster(nCores)
registerDoParallel(cluster)

############# Implementation of the conditional permutation test ##################

generate_X_CPT_gaussian = function(nstep,M,X0,mu,sig2){
  # Runs the conditional permutation test using the distribution X | Z=Z[i] ~ N(mu[i],sig2[i])
  
  # For X[i] | Z=Z[i,], we have: p(X[i]|Z) \propto exp(-1/2(x-mu)^2/sigma^2) \propto exp(-1/(2sigma^2)x^2+x*mu/\sigma^2)
  log_lik_mat = -(X0^2)%*%t(1/2/sig2) + X0%*%t(mu/sig2)
  # log_lik_mat[i,j] = density at X=X0[i] when Z=Z[j]
  Pi_mat = generate_X_CPT(nstep,M,log_lik_mat)
  X_mat = X0[Pi_mat]
  dim(X_mat) = c(M,length(X0))
  return(t(X_mat))
}

generate_X_CPT_logistic = function(nstep,M,X0,mu,sig2){
  # Runs the conditional permutation test using the distribution X | Z=Z[i] ~ Ber(1/(1+exp(-mu)))
  
  log_lik_mat = -X0%*%t(log(1+exp(-mu)))-(1-X0)%*%t(log(exp(mu)+1))
  print(exp(log_lik_mat))
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
generate_XYZ = function(type,param,a){
  #a = rnorm(p)/p; b = rnorm(p)
  b=rnorm(p)
  Z = matrix(rnorm(n*p),n,p)
  if(type=='gaussian'){
    X = Z%*%b + rnorm(n)
    Y = Z%*%a + X*param + rnorm(n)
    X_CPT = generate_X_CPT_gaussian(nstep,M,X,Z%*%b,rep(1,n))
  } else if(type=='logistic') {
    # Problem that Y here is a vector and above matrix?
    X = rbinom(n, 1, 1/(1+exp(-Z%*%b)))
    Y = rbinom(n, 1, 1/(1+exp(-Z%*%a-param*X)))
    X_CPT = generate_X_CPT_logistic(nstep,M,X,Z%*%b)
  } else {
    stop("Specify a correct data generation type!") 
  }
  # Below is just a fancy way of writing (X = Z%*%b + rnorm(n)) M times.
  # That is, (Z%*%b)%*%t(rep(1,M)) gives a matrix with Z%*%b in each column. 
  # Then we add a random uniform matrix for the randomisation.
  # X_CRT = (Z%*%b)%*%t(rep(1,M)) + matrix(rnorm(n*M),n,M) 
  
  XYZ = list()
  XYZ$X = X; XYZ$Y = Y; XYZ$Z = Z; XYZ$X_CPT = X_CPT#; XYZ$X_CRT = X_CRT
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
#alpha = 0.05; nrep = 1000; n = 50; p = 20; M = 1000; nstep = 50
alpha = 0.05; nrep = 10; n = 100; p = 20; M = 1000; nstep = 50
density_type = "logistic"

# This function can/should be optimized:
if (density_type == "gaussian") {
  dens = function(x, Y, Z, a, c) { 
    return(dmvnorm(Y[,1], mean=(Z%*%a+c*x)[,1],sigma=diag(dim(Y)[1])))
  } 
} else if (density_type == "logistic") {
  dens = function(x, Y, Z, a, c) {
    pr = 1/(1+exp(-Z%*%a-c*x))
    return(prod(pr^Y*(1-pr)^(1-Y)))
  } 
} else {
  stop("Specify a correct density type.")
}

#cs = seq(0, 0.7, 0.1)
cs = seq(0,1,0.1) 

cl <- length(cs)
true_cor <- matrix(nrow = cl, ncol = nrep)
cpt_cor <- array(dim = c(nrep, cl, M))
true_dens <- matrix(nrow = cl, ncol = nrep)
cpt_dens <- array(dim = c(nrep, cl, M))
#true_codec <- matrix(nrow = cl, ncol = nrep)
#cpt_codec <- array(dim = c(nrep, cl, M))

pb <- txtProgressBar(max = length(cs))
for(ic in 1:length(cs)){
  setTxtProgressBar(pb, ic)
  c = cs[ic]
  cordens <-
    foreach(irep = 1:nrep, .inorder = FALSE, .packages = c('mvtnorm', 'FOCI' )) %dopar% {
      a = rnorm(p)/p
      XYZ = generate_XYZ(density_type,c,a)
      X = XYZ$X; Y = XYZ$Y; Z = XYZ$Z; X_CPT = XYZ$X_CPT
      t_c <- c(abs(cor(X,Y)))
      cpt_c <- c(abs(cor(X_CPT,Y)))
    
      t_d <- dens(X, Y, Z, a, c)
      cpt_d <- apply(X_CPT,MARGIN=2,FUN=function(x) dens(x,Y,Z,a,c))
      
      #t_codec <- codec(Y,X,Z)
      #cpt_codec <- apply(X_CPT, MARGIN=2, FUN=function(x) codec(Y,x,Z))
      list(t_c, cpt_c, t_d, cpt_d) #t_codec, cpt_codec)
    }
  true_cor[ic, ] <- sapply(cordens, function(x) x[[1]])
  cpt_cor[, ic, ] <- do.call(rbind, lapply(cordens, function(x) x[[2]]))
  true_dens[ic, ] <- sapply(cordens, function(x) x[[3]])
  cpt_dens[, ic, ] <- do.call(rbind, lapply(cordens, function(x) x[[4]]))
  #true_codec[ic, ] <- sapply(cordens, function(x) x[[5]])
  #cpt_codec[, ic, ] <- do.call(rbind, lapply(cordens, function(x) x[[6]]))
}
close(pb)
stopCluster(cluster)

# Save results for later use...
save(
  list = c("true_cor", "cpt_cor", "true_dens", "cpt_dens","alpha", "nrep", "n", "p", "cs", "nstep", "M"),
  file = sprintf("Simulation_%s.rda",density_type)
)