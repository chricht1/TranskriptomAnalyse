setwd("Serien")

em_alg_multivar_gauss <- function(data, K, n_start, pis=NULL, mus=NULL, sigmas=NULL){
  set.seed(1)
  
  log_likelih_all_runs <- list()
  clusterings <- list()
  tolerance <- 1e-3
  
  N <- nrow(data)
  L <- ncol(data)
  
  infer_pis <- FALSE
  infer_mus <- FALSE
  infer_sigmas <- FALSE
  
  if(is.null(pis)){
    infer_pis <- TRUE
  }
  if(is.null(mus)){
    infer_mus <- TRUE
  }
  if(is.null(sigmas)){
    infer_sigmas <- TRUE
  }
  
  for(runs in 1:n_start){
    if((runs%%10)==0){
      cat("K=", K, " run:", runs, "\n")
    }
    # initialize parameters
    pis_ <- pis
    mus_ <- mus
    sigmas_ <- sigmas
    
    rnd_classes <- sample(K,N, replace=TRUE)
    if(infer_pis){
      pis_ <- double(K)
      for(k in 1:K){
        pis_[k] <- sum(rnd_classes==k)/n
      }
    }
    if(infer_mus){
      mus_ <- matrix(nrow=K,ncol=L)
      for(k in 1:K){
        for(l in 1:L){
          mus_[k,l] <- mean(data[rnd_classes==k,][,l])
        }
      }
    }
    if(infer_sigmas){  
      sigmas_ <- matrix(nrow=K,ncol=L)
      for (k in 1:K){
        for(l in 1:L){
          sigmas_[k,l] <- sd(data[rnd_classes==k,][,l])  
        }
      }
    }
    
    log_likelihoods <- c()
    log_likelihood_old <- -Inf
    
    repeat{
      # E-step:
      gammas_t <- matrix(nrow=N, ncol=K)
      for(k in 1:K){
        gammas_t[,k] <- pis_[k]
        for(l in 1:L){
          gammas_t[,k] <- gammas_t[,k] * dnorm(data[,l], mean = mus_[k,l], sd = sigmas_[k,l])
        }
      }
      likelihood_xn <- rowSums(gammas_t)
      gammas_t <- gammas_t/likelihood_xn
      
      log_likelihood <- sum(log(likelihood_xn))
      log_likelihoods <- c(log_likelihoods, log_likelihood)
      
      # M-step:
      for(k in 1:K){
        if(infer_pis){
          pis_[k] <- sum(gammas_t[,k])/N
        }
        for(l in 1:L){
          if(infer_mus){
            mus_[k,l] <- sum(gammas_t[,k] * data[,l])/sum(gammas_t[,k])
          }
          if(infer_sigmas){
            sigmas_[k,l] <- sqrt(sum(gammas_t[,k]*(data[,l]-mus_[k,l])^2)/sum(gammas_t[,k]))
          }
        }
      }
      # termination criterion
      if(abs(log_likelihood - log_likelihood_old) < tolerance){
        clusters <- max.col(gammas_t)
        clustering <- integer(K)
        for(k in 1:K){
          clustering[k] <- sum(clusters==k)
        }
        break
      }
      
      log_likelihood_old <- log_likelihood
    }
    log_likelih_all_runs[[length(log_likelih_all_runs)+1]] <- log_likelihoods
    clusterings[[length(clusterings)+1]] <- clustering
  }
  
  return(list(log_likelih_all_runs, clusterings))
}


eval_multivar_em_alg <- function(log_likelih_all_runs, case, N, K, L){
  n_start <- length(log_likelih_all_runs)
  iteration_counts <- integer(n_start)
  max_log_likelihoods <- double(n_start)
  for (i in 1:n_start) {
    log_likelihoods_one_run <- log_likelih_all_runs[[i]]
    iteration_counts[i] <- length(log_likelihoods_one_run)
    max_log_likelihoods[i] <- log_likelihoods_one_run[length(log_likelihoods_one_run)]
    if (!all(diff(log_likelihoods_one_run) >= 0)) {
      warning(paste("log likelihood ist in Fall: ", case, " nicht monoton steigend! Test:", i))
    }
  }
  n_free_param = 0
  if((case %% 2) == 1){
    n_free_param = K
  }
  if(case == 2|case == 3|case == 6|case == 7){
    num_free_param <- n_free_param + K * L
  }
  if(case>=4){
    num_free_param <- n_free_param + K * L
  }
  bic_score <- max(max_log_likelihoods) - n_free_param/2 * log(N)
  cat("bic-score K:", K, "= ", bic_score, "\n")
  plot(1, type="n", xlim=c(1, max(iteration_counts)), ylim=range(as.double(unlist(log_likelih_all_runs))),
       xlab="Iteration", ylab="log-likelihood", main=paste("log-likelihoods in 100 Läufen\nK= ", K))
  for (i in 1:n_start) {
    log_likelihoods_one_run <- as.double(log_likelih_all_runs[[i]])
    lines(seq_along(log_likelihoods_one_run),log_likelihoods_one_run, col=rgb(0,0,1,0.2))
  }
  
  plot(iteration_counts, as.numeric(max_log_likelihoods),
       xlab = "Anzahl Iterationen",
       ylab = "maximale log-Likelihood",
       main = paste("Anzahl Iterationen vs. max log-Likelihood\nK= ", K),
       pch = 19, col = "darkgreen")
  hist(iteration_counts, main=paste("Histogram der Iterationszahlen\nK= ", K))
  options(digits = 16)
  hist(max_log_likelihoods, 
       main=paste("Histogram der maximalen log-Likelihoods\nK= ", K),
       xaxt = 'n')
  axis(1, at = c(min(max_log_likelihoods), max(max_log_likelihoods)),
       labels = c(round(min(max_log_likelihoods),8), round(max(max_log_likelihoods),8)))
  options(digits = 7)
}

data <- read.csv("mda.csv")[,-1]

case=7
N=nrow(data)
L=ncol(data)

for (K in 2:5) {
  log_likelih_all_runs <- em_alg_multivar_gauss(data, K, n_start)
  eval_multivar_em_alg(log_likelih_all_runs, case, N, K, L)
}
