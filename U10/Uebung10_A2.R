setwd("Serien")

em_alg_univar_gauss <- function(data, K, n_start, pis=NULL, mus=NULL, sigmas=NULL){
  set.seed(1)
  
  log_likelih_all_runs <- list()
  tolerance <- 1e-6
  
  x <- data$V2
  N <- length(x)
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
  
  pis_ <- double(K)
  mus_ <- double(K)
  sigmas_ <- double(K)
  
  for(runs in 1:n_start){
    # initialize parameters
    pis_ <- pis
    mus_ <- mus
    sigmas_ <- sigmas
    
    rnd_classes <- sample(K, N, replace=TRUE)
    if(infer_pis){
      pis_ <- double(K)
      for(k in 1:K){
        pis_[k] <- sum(rnd_classes==k)/N
      }
    }
    if(infer_mus){
      mus_ <- double(K)
      for (k in 1:K){
        mus_[k] <- mean(data[rnd_classes==k,][,2])
      }
    }
    if(infer_sigmas){  
      sigmas_ <- double(K)
      for (k in 1:K){
        sigmas_[k] <- sd(data[rnd_classes==k,][,2]) 
      }
    }
    
    log_likelihoods <- c()
    log_likelihood_old <- -Inf
    
    repeat{
      # E-step:
      gammas_t <- matrix(nrow=N, ncol=K)
      for(k in 1:K){
        gammas_t[,k] <- pis_[k] * dnorm(x, mean = mus_[k], sd = sigmas_[k])
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
        if(infer_mus){
          mus_[k] <- sum(gammas_t[,k] * x)/sum(gammas_t[,k])
        }
        if(infer_sigmas){
          sigmas_[k] <- sqrt(sum(gammas_t[,k]*(x-mus_[k])^2)/sum(gammas_t[,k]))
        }
      }
      # termination criterion
      if (abs(log_likelihood - log_likelihood_old) < tolerance) {
        break
      }
      
      log_likelihood_old <- log_likelihood
    }
    log_likelih_all_runs[[length(log_likelih_all_runs)+1]] <- log_likelihoods
  }
  return(log_likelih_all_runs)
}


eval_univar_em_alg <- function(log_likelih_all_runs, case, N, K){
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
    num_free_param <- n_free_param + K
  }
  if(case>=4){
    num_free_param <- n_free_param + K-1
  }
  bic_score <- max(max_log_likelihoods) - n_free_param/2 * log(N)
  cat("bic-score case", case, "= ", bic_score, "\n")
  plot(1, type="n", xlim=c(1, max(iteration_counts)), ylim=range(as.double(unlist(log_likelih_all_runs))),
       xlab="Iteration", ylab="log-likelihood", main=paste("log-likelihoods in 100 Läufen\nFall: ", case))
  for (i in 1:n_start) {
    log_likelihoods_one_run <- as.double(log_likelih_all_runs[[i]])
    lines(seq_along(log_likelihoods_one_run),log_likelihoods_one_run, col=rgb(0,0,1,0.2))
  }
  options(digits = 16)
  plot(iteration_counts, as.numeric(max_log_likelihoods),
       xlab = "Anzahl Iterationen",
       ylab = "maximale log-Likelihood",
       main = paste("Anzahl Iterationen vs. max log-Likelihood\nFall: ", case),
       pch = 19, col = "darkgreen",
       yaxt = 'n')
  axis(2, at = c(min(max_log_likelihoods), max(max_log_likelihoods)),
       labels = c(round(min(max_log_likelihoods),8), round(max(max_log_likelihoods),8)))
  hist(iteration_counts, main=paste("Histogram der Iterationszahlen\nFall: ", case))
  hist(max_log_likelihoods, 
       main=paste("Histogram der maximalen log-Likelihoods\nFall: ", case),
       xaxt = 'n')
  axis(1, at = c(min(max_log_likelihoods), max(max_log_likelihoods)),
       labels = c(round(min(max_log_likelihoods),8), round(max(max_log_likelihoods),8)))
  options(digits = 7)
}

data <- read.table("gauss_2.txt")

N=nrow(data)
K=2
n_start=100

#ground truth parameters
pis <- double(K)
sigmas <- double(K)
mus <- c(-1,1)
sigmas <- c(1,1)
pis <- c(0.5,0.5)
#for(k in 1:K){
#  pis[k] <- sum(data$V1==k)/nrow(data)
  #mus[k] <- mean(data[data$V1==k,][,2])
#  sigmas[k] <- sd(data[data$V1==k,][,2])
#}

log_likelih_all_runs <- em_alg_univar_gauss(data, K, n_start, pis=pis, sigmas=sigmas)
eval_univar_em_alg(log_likelih_all_runs, 1, N, K)
log_likelih_all_runs <- em_alg_univar_gauss(data, K, n_start, pis=pis, mus=mus)
eval_univar_em_alg(log_likelih_all_runs, 2, N, K)
log_likelih_all_runs <- em_alg_univar_gauss(data, K, n_start, pis=pis)
eval_univar_em_alg(log_likelih_all_runs, 3, N, K)
log_likelih_all_runs <- em_alg_univar_gauss(data, K, n_start, mus=mus, sigmas=sigmas)
eval_univar_em_alg(log_likelih_all_runs, 4, N, K)
log_likelih_all_runs <- em_alg_univar_gauss(data, K, n_start, sigmas=sigmas)
eval_univar_em_alg(log_likelih_all_runs, 5, N, K)
log_likelih_all_runs <- em_alg_univar_gauss(data, K, n_start, mus=mus)
eval_univar_em_alg(log_likelih_all_runs, 6, N, K)
log_likelih_all_runs <- em_alg_univar_gauss(data, K, n_start)
eval_univar_em_alg(log_likelih_all_runs, 7, N, K)

