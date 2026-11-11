setwd("Serien")

em_alg <- function(data, K, n_start){
  set.seed(1)
  
  final_mus <- list()
  log_likelihoods_all_runs <- list()
  
  x <- data$V2
  n <- length(x)
  
  # Initialize parameters
  pis <- double(K)
  sigmas <- double(K)
  mus <- double(K)
  for(k in 1:K){
    pis[k] <- sum(data$V1==k)/n
    sigmas[k] <- sd(data[data$V1==k,][,2])
  }
  
  for(runs in 1:n_start){
    log_likelihoods <- c()
    # initialize parameters
    for (k in 1:K){
      #mus[k] <- mean(sample(data[data$V1==k,][,2], n/2, replace = TRUE))
      mus[k] <- rnorm(1, mean = mean(data[data$V1==k,][,2]), sd=2)
    }
    
    log_likelihood_old <- -Inf
    tolerance <- 1e-6
    
    repeat{
      # E-step:
      gammas_t <- matrix(nrow=n, ncol=K)
      for(k in 1:K){
          gammas_t[,k] <- pis[k] * dnorm(x, mean = mus[k], sd = sigmas[k])
      }
      likelihood_x <- rowSums(gammas_t)
      gammas_t <- gammas_t/likelihood_x
      
      # M-step:
      mus_new <- double(K)
      for(k in 1:K){
        mus_new[k] <- sum(gammas_t[,k] * x)/sum(gammas_t[,k])
      }
      
      log_likelihood <- sum(log(likelihood_x))
      
      # termination criterion
      if (abs(log_likelihood - log_likelihood_old) < tolerance) {
        break
      }
      
      mus <- mus_new
      log_likelihood_old <- log_likelihood
      log_likelihoods <- c(log_likelihoods, log_likelihood)
    }
    final_mus[[length(final_mus)+1]] <- mus
    log_likelihoods_all_runs[[length(log_likelihoods_all_runs)+1]] <- log_likelihoods
  }
  return(list(final_mus,log_likelihoods_all_runs))
}

data <- read.table("gauss_2.txt")

K=2
n_start=100

res <- em_alg(data,K, n_start)

iteration_counts <- integer(n_start)
max_log_likelihoods <- double(n_start)
log_likelihoods_all_runs <- res[[2]]
for (i in 1:n_start) {
  log_likelihoods_one_run <- log_likelihoods_all_runs[[i]]
  iteration_counts[i] <- length(log_likelihoods_one_run)
  max_log_likelihoods[i] <- log_likelihoods_one_run[length(log_likelihoods_one_run)]
  if (!all(diff(log_likelihoods_one_run) >= 0)) {
    warning(paste("log likelihood ist nicht monoton steigend! Test:", i))
  }
}

plot(1, type="n", xlim=c(1, max(iteration_counts)), ylim=range(as.double(unlist(log_likelihoods_all_runs))),
     xlab="Iteration", ylab="log-likelihood", main="log-likelihoods in 100 Läufen")
for (i in 1:n_start) {
  log_likelihoods_one_run <- as.double(log_likelihoods_all_runs[[i]])
  lines(seq_along(log_likelihoods_one_run),log_likelihoods_one_run, col=rgb(0,0,1,0.2))
}

plot(iteration_counts, as.numeric(max_log_likelihoods),
     xlab = "Anzahl Iterationen",
     ylab = "maximale log-Likelihood",
     main = "Scatterplot: Iterationen vs. max log-Likelihood",
     pch = 19, col = "darkgreen")

hist(iteration_counts)

hist(max_log_likelihoods,
     xaxt = 'n')
axis(1, at = c(min(max_log_likelihoods), max(max_log_likelihoods)), 
     labels = c(round(min(max_log_likelihoods),4), round(max(max_log_likelihoods),4)))

final_mus <- matrix(unlist(res[[1]]),ncol=2,byrow=TRUE)
plot(final_mus[,1],final_mus[,2])

K=2
x <- data$V2
n <- length(x)
gammas_t <- matrix(nrow=n, ncol=K)
pis <- double(K)
sigmas <- double(K)
mus <- double(K)
for(k in 1:K){
  pis[k] <- sum(data$V1==k)/n
  sigmas[k] <- sd(data[data$V1==k,][,2])
  mus[k] <- mean(data[data$V1==k,][,2])
  gammas_t[,k] <- pis[k] * dnorm(x, mean = mus[k], sd = sigmas[k])
}
cat("true means k=1: ", mus[1], "|k=2: ",mus[2], "\n")

likelihood_x <- rowSums(gammas_t)
log_likelihood <- sum(log(likelihood_x))
cat("true log-likelihood: ", log_likelihood, "\n")