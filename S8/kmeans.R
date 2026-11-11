setwd("Serien")
library(ggplot2)
library(wordspace)

k_means_ <- function(data_df, k, nstart){
  w_tilde_curves_tials <- list()
  for(trial in 1:nstart){
    cat("trial: ", trial, "\n")
    data_mat <- as.matrix(data_df)
    cluster_ctrs <- data_mat[as.vector(sample(1:nrow(data_mat),k)),]
    iter <- 0
    w_tildes <- double()
    
    while(TRUE){
      dists_clusters <- matrix(nrow=0,ncol=nrow(data_mat))
      for(i in 1:k){
        dists_cluster <- rowSums(t((t(data_mat) - as.double(cluster_ctrs[i,]))^2))
        dists_clusters <- rbind(dists_clusters,dists_cluster)
      }
      clusters <- apply(dists_clusters, 2, which.min)
      
      data_mat_with_cl_ids <- data.frame(data_mat, clusters)
      
      new_cluster_ctrs <- aggregate(. ~ clusters, data = data_mat_with_cl_ids, FUN = mean)[,-1]
      
      intra_var <- 0
      for (i in 1:k) {
        x_l <- as.matrix(data_mat_with_cl_ids[data_mat_with_cl_ids$clusters==i,])
        N_k <- nrow(x_l)
        if (nrow(x_l) > 1) {
          intra_var <- intra_var + N_k * sum(apply(x_l, 2, var, na.rm = TRUE))
        }
      }
      #cat("w_tilde ", 2* intra_var,"\n")
      w_tildes[[length(w_tildes)+1]] <- 2* intra_var
      
      if(identical(new_cluster_ctrs, cluster_ctrs)){
        #cat("iterations: ", length(w_tildes), "\n")
        # plot for simple data example
        #for(i in 1:nrow(cluster_ctrs)){
        #  data_mat_with_cl_ids[nrow(data_mat_with_cl_ids) + 1,] = append(as.vector(cluster_ctrs[i,]),i)
        #}
        #print(ggplot(data_mat_with_cl_ids, aes(x = X1, y = X2, color = as.factor(clusters))) +
        #        geom_point(size = 3) +
        #        labs(color = "Cluster") +
        #        theme_minimal())
        
        break
      }
      cluster_ctrs <- new_cluster_ctrs
    }
    w_tilde_curves_tials[[length(w_tilde_curves_tials)+1]] <- unlist(w_tildes)
  }
  return(w_tilde_curves_tials)
}


# test k-means on simple example
#data <- data.frame(matrix(runif(20,0,5),ncol=2))
#data <- data.frame(matrix(c(1,1.5,1,1.5,1.5,1.5,0.5,0.5),ncol=2))
#k_means_(data,2)

mda <- read.csv("mda.csv")
#rownames(mda) <- mda$X
mda <- mda[,-1]
tests = 100
all_fW_traces <- k_means_(mda,2,tests)
#generate_plots(all_fW_traces, tests)

min_fW <- double(tests)
iteration_counts <- integer(tests)
for (i in 1:tests) {
  fW_trace <- all_fW_traces[[i]]
  min_fW[i] <- tail(fW_trace, 1)
  iteration_counts[i] <- length(fW_trace)
  if (!all(diff(fW_trace) <= 0)) {
    warning(paste("fW ist nicht monoton fallend! Test:", i))
  }
}

plot(1, type="n", xlim=c(1, max(iteration_counts)), ylim=range(as.double(unlist(all_fW_traces))),
     xlab="Iteration", ylab="fW", main="Intra-Cluster-Varianz (fW) in 100 Läufen")
for (i in 1:tests) {
  all_fW_traces <- as.double(unlist(all_fW_traces[i]))
  lines(seq_along(all_fW_traces),all_fW_traces, col=rgb(0,0,1,0.2))
}

## Histogramme:
hist(as.numeric(unlist(min_fW)), breaks = 20,
     main = "Histogramm der minimalen modifizierten Intra-Cluster-Varianzen",
     xlab = "Minimale fW",
     col = "skyblue",
     border = "white",
     xlim = c(593632.29, 593632.35)
)

hist(iteration_counts, breaks=20, main="Histogramm der Iterationsanzahl",
     xlab="Anzahl Iterationen", col="salmon", border="white")

## Scatterplot:
plot(iteration_counts, as.numeric(unlist(min_fW)),
     xlab = "Anzahl Iterationen",
     ylab = "Minimale fW",
     main = "Scatterplot: Iterationen vs. minimale fW",
     pch = 19, col = "darkgreen")


# b)

results_all <- list()

for (K in c(3,4,5)) {
  message("Starte K = ", K)
  all_fW_traces <- k_means_(mda, K, tests)
  
  min_fW <- double(tests)
  iteration_counts <- integer(tests)
  for (i in 1:tests) {
    fW_trace <- all_fW_traces[[i]]
    min_fW[i] <- tail(fW_trace, 1)
    iteration_counts[i] <- length(fW_trace)
    if (!all(diff(fW_trace) <= 0)) {
      warning(paste("fW ist nicht monoton fallend! Test:", i))
    }
  }
  plot(1, type = "n", xlim = c(1, max(iteration_counts)), ylim = range(as.double(unlist(all_fW_traces))),
       xlab = "Iteration", ylab = "fW", main = paste("fW-Traces für K =", K))
  for (i in 1:tests) {
    lines(all_fW_traces[[i]], col = rgb(0, 0, 1, 0.2))
  }
  hist(as.numeric(unlist(min_fW)), breaks = 20,
       main = paste("Histogramm der min. fW (K =", K, ")"),
       xlab = "Minimale fW", col = "skyblue", border = "white")
  hist(iteration_counts, breaks = 20,
       main = paste("Histogramm Iterationsanzahl (K =", K, ")"),
       xlab = "Anzahl Iterationen", col = "salmon", border = "white")
  plot(iteration_counts, min_fW,
       xlab = "Anzahl Iterationen", ylab = "Minimale fW",
       main = paste("Scatterplot: Iterationen vs. fW (K =", K, ")"),
       pch = 19, col = "darkgreen")
}




# c)
K <- 5
message("Starte clustering der proben", K)
all_fW_traces <- k_means_(t(mda),K,tests)
#generate_plots(all_fW_traces, tests)

min_fW <- double(tests)
iteration_counts <- integer(tests)
for (i in 1:tests) {
  fW_trace <- all_fW_traces[[i]]
  min_fW[i] <- tail(fW_trace, 1)
  iteration_counts[i] <- length(fW_trace)
  if (!all(diff(fW_trace) <= 0)) {
    warning(paste("fW ist nicht monoton fallend! Test:", i))
  }
}


plot(1, type="n", xlim=c(1, max(iteration_counts)), ylim=range(as.double(unlist(all_fW_traces))),,
     xlab="Iteration", ylab="fW", main="Intra-Cluster-Varianz (fW) in 100 Läufen")
for (i in 1:tests) {
  lines(all_fW_traces[[i]], col=rgb(0,0,1,0.2))
}

hist(as.numeric(unlist(min_fW)), breaks = 20,
     main = "Histogramm der minimalen modifizierten Intra-Cluster-Varianzen (Proben)",
     xlab = "Minimale fW",
     col = "skyblue",
     border = "white")
hist(iteration_counts, breaks=20, main="Histogramm der Iterationsanzahl (Proben)",
     xlab="Anzahl Iterationen", col="salmon", border="white")
plot(iteration_counts, min_fW,
     xlab = "Anzahl Iterationen",
     ylab = "Minimale fW",
     main = "Scatterplot: Iterationen vs. minimale fW (Proben)",
     pch = 19, col = "darkgreen")