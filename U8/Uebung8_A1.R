setwd("Serien")
library(ggplot2)
library(wordspace)

#library(factoextra)
#library(FactoMineR)

k_means_ <- function(data_df, k){
  data_mat <- as.matrix(data_df)
  start_cluster_idx <- sample(1:nrow(data_df),k)
  cluster_ctrs <- data_mat[as.vector(start_cluster_idx),]
  #print("cluster_ctrs")
  #print(cluster_ctrs)
  iter <- 0
  w_tildes <- double()
  repeat{
    iter <- iter+1
    clusters <- integer(nrow(data_mat))
    data_mat_with_ctrs <- rbind(data_mat,cluster_ctrs)
    #print(data_df_with_ctrs)
    dist_mat <- as.matrix(dist(data_mat_with_ctrs))
    #print("dist mat")
    #print(dist_mat)
    subset_m <- dist_mat[(nrow(dist_mat)-k+1):nrow(dist_mat),]
    subset_m <- subset_m[,1:(ncol(subset_m)-k)]
    #print("subset_m")
    #print(subset_m)
    clusters <- as.vector(mapply(which.min,split(subset_m, 
                                    rep(1:ncol(subset_m), each = nrow(subset_m)))))

    data_df_with_cl_ids <- data.frame(data_mat, clusters)
    
    
    old_cl_ctrs <- cluster_ctrs
    cluster_ctrs <- aggregate(. ~ clusters, data = data_df_with_cl_ids, FUN = mean)[,-1]
    colnames(cluster_ctrs) <- colnames(data_mat)
    # calc mod intra-cluster-var
    # this doesn't work
    cluster_sizes <- integer(k)
    diff_genes_means <- matrix(nrow=0, ncol=length(data_df))
    for(i in 1:k){
      cluster_sizes[i] <- sum(clusters==i)
      diff_data_cl_mean <- data_df_with_cl_ids[data_df_with_cl_ids$clusters==i,][,-1] - 
        as.vector(cluster_ctrs[i,])
      
      diff_data_cl_mean <- as.matrix(diff_data_cl_mean)
      colnames(diff_data_cl_mean) <- NULL
      colnames(diff_genes_means) <- NULL
      print(colnames(diff_genes_means))
      print(colnames(diff_data_cl_mean))
      diff_genes_means <- rbind(diff_genes_means, diff_data_cl_mean)
      
    }
    w_tilde2 <- 2*sum(rowNorms(diff_genes_means)^2)
    cat("w_tilde2 ", w_tilde2, "\n")
    
    #cluster_probe_vars <- aggregate(. ~ clusters, data = data_df_with_cl_ids, FUN = var)[,1:length(d)-1]
    #w_tilde_3 <- 2*sum(cluster_sizes* rowSums(cluster_probe_vars))
    #cat("w_tilde_3: ", w_tilde_3, "\n")
  
    # this works
    intra_var <- 0
    for (i in 1:k) {
      x_l <- as.matrix(data_mat_with_cl_ids[data_mat_with_cl_ids$clusters==i,])
      N_k <- nrow(x_l)
      if (nrow(x_l) > 1) {
        intra_var <- intra_var + N_k * sum(apply(x_l, 2, var, na.rm = TRUE))
      }
    }
    cat("w_tilde ", 2* intra_var,"\n")
    w_tildes[[length(w_tildes)+1]] <- 2* intra_var
  
    if(identical(old_cl_ctrs, cluster_ctrs)){
      plot(seq(1:length(w_tildes)), w_tildes, type='l')
      cat("iterations: ", iter-1, "\n")
      break
    }
  }
}

#data <- data.frame(matrix(runif(20,0,5),ncol=2))
#data <- data.frame(matrix(c(1,2,1,2,1.5,1.5,0,1),ncol=2))
#k_means_(data,2)

mda <- read.csv("mda.csv")

#rownames(mda) <- mda$X
mda <- mda[,-1]

k_means_(mda,2)

#res <- kmeans(mda, 2, nstart = 100)

#kmeans(mda, 3, nstart = 100)

#kmeans(mda, 4, nstart = 100)

#clstr_probes <- kmeans(t(mda), 5, nstart=100)

#fviz_cluster(clstr_probes, data = t(mda), geom = "point")