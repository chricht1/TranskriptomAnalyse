setwd("Serien")

mda <- read.csv("mda.csv")
mda <- mda[,-1]

###
### not used because to inefficient
###
create_dist_mat <- function(data, dist_measure){
  n <- length(colnames(data))
  dist_mat <- matrix(0,n,n)
  for(i in 1:n){
    for(j in i:n){
      if(i==j){
        dist_mat[i,j] <- 0
      }else{
        if(dist_measure=="euclidean"){
          dist_mat[i,j] <- dist(as.double(rbind(data[[i]],data[[j]])))[1]
        }else if(dist_measure=="inverse-sqrt-corr"){
          dist_mat[i,j] <- sqrt(1-cor(data[i],data[j]))
        }
      }
    }
  }
  return(as.dist(dist_mat))
}

#d_euc <- create_dist_mat(mda, "euclidean") 
d_euc <- dist(mda, method = "euclidean")

methods <- c("single", "complete","average")
for(method in methods){
  hc <- hclust(d_euc, method = method)
  plot(hc, main=paste0("dendogram mda euclidean distance ", method, "-linkage"))
}

#d_corr <- create_dist_mat(mda, "inverse-sqrt-corr")
cor_mat <- cor(mda)
d_corr <- as.dist(sqrt(1 - cor_mat))
for(method in methods){
  hc <- hclust(d_corr, method = method)
  plot(hc, main=paste0("dendogram mda inverse-sqrt-corr distance ", method, "-linkage"))
}