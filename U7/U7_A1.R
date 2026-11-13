setwd("Serien")

#a.)

read_dist_mat <- function(filename){
  dist_mat <- read.table(filename, header=TRUE)
  dist_mat_df <- data.frame(dist_mat, row.names=dist_mat$node)
  dist_mat_df <- dist_mat_df[,-1] # remove fist row
  colnames(dist_mat_df) <- dist_mat$node
  dist_mat <- as.dist(dist_mat_df)
  return (dist_mat)
}

d1_S6_A2_a <- read_dist_mat("d1.txt")
d2_S6_A2_b <- read_dist_mat("d2.txt")

hc <- hclust(d1_S6_A2_a, method = "average")
plot(hc)

hc <- hclust(d2_S6_A2_b, method = "average")
plot(hc)



#b.)

methods <- c("single", "complete")
for(method in methods){
  hc <- hclust(d1_S6_A2_a, method = method)
  plot(hc)
}

for(method in methods){
  hc <- hclust(d2_S6_A2_b, method = method)
  plot(hc)
}