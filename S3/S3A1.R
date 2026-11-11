# load and preprocess data
data <- read.csv("mda.csv", header = TRUE)
data <- data[, -1] # remove first col
dmso <- data[, 1:4]
ba <- data[, 5:8]
means_dmso <- rowMeans(dmso)
means_ba <- rowMeans(ba)
lg_fc <- means_ba - means_dmso
lg_fc_abs <- abs(lg_fc)

# calc p-values
calc_t_test <- function(x,y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  L1 <- length(x)
  L2 <- length(y)
  v_x <- var(x)
  v_y <- var(y)
  pooled_var <- ((L1-1)*v_x+(L2-1)*v_y)/(L1+L2-2)
  t_val <- (mean(x)-mean(y))*sqrt(L1*L2/(L1+L2))/sqrt(pooled_var)
  if(t_val <= 0){
    p_val <- 2*pt(t_val, L1+L2-2)
  }
  else{
    p_val <- 2*pt(-t_val, L1+L2-2)
  }
  return(p_val)
}
p_vals <- mapply(calc_t_test, 
                 split(dmso, 1:nrow(dmso)), 
                 split(ba, 1:nrow(ba)))

p_vals <- unlist(p_vals)
p_vals_sorted <- sort(p_vals)

n <- length(p_vals)
print(n)
print("number of p-values less than 0.05")
print(sum(p_vals < 0.05))


# FWER and FDR corrections
sidak <- (p_vals < 1-(1-0.05)^(1/n))
de_sidak <- sum(sidak)

bonf <- (p_vals < 0.05/n)
de_bonf <- sum(bonf)



fp_corr <- function(p_vals_sorted, method, alpha){
  n <- length(p_vals_sorted)
  is_de_list <- vector(mode = "logical", length = n)
  flag <- FALSE
  if (method == 'hochberg' | method == 'benj-hochberg' | method== 'benj-yekutieli'){
    flag <- TRUE
  }
  for (k in 1:n) {
    if (method == 'holm-sidak'){
      threshold <- 1 - (1 - alpha)^(1 / (n + 1 - k))
    }else if (method == 'holm-bonf' | method == 'hochberg'){
      threshold <- alpha/(n+1-k)
    }else if (method == 'benj-hochberg'){
      threshold <- alpha*k/n
    }else if (method == 'benj-yekutieli'){
      f <- sum(1/(1:n))
      threshold <- alpha*k/(f*n)
    }
    
    if (flag == TRUE & 
        (method == 'hochberg' | method == 'benj-hochberg' | method == 'benj-yekutieli') & 
        p_vals_sorted[k] <= threshold){
      flag <- FALSE #declare all following Positions as differentially expressed 
    }
    else if (flag == FALSE & (method == 'holm-sidak' | method == 'holm-bonf') & p_vals_sorted[k] > threshold) {
      flag <- TRUE
    }
    
    if (flag) {
      is_de_list[k] <- FALSE
    } else {
      is_de_list[k] <- TRUE
    }
  }
  return(is_de_list)
}

holm_sidak <- fp_corr(p_vals_sorted, method='holm-sidak', alpha=0.05)
de_holm_sidak <- sum(holm_sidak)

holm_bonf <-  fp_corr(p_vals_sorted, method='holm-bonf', alpha=0.05)
de_holm_bonf <- sum(holm_bonf)

p_vals_sorted_rev <- rev(p_vals_sorted)
hochberg <- fp_corr(p_vals_sorted_rev, method='hochberg', alpha=0.05)
hochberg <- rev(hochberg)
de_hochberg <- sum(hochberg)

benj_hb <- fp_corr(p_vals_sorted_rev, method='benj-hochberg', alpha=0.05)
benj_hb <- rev(benj_hb)
de_benj_hb = sum(benj_hb)

benj_yk <- fp_corr(p_vals_sorted_rev, method='benj-yekutieli', alpha=0.05)
benj_yk <- rev(benj_yk)
de_benj_yk = sum(benj_yk)

print(de_sidak)
print(de_bonf)
print(de_holm_sidak)
print(de_holm_bonf)
print(de_hochberg)
print(de_benj_hb)
print(de_benj_yk)