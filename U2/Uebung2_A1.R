#install.packages("VennDiagram")
#install.packages("gridExtra")
#install.packages("tidyverse")
#install.packages("ggrepel")

#a.)
data <- read.csv("mda.csv", header = TRUE)
data <- data[, -1] # remove first col
dmso <- data[, 1:4]
ba <- data[, 5:8]
means_dmso <- rowMeans(dmso)
means_ba <- rowMeans(ba)
lg_fc <- means_ba - means_dmso
lg_fc_abs <- abs(lg_fc)


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
print(p_vals[1:5])

#p_vals_ <- mapply(function(x,y) {t.test(as.numeric(x), as.numeric(y), 
#                                     paired=TRUE, alternative="two.sided")$p.value}, 
#    split(dmso, 1:nrow(dmso)), 
#    split(ba, 1:nrow(ba)))
#print(p_vals_[1:5])

p_vals <- unlist(p_vals)
p_vals_sorted <- sort(p_vals)

print(p_vals_sorted[1:3])


#b.)
library("VennDiagram")

lg_fc_abs <- unlist(lg_fc_abs)
lg_fc_abs_gr_1 <- lg_fc_abs > 1
p_vals_signific <- p_vals < 0.05
n_lg_fc_abs_gr_1 <- sum(lg_fc_abs_gr_1)
n_p_vals_signific <- sum(p_vals_signific)
cat("number of log fold change absolutes greater than 1: ", n_lg_fc_abs_gr_1, "\n")
cat("number of p values less than 0.05: ", n_p_vals_signific, "\n")


n_inter <- sum(lg_fc_abs_gr_1 == p_vals_signific & lg_fc_abs_gr_1 == TRUE)


grid.newpage()                                        

diagram <- draw.pairwise.venn(area1=n_lg_fc_abs_gr_1,area2=n_p_vals_signific,cross.area=n_inter,
                   category=c("|log-fc|>1","p-value<0.05"), 
                   cex=1.5, cat.cex=1.5, cat.just = list(c(-0.8, 0.8), c(1.5, 1)))
require(gridExtra)
grid.arrange(gTree(children=diagram), top=textGrob(
  "Venn Diagram of p-Value and log-fc of transcripts",
  gp = gpar(fontsize = 15)
))


#c.)
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.7), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.7), color = 'black'),
              plot.title = element_text(hjust = 0.5, size=15),
              panel.background = element_rect(fill = NA),
              panel.grid.major = element_line(colour = "grey50", size=0.2)
            ))

log_p_val <- -log10(p_vals)
color <- ifelse(lg_fc_abs <= 1, "gray", "black")
vulc_df <- data.frame(lg_fc, log_p_val, color=factor(color))
print(str(vulc_df))

ggplot(data = vulc_df, aes(x = lg_fc, y = log_p_val)) +
  geom_point() + xlab("log-fc") + ylab("-log(p-value)") + 
  ggtitle("Vulcano Plot BA vs DMSO") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(aes(-4,-log10(0.05), label = "-log10(0.05)", vjust = -1)) +
  geom_point(aes(color = color)) + 
  scale_color_identity() + # use the color values directly
  scale_x_continuous(breaks = seq(-13, 13, by = 1)) +
  scale_y_continuous(breaks = sort(seq(-1, 16, by = 1)))
  #scale_color_manual(
  #  values = c("gray" = "|log2FC| < 1", "black" = "|log2FC| > 1"),
  #  name = "Legend")


#d.)
n <- length(p_vals)

sidak <- (p_vals < 1-(1-0.05)^(1/n))
de_sidak <- sum(sidak)

bonf <- (p_vals < 0.05/n)
de_bonf <- sum(bonf)


fp_corr <- function(p_vals_sorted, method, alpha){
  n <- length(p_vals_sorted)
  start <- 1
  end <- n
  if (method == 'hochberg' | method == 'benj-hochberg' | method== 'benj-yekutieli'){
    start <- n
    end <- 1
  }
  
  for (k in start:end) {
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
    
    if ((method == 'hochberg' | method == 'benj-hochberg' | method == 'benj-yekutieli') & 
        p_vals_sorted[k] <= threshold){
      return(k)
    }
    else if ((method == 'holm-sidak' | method == 'holm-bonf') & 
        p_vals_sorted[k] > threshold) {
      return(k-1)
    }
  }
}

de_holm_sidak <- fp_corr(p_vals_sorted, method='holm-sidak', alpha=0.05)

de_holm_bonf <-  fp_corr(p_vals_sorted, method='holm-bonf', alpha=0.05)

de_hochberg <- fp_corr(p_vals_sorted, method='hochberg', alpha=0.05)

de_benj_hb <- fp_corr(p_vals_sorted, method='benj-hochberg', alpha=0.05)

de_benj_yk <- fp_corr(p_vals_sorted, method='benj-yekutieli', alpha=0.05)

print(de_sidak)
print(de_bonf)
print(de_holm_sidak)
print(de_holm_bonf)
print(de_hochberg)
print(de_benj_hb)
print(de_benj_yk)
