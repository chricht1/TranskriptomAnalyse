#a.)
library(DESeq2)

data <- read.csv("mda-raw.csv", header = TRUE)
row.names(data) <- data$X # set first column of dataframe as rownames
data <- data[, -1] # remove first column

split_vec <- strsplit(colnames(data), "\\.") # split the treatment names and sample names 
                                             # from column names of original dataframe
treatment <- sapply(split_vec, `[`, 1)
sample <- sapply(split_vec, `[`, 2)

col_df <- data.frame(treatment, sample)
str(col_df)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = col_df,
                              design= ~ treatment) # to compare DMSO and BA later

print(dds)
dds <- DESeq(dds, useT = TRUE)
resultsNames(dds) # lists the coefficients
res <- results(dds, contrast=c("treatment","BA","DMSO"), alpha=0.05)#, lfcThreshold = 1, altHypothesis = "greaterAbs")
resOrdered <- res[order(res$pvalue),]
summary(resOrdered)

resOrdered$pvalue[1:3]
resOrdered$log2FoldChange[1:3]



#b.)
p_vals_signific <- resOrdered$pvalue < 0.05
n_p_vals_signific <- sum(p_vals_signific)

lg_fc_abs_gr_1 <- abs(resOrdered$log2FoldChange) > 1
n_lg_fc_abs_gr_1 <- sum(lg_fc_abs_gr_1)

cat("number of log fold change absolutes > 1: ", n_lg_fc_abs_gr_1, "\n")
cat("number of p values < 0.05: ", n_p_vals_signific, "\n")


library("VennDiagram")

n_inter <- sum(lg_fc_abs_gr_1 == p_vals_signific & lg_fc_abs_gr_1 == TRUE)

grid.newpage()                                        

diagram <- draw.pairwise.venn(area1=n_lg_fc_abs_gr_1,area2=n_p_vals_signific,cross.area=n_inter,
                              category=c("|log-fc|>1","p-value<0.05"), 
                              cex=1.5, cat.cex=1.5, cat.just = list(c(-0.8, 0.8), c(1.5, 1)))
require(gridExtra)
grid.arrange(gTree(children=diagram), top=textGrob(
  "Venn Diagram of DESeq2 p-Values and log-fc's",
  gp = gpar(fontsize = 15))
  )



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

log_p_val <- -log10(res$pvalue)
color <- ifelse(abs(res$log2FoldChange) <= 1, "gray", "black")
vulc_df <- data.frame(res$log2FoldChange, log_p_val, color=factor(color))

ggplot(data = vulc_df, aes(x = res.log2FoldChange, y = log_p_val)) +
  geom_point() + xlab("log-fc") + ylab("-log(p-value)") + 
  ggtitle("Vulcano Plot of DESeq2 results") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(aes(-9,-log10(0.05), label = "-log10(0.05)", vjust = -1)) +
  geom_point(aes(color = color)) + 
  scale_color_identity() + # use the color values directly
  scale_x_continuous(breaks = seq(-13, 6, by = 1)) +
  scale_y_continuous(breaks = sort(seq(0, 20, by = 1)))
#scale_color_manual(
#  values = c("gray" = "|log2FC| < 1", "black" = "|log2FC| > 1"),
#  name = "Legend")


#d.)
n <- length(res$pvalue)

sidak <- (res$pvalue < 1-(1-0.05)^(1/n))
de_sidak <- sum(sidak)

bonf <- (res$pvalue < 0.05/n)
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

de_holm_sidak <- fp_corr(resOrdered$pvalue, method='holm-sidak', alpha=0.05)

de_holm_bonf <-  fp_corr(resOrdered$pvalue, method='holm-bonf', alpha=0.05)

de_hochberg <- fp_corr(resOrdered$pvalue, method='hochberg', alpha=0.05)

de_benj_hb <- fp_corr(resOrdered$pvalue, method='benj-hochberg', alpha=0.05)

de_benj_yk <- fp_corr(resOrdered$pvalue, method='benj-yekutieli', alpha=0.05)

print(de_sidak)
print(de_bonf)
print(de_holm_sidak)
print(de_holm_bonf)
print(de_hochberg)
print(de_benj_hb)
print(de_benj_yk)


#e.)
data <- read.csv("mda.csv", header = TRUE)
data <- data[, -1]
dmso <- data[, 1:4]
ba <- data[, 5:8]
means_dmso <- rowMeans(dmso)
means_ba <- rowMeans(ba)
lg_fc_s2 <- means_ba - means_dmso
lg_fc_abs_s2 <- abs(lg_fc_s2)


plot(lg_fc_s2, res$log2FoldChange, xlab="log-fc Vorlesung", ylab="log-fc DESeq2", 
     main="Korrelation zwischen log-fc's mit DESeq2 und Vorlesung\nLinie zeigt perfekte Korrelation")
abline(0,1)

#(lg_fc_s2, col=2,pch=3)
#legend(-5,-7,c("DESeq2","Vorlesung"), col=c(1,2), pch=c(1,3))
