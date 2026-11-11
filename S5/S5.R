library("VennDiagram")

setwd("Serien")

#a.)
DE_data <- read.csv("res.csv")

DE_data['scores'] <- pmin(DE_data$X.log2FoldChange/2, -log10(DE_data$padj)/-log10(0.05))

DE_data <- DE_data[order(DE_data$scores, decreasing = TRUE), ]
DE_data$scores[1:3]

#b.)
pathways_test <- readRDS("test-anno.rds")

for (i in seq_along(pathways_test)) {
  pathway <- pathways_test[[i]]
  pathway_name <- names(pathways_test)[i]
  cat("Pathway:", pathway_name, " mit ", length(pathway), " Genen\n")
  scores_Pw <- DE_data$scores[DE_data$X %in% pathway]
  scores_notPw <- DE_data$scores[!(DE_data$X %in% pathway)]
  
  p_val <- ks.test(scores_Pw, scores_notPw, alternative="two.sided", exact=FALSE)$p.value
  cat("p = ", p_val, "\n")
}



#c.)
pathways <- readRDS("anno_s4.rds")
pw_names <- character(length(pathways))
p_vals <- double(length(pathways))

for (i in seq_along(pathways)) {
  pathway <- pathways[[i]]
  pathway_name <- names(pathways)[i]
  scores_Pw <- DE_data$scores[DE_data$X %in% pathway]
  scores_notPw <- DE_data$scores[!(DE_data$X %in% pathway)]
  if(length(scores_Pw)==0 | length(scores_notPw)==0){
    pw_names[i] <- pathway_name
    p_vals[i] <- NA_real_
    next
  }
  p_val <- ks.test(scores_Pw, scores_notPw, alternative="two.sided", exact=FALSE)$p.value
  pw_names[i] <- pathway_name
  p_vals[i] <- p_val
}

pws_enrichment <- data.frame(pathways <- pw_names, p_val <- p_vals)
colnames(pws_enrichment) <- c("pathway", "p_val")
rows_ks_test_NA <- pws_enrichment$p_val %in% NA_real_
pws_enrichment <- pws_enrichment[complete.cases(pws_enrichment),] # remove rows with "NA"

#(i)
pws_enrichment_ordered <- pws_enrichment[order(pws_enrichment$p_val),]
print(pws_enrichment_ordered[1:3,])



#(ii)
cat("num p_vals ks-test < 0.05: ", sum(pws_enrichment$p_val < 0.05), "\n")

calc_fisher_test <- function(gene_list, pathway, n_all_genes){
  
  vierfelder_tafel <- matrix(, nrow = 3, ncol = 3)
  vierfelder_tafel[3,3] <- n_all_genes
  vierfelder_tafel[1,1] <- sum(gene_list %in% pathway)
  vierfelder_tafel[2,1] <- sum(!(gene_list %in% pathway))
  vierfelder_tafel[1,2] <- sum(!(pathway %in% gene_list))
  vierfelder_tafel[2,2] <- n_all_genes - length(union(gene_list, pathway))
  vierfelder_tafel[1,3] <- vierfelder_tafel[1,1] + vierfelder_tafel[1,2]
  vierfelder_tafel[3,1] <- vierfelder_tafel[1,1] + vierfelder_tafel[2,1]
  vierfelder_tafel[2,3] <- vierfelder_tafel[3,3] - vierfelder_tafel[1,3]
  vierfelder_tafel[3,2] <- vierfelder_tafel[3,3] - vierfelder_tafel[3,1]
  
  odds_ratio <- (vierfelder_tafel[1,1]/vierfelder_tafel[1,2])/(vierfelder_tafel[2,1]/vierfelder_tafel[2,2])
  log_odds_ratio <- log2(odds_ratio)
  
  first_probability <- exp(lchoose(vierfelder_tafel[3,1], vierfelder_tafel[1,1]) +
                             lchoose(vierfelder_tafel[3,2], vierfelder_tafel[1,2]) -
                             lchoose(vierfelder_tafel[3,3], vierfelder_tafel[1,3])
  )
  #                    x = vierfelder_tafel[1,1],    # Gene in Pathway und d.e.
  #                    m = vierfelder_tafel[3,1],     #Gene im Pathway insgesamt
  #                    n = vierfelder_tafel[3,2],     # Gene nicht im Pathway
  #                    k = vierfelder_tafel[1,3]   # Gene d.e.
  
  #Vierfeldertafeln für die festen Randwerte sind nur von k abhängig!
  kmin <- max(0, vierfelder_tafel[1,3]+ vierfelder_tafel[3,1]-vierfelder_tafel[3,3])
  kmax <- min(vierfelder_tafel[1,3], vierfelder_tafel[3,1])
  if (kmin<=kmax){
    all_probs <- sapply(seq(kmin, kmax), function(x) {
      exp(
        lchoose(vierfelder_tafel[3,1], x) +
          lchoose(vierfelder_tafel[3,2], vierfelder_tafel[1,3] - x) -
          lchoose(vierfelder_tafel[3,3], vierfelder_tafel[1,3])
      )
    })
    return(list(
      log_odds_ratio = log_odds_ratio,
      p_val = sum(all_probs[all_probs <= first_probability])
    ))
  }
  return(list(
    log_odds_ratio = NA,
    p_val = NA
  ))
}

results <- data.frame(
  pathway = rep(NA, length(pathways)),
  log_odds_ratio = rep(NA, length(pathways)),
  abs_log_odds_ratio = rep(NA, length(pathways)),
  p_val = rep(NA, length(pathways)),
  stringsAsFactors = FALSE
)

pathways <- readRDS("anno_s4.rds") 
gene_list <- read.csv("DEG-list.txt")[[1]]
n_all_genes <- length(union(gene_list, unique(unlist(pathways))))


for (i in seq_along(pathways)) {
  l <- calc_fisher_test(gene_list, pathways[[i]], n_all_genes)
  results$pathway[i] <- names(pathways)[i]
  results$log_odds_ratio[i] <- l$log_odds_ratio
  results$abs_log_odds_ratio[i] <- abs(l$log_odds_ratio)
  results$p_val[i] <- l$p_val
}

cat("num p_vals fisher-test < 0.05: ", sum(results$p_val < 0.05), "\n")

cat("num all p_vals: ", nrow(pws_enrichment), "\n")

n_inter <- sum(pws_enrichment$p_val < 0.05 & results$p_val < 0.05)

grid.newpage()                                        
diagram <- draw.pairwise.venn(area1=sum(pws_enrichment$p_val < 0.05),area2=sum(results$p_val < 0.05),cross.area=n_inter,
                              category=c("p ks-test < 0.05","p fisher-test<0.05"), 
                              cex=1.5, cat.cex=1.5, cat.just = list(c(-0.8, 0.8), c(1.5, 1)))



#(iii)
results_aligned <- results[rows_ks_test_NA==FALSE,] # remove rows where ks-test was "NA"
plot(-log10(pws_enrichment$p_val), -log10(results_aligned$p_val), xlab="-log10(p) kolmogorov-smirnov-test", 
     ylab="-log10(p) fisher-test")
#flag log für logarithmische achsen
abline(0,1)



#(iv)
fdr <- 0.05
n <- nrow(pws_enrichment)

F_val <- sum(1/c(1:n))
de_benjamini_yekutieli <- 0
for(k in n:1){
  if(pws_enrichment_ordered$p_val[k] <= fdr * k/(n*F_val)){
    de_benjamini_yekutieli_ks <- k
    break
  }
}
cat("num significant of corrected pvalues ks-test: ", de_benjamini_yekutieli_ks, "\n")

results_ordered <- results[order(results$p_val),]
n <- nrow(results)
de_benjamini_yekutieli <- 0
for(k in n:1){
  if(results_ordered$p_val[k] <= fdr * k/(n*F_val)){
    de_benjamini_yekutieli_fisher <- k
    break
  }
}

by_kc_list <- numeric(nrow(pws_enrichment))
by_kc_list[1:de_benjamini_yekutieli_ks] = 1
by_fisher_list <- numeric(nrow(results))
by_fisher_list[1:de_benjamini_yekutieli_fisher] = 1
by_fisher_list <- by_fisher_list[rows_ks_test_NA==FALSE]

n_inter = sum(by_kc_list == by_fisher_list & by_kc_list==TRUE)

cat("num significant of corrected pvalues fisher-test: ", de_benjamini_yekutieli_fisher, "\n")

grid.newpage()                                    
diagram <- draw.pairwise.venn(area1=de_benjamini_yekutieli_ks,
                              area2=de_benjamini_yekutieli_fisher,
                              cross.area=n_inter,
                              category=c("corr p ks < 0.05","corr p fisher < 0.05"), 
                              cex=1.5, cat.cex=1.5, cat.just = list(c(1, -8), c(-0.2, -15)))

