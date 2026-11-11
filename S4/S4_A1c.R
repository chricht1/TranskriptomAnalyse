library("VennDiagram")
library(ggplot2)

#a.)
fisher_test <- function(transcripts_DE, transcripts_Pws){
  results <- data.frame(pathway = character(), log_odds = double(), abs_log_odds = double(), fisher_p_val = double())
  #results <- data.frame(log_odds = double(), abs_log_odds = double(), fisher_p_val = double())
  transcripts_all_pws <- vector()
  for(tr_Pw_name in names(transcripts_Pws)){

    transcripts_all_pws <- union(transcripts_all_pws, transcripts_Pws[tr_Pw_name][[1]])
  }
  
  transcripts_DE_ <- unlist(transcripts_DE$V1)
  total <- length(union(transcripts_DE_, transcripts_all_pws))
  
  for(tr_Pw_name in names(transcripts_Pws)){
    transcripts_Pw <- transcripts_Pws[tr_Pw_name][[1]]
    
    DE <- length(transcripts_DE_)
    Pw <- length(transcripts_Pw)
    notDE <- total - DE
    notPw <- total - Pw
    
    DE_Pw <- length(intersect(transcripts_DE_, transcripts_Pw))
    DE_notPw <- DE - DE_Pw
    notDE_Pw <- Pw - DE_Pw
    notDE_notPw <- notDE - notDE_Pw
    
    fourfold <- matrix(c(DE_Pw, DE_notPw,  DE, 
                         notDE_Pw, notDE_notPw, notDE, 
                         Pw, notPw, total), nrow=3, byrow=TRUE)
    rownames(fourfold) <- c("DE", "not DE","")
    colnames(fourfold) <- c("Pw", "not Pw", "")
    
    #print(fourfold)
    #cat("\n")
    
    log_odds <- log(DE_Pw * notDE_notPw / (notDE_Pw * DE_notPw), base=2)
    
    #cat("log odds ratio:", log_odds, "and abs: ", abs(log_odds), "\n\n")
    
    
    p_hyp <- exp((lchoose(Pw, DE_Pw) + lchoose(notPw, DE_notPw)) - lchoose(total,DE))
    #cat("dhyper: ", p_hyp, "\n")
    
    p_hyp_arr <- as.list(numeric(DE+1))
    for(k in 0:DE){
      DE_Pw_new <- DE-k
      DE_notPw_new <- DE - DE_Pw_new
      p_hyp_arr[k] <- exp((lchoose(Pw, DE_Pw_new) + lchoose(notPw, DE_notPw_new)) - lchoose(total,DE))
    }
    
    p_hyp_arr <- unlist(p_hyp_arr)
    #cat("asserting sum probs of all generated outcomes == 1:", round(sum(p_hyp_arr),0) == 1, "\n")
    p_val <- sum(p_hyp_arr[p_hyp_arr <= p_hyp])
    #cat("p value of 2-sided fisher test: ", p_val, "\n")
    
    
    #test <- fisher.test(matrix(c(DE_Pw, DE_notPw, 
    #         notDE_Pw, notDE_notPw), nrow=2, byrow=TRUE))
    #print(test)
    results <- rbind(results, list(tr_Pw_name, log_odds, abs(log_odds), p_val))
    #results <- rbind(results, list(log_odds, abs(log_odds), p_val))
  }
  colnames(results) <- c("pathway", "log_odds", "abs_log_odds", "fisher_p_val")
  return(results)
}

#b.)
setwd("Serien")
transcripts_DE <- read.table("DEG-list.txt")
transcripts_Pws <- readRDS("anno.rds")

results <- fisher_test(transcripts_DE, transcripts_Pws)

#i
results <- results[order(results$fisher_p_val),]
print(results[1:3,])


#ii
abs_log_odds_gr_1 <- results$abs_log_odds > 1
n_abs_log_odds_gr_1 <- sum(abs_log_odds_gr_1)
cat("num of abs(log(odds)) > 1: ", n_abs_log_odds_gr_1, "\n")

p_vals_signific <- results$fisher_p_val < 0.05
n_p_vals_signific <- sum(p_vals_signific) 
cat("num of p-values < 0.05: ", n_p_vals_signific, "\n")


n_inter <- sum(abs_log_odds_gr_1 == p_vals_signific & abs_log_odds_gr_1 == TRUE)
print(n_inter)

grid.newpage()                                        
diagram <- draw.pairwise.venn(area1=n_abs_log_odds_gr_1, area2=n_p_vals_signific, cross.area=n_inter,
                              category=c("|log(odds)| > 1","p-value<0.05"), 
                              cex=2, cat.cex=1.7, cat.just = list(c(-0.8, 0.8), c(1.175, 1.9)))
require(gridExtra)
grid.arrange(gTree(children=diagram), top=textGrob(
  "Venn Diagram of fisher-test sign.-p-Values and |log2(odds)| > 1",
  gp = gpar(fontsize = 15))
)


#iii
results$group <- ifelse(
  (results$abs_log_odds <= 1 | results$fisher_p_val >= 0.05),
  "Not Significant",
  "Significant")

ggplot(data=results, aes(x=log_odds, y=-log10(fisher_p_val), color = group)) +
  geom_point() +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "lightgreen")) +
  xlab("log2(odds)") +
  ylab("-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(aes(4.5,-log10(0.05), label = "-log10(0.05)", vjust = -1)) +
  scale_x_continuous(breaks = seq(-4, 7, by = 1)) +
  scale_y_continuous(breaks = sort(seq(0, 8, by = 1))) +
  ggtitle("Vulcano Plot of fisher-tests and log-odds-ratios between \nDE-transcripts and different pathways")
  #geom_vline(xintercept= -1, linetype = "dashed") +
  #geom_vline(xintercept= 1, linetype = "dashed") +


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

de_benj_yk <- fp_corr(results$fisher_p_val, method='benj-yekutieli', alpha=0.05)

cat("num of sign. p-values after benjamini-yekutieli correction: ", de_benj_yk, "\n")