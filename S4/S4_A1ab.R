#a.)
fisher_test <- function(DEG_list, genes_Pws_list){
  #i
  DEG_vec <- as.vector(DEG_list)
  
  n_all_genes <- length(union(DEG_vec, unlist(genes_Pws_list)))
  n_DE <- length(DEG_list)
  
  for(tr_Pw_name in names(genes_Pws_list)){
    print(tr_Pw_name)
    transcripts_Pw <- genes_Pws_list[tr_Pw_name][[1]]
    
    n_Pw <- length(transcripts_Pw)
    n_notDE <- n_all_genes - n_DE
    n_notPw <- n_all_genes - n_Pw

    n_DE_Pw <- length(intersect(DEG_vec, transcripts_Pw))
    n_DE_notPw <- n_DE - n_DE_Pw
    n_notDE_Pw <- n_Pw - n_DE_Pw
    n_notDE_notPw <- n_notDE - n_notDE_Pw
    
    fourfold <- matrix(c(n_DE_Pw, n_DE_notPw,  n_DE, 
                         n_notDE_Pw, n_notDE_notPw, n_notDE, 
                         n_Pw, n_notPw, n_all_genes), nrow=3, byrow=TRUE)
    rownames(fourfold) <- c("DE", "not DE","")
    colnames(fourfold) <- c("Pw", "not Pw", "")
    
    print(fourfold)
    cat("\n")
    
    #ii
    odds <- n_DE_Pw * n_notDE_notPw / (n_notDE_Pw * n_DE_notPw)
    cat("odds ratio: ", odds, "\n")
    
    cat("log2 odds ratio:", log(odds, base=2), "abs log2 odds ratio", abs(log(odds, base=2)),"\n\n")
    
    
    #iii
    p_hyp <- exp((lchoose(n_Pw, n_DE_Pw) + lchoose(n_notPw, n_DE_notPw)) - lchoose(n_all_genes,n_DE))
    cat("dhyper: ", p_hyp, "\n")
    
    
    #iv
    p_hyp_arr <- as.list(numeric(n_DE+1))
    for(k in 0:n_DE){
      n_DE_Pw_new <- n_DE-k
      n_DE_notPw_new <- n_DE - n_DE_Pw_new
      p_hyp_arr[k] <- exp((lchoose(n_Pw, n_DE_Pw_new) + lchoose(n_notPw, n_DE_notPw_new)) - lchoose(n_all_genes,n_DE))
    }
    p_hyp_arr <- unlist(p_hyp_arr)
    cat("asserting sum probs of all generated outcomes == 1:", round(sum(p_hyp_arr),0) == 1, "\n")
    
    #v
    p_val <- sum(p_hyp_arr[p_hyp_arr <= p_hyp])
    cat("p value of 2-sided fisher test: ", p_val, "\n\n")
    
    
    #test <- fisher.test(matrix(c(DE_Pw, DE_notPw, 
    #         notDE_Pw, notDE_notPw), nrow=2, byrow=TRUE))
    #print(test)
  }
}

#b.)
setwd("Serien")
DEG_df <- read.table("DEG-list.txt")
genes_Pws_list <- readRDS("test-anno.rds")

DEG_list <- as.list(DEG_df[[1]])

fisher_test(DEG_list, genes_Pws_list)

