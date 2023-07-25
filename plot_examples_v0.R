library(patchwork)

geno <- as.data.frame(t(read.delim("snps.genotypes.txt", as.is = T, sep = "\t", check.names = F, row.names = 1)))
ids <- row.names(geno)
geno <- as.data.frame(sapply(geno, function(x) as.factor(as.character(x))))
geno$"IID" <- ids

pdf("quantitative_top_hits.pdf")
plots <- list()
cnt <- 1
results <- c("W2-Escherichia_coli-17:56518534", "M1-Streptococcus_parasanguinis-5:52285098", "M1-Streptococcus_parasanguinis-4:189900854", "M2-Bifidobacterium_bifidum-8:13938718", "M3-Parabacteroides_distasonis-7:77836451", "M6-Bifidobacterium_bifidum-10:124514229", "M12-Phocaeicola_vulgatus-20:41719099")
for (res in results){
  tp = unlist(strsplit(res,"-"))[1]
  sp = unlist(strsplit(res,"-"))[2]
  snp = unlist(strsplit(res,"-"))[3]
  
  meta_tp <- meta2[meta2$Timepoint_categorical == tp, ]
  duplicate_ids <- remove_one_of_duplicates(meta2, duplicates, tp)
  tp_ids <- setdiff(row.names(meta_tp), duplicate_ids)
  tp_taxa <- taxa_clr[tp_ids,]
  all(row.names(tp_taxa) == row.names(meta_tp[tp_ids,]))
  row.names(tp_taxa) <- meta_tp[tp_ids,"NEXT_ID"]
  tp_taxa <- tp_taxa %>% 
    rownames_to_column(var = "IID")
  
  merged <- inner_join(tp_taxa, geno, by = c("IID" = "IID"))[,c("IID", sp, snp)]
  colnames(merged) <- c("IID", "sp", "snp")
  
  #ggplot(merged, aes(x = sp, fill = snp)) + geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')
  table(merged[,3])
  
  print(ggplot(merged, aes(x = snp, color = snp, y = sp)) +geom_point() + geom_boxplot(alpha = 0.2) + theme_bw())
  cnt <- cnt + 1
}
dev.off()
