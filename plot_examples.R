library(patchwork)
setwd("/Users/Dasha/work/UMCG/data/NEXT/mbQTLs/plots/")

geno <- as.data.frame(t(read.delim("snps.genotypes.txt", as.is = T, sep = "\t", check.names = F, row.names = 1)))
ids <- row.names(geno)
geno <- as.data.frame(sapply(geno, function(x) as.factor(as.character(x))))
geno$"IID" <- ids

toplot <- read.delim("toplot.txt", as.is = T, sep = "\t", check.names = F)
toplot$SNP <- paste0(toplot$CHR, ":", toplot$POS)

plots_hist <- list()
plots <- list()
cnt <- 1
pdf("v2_top_hits.pdf")
#results <- c("W2-Escherichia_coli-17:56518534", "M1-Streptococcus_parasanguinis-5:52285098", "M1-Streptococcus_parasanguinis-4:189900854", "M2-Bifidobacterium_bifidum-8:13938718", "M3-Parabacteroides_distasonis-7:77836451", "M6-Bifidobacterium_bifidum-10:124514229", "M12-Phocaeicola_vulgatus-20:41719099")
for (i  in 1:nrow(toplot)){
  tp = toplot$Timepoint[i]
  sp = toplot$Species[i]
  snp = toplot$SNP[i]
  
  meta_tp <- infant_metadata[infant_metadata$Timepoint_categorical == tp, ]
  duplicate_ids <- remove_one_of_duplicates(infant_metadata, duplicates, tp)
  tp_ids <- setdiff(row.names(meta_tp), duplicate_ids)
  
    if (toplot$type[i] == 'clr'){
    tp_taxa <- extract_timepoint_convert_ids(infant_taxa_SGB_all_filt_CLR, tp_ids, tp)
  } else if (toplot$type[i] == 'quant'){
    tp_taxa <- extract_timepoint_convert_ids(infant_taxa_SGB_all_filt_CLR_quant, tp_ids, tp)
  }
  
     
  merged <- na.omit(inner_join(tp_taxa, geno, by = c("IID" = "IID"))[,c("IID", sp, snp)])
  colnames(merged) <- c("IID", "sp", "snp")
  
  colors = c("#1D91C0", "#EC7014", "#00A072")
  
  merged <- merged %>% 
    dplyr::count(snp) %>% 
    mutate(snp_with_n = paste0(snp,' (N = ', n, ")") ) %>%
    left_join(merged, ., by = "snp")
  p1 <- ggplot(merged, aes(x = sp, fill = snp_with_n)) + 
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') + 
    scale_fill_manual(values = colors, name = toplot$rsid[i]) + 
    xlab(sp) +
    theme_bw()

  p2 <- ggplot(merged, aes(x = snp, color = snp_with_n, y = sp)) +
          geom_jitter(alpha = 0.5, width = 0.2) + 
          geom_boxplot(alpha = 0.2) + 
          scale_color_manual(values = colors, name = toplot$rsid[i]) +
          ylab(sp) +
          theme_bw() +
          ggtitle(paste0("tp = ", tp, "; N = ",toplot$N, "; P = ", formatC(toplot$pval[i], digits = 3)))
  print(p1/p2)
}


dev.off()





