args <- commandArgs(trailingOnly = TRUE)

library(tidyr)
library(dplyr)
library(ieugwasr)

infile <- args[1]

d <- read.delim(infile, header = T, sep  ="\t", as.is = T, check.names = F)
d <- tibble(d %>% 
  rename(P = pval) %>% 
  unite(col = "id", Timepoint:Species, sep = ":", remove = F ))


res <-  ld_clump(d, plink_bin = genetics.binaRies::get_plink_binary(), bfile = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/1000G_snps_nodup.EUR", clump_r2 = 0.1, clump_kb = 250000)
write.table(res, file = paste0(infile, ".clumped_0.1.txt"), sep = "\t", quote = F, row.names = F)



