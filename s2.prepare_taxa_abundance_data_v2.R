library(tidyverse)

setwd("~/work/UMCG/data/NEXT/mbQTLs/data/")
taxa <- read.delim("LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt", header = T, row.names = 1, as.is = T, check.names = F)
metadata <- read.delim("LLNEXT_metadata_CLEAN_10_07_2023.txt", header = T,  as.is = T, check.names = F, row.names = 1)

# Filtering of infant taxa 

#1) Infants 
infant_metadata<-metadata[metadata$Type=="infant",]
infant_taxa<-taxa[row.names(taxa)%in% rownames(infant_metadata),] 
infant_taxa<-infant_taxa[match(row.names(infant_taxa),row.names(infant_metadata)),]
infant_NEXT_ID<-infant_metadata %>%
  select(NEXT_ID)
infant_taxa_all<-merge(infant_NEXT_ID,infant_taxa, by="row.names" )
row.names(infant_taxa_all)<-infant_taxa_all$Row.names
infant_taxa_all$Row.names=NULL
unique_counts <- sapply(infant_taxa_all, function(x) length(unique(infant_taxa_all$NEXT_ID[x >0.1]))) # Here I am counting the unique elements in the NEXT_ID column where the corresponding value in each column (i.e., x) is greater than the given cut-off. 
infant_taxa_all_filt <- infant_taxa_all[, unique_counts >= 0.1*length(unique(infant_taxa_all$NEXT_ID)) ] # Setting a 10% cut-off on prevalence 
infant_taxa_all_filt$NEXT_ID=NULL
infant_taxa_all_filt$UNCLASSIFIED=NULL

# Visualize number of samples
table(infant_metadata$Timepoint_categorical)
counts_per_timepoint <- dcast(infant_metadata[,c("NEXT_ID", "Timepoint_categorical")], formula = NEXT_ID~Timepoint_categorical, fun.aggregate = length)
duplicates <- counts_per_timepoint %>% filter_all(any_vars(. %in% c(2)))
counts_per_timepoint[counts_per_timepoint == 2] <- 1
upset(counts_per_timepoint)



# Only SGB level 
# Note: as opposed to previously, we are doing this on SGB level as there are some species that have mutiple SGB's like Faecalibacterium_prausnitzii, Sutterella_wadsworthensis
# For example: 
## cols_with_fae <- grep("Fae", colnames(infant_taxa_SGB), value = TRUE)
##cols_with_fae

infant_taxa_SGB_filt=infant_taxa_all_filt[,grep("t__",colnames(infant_taxa_all_filt))]

#Simplify names
colnames(infant_taxa_SGB_filt)=gsub(".*s__","",colnames(infant_taxa_SGB_filt))
names (infant_taxa_SGB_filt)
# To check this retains species of interest like those having two or more SGB's 
fae<-grep("Fae", colnames(infant_taxa_SGB_filt), value = TRUE) #[1] "Faecalibacterium_prausnitzii.t__SGB15316_group" "Faecalibacterium_prausnitzii.t__SGB15342" 
sut<-grep("Sut", colnames(infant_taxa_SGB_filt), value = TRUE) #"Sutterella_wadsworthensis.t__SGB9283" "Sutterella_wadsworthensis.t__SGB9286"


# CLR Transformation 
do_clr_externalWeighting = function(interest_matrix, core_matrix){
  interest_matrix = interest_matrix + min(core_matrix[core_matrix>0])/2
  core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  # estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

extract_timepoint_convert_ids <- function(taxa_clr, tp_ids, tp){
  
  # exclude one of the duplicates
  
  # subset the full sp abundance table to samples rom this timepoint
  tp_taxa <- taxa_clr[tp_ids,]
  
  #change the NG_ID to NEXT_ID for the GWAS and write it into a IID column
  all(row.names(tp_taxa) == row.names(meta_tp[tp_ids,]))
  row.names(tp_taxa) <- meta_tp[tp_ids,"NEXT_ID"]
  tp_taxa <- tp_taxa %>% 
    rownames_to_column(var = "IID")
  
  # add a 0-filled FID column required by fastGWA
  tp_taxa <- cbind(a=0, tp_taxa)
  colnames(tp_taxa)[1] <- "#FID"
  return(tp_taxa)
}

remove_one_of_duplicates <- function(meta2, duplicates, tp){
  tp_dup <- duplicates[duplicates[,tp] > 1, "NEXT_ID"]
  ids_to_rm <- c()
  for (dup in tp_dup){
    ids_to_rm <- c(ids_to_rm, infant_metadata[infant_metadata$NEXT_ID %in% dup & infant_metadata$Timepoint_categorical == tp,"NG_ID"][2])
  }
  return(ids_to_rm)
}



# SGB level transformed matx 
infant_taxa_SGB_all_filt_CLR<-do_clr_externalWeighting (infant_taxa_SGB_filt, infant_taxa_SGB_filt)
infant_taxa_SGB_all_filt_CLR<-as.data.frame(infant_taxa_SGB_all_filt_CLR)

# Binary SGB level table
infant_taxa_SGB_filt_binary <- as.data.frame(ifelse(infant_taxa_SGB_filt > 0, 1, 0))
infant_taxa_SGB_filt_binary <- infant_taxa_SGB_filt_binary[row.names(infant_taxa_SGB_all_filt_CLR), colnames(infant_taxa_SGB_all_filt_CLR)] # use only species that pass filtering

# quantitative SGB level table, replace 0s prior to CLR with NA
infant_taxa_SGB_all_filt_CLR_quant <- infant_taxa_SGB_all_filt_CLR
infant_taxa_SGB_all_filt_CLR_quant[infant_taxa_SGB_filt_binary == 0] <- NA

hist(infant_taxa_SGB_all_filt_CLR[,"Faecalibacterium_prausnitzii.t__SGB15316_group"], breaks = 100)
hist(infant_taxa_SGB_all_filt_CLR_quant[,"Faecalibacterium_prausnitzii.t__SGB15316_group"], breaks = 100)
hist(infant_taxa_SGB_filt_binary[,"Faecalibacterium_prausnitzii.t__SGB15316_group"], breaks = 100)



# write per timepoint tables
for (tp in unique(infant_metadata$Timepoint_categorical)){
  meta_tp <- infant_metadata[infant_metadata$Timepoint_categorical == tp, ]
  duplicate_ids <- remove_one_of_duplicates(infant_metadata, duplicates, tp)
  tp_ids <- setdiff(row.names(meta_tp), duplicate_ids)
  
  tp_taxa <- extract_timepoint_convert_ids(infant_taxa_SGB_all_filt_CLR, tp_ids, tp)
  write.table(tp_taxa, file = paste0("per_timepoint_v2/", tp, "_CLR.txt"), sep = "\t", quote = F, row.names = F)
  write.table(tp_taxa, file = paste0("per_timepoint_v2/", tp, "_CLR.noheader.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  tp_taxa_bin <- extract_timepoint_convert_ids(infant_taxa_SGB_filt_binary, tp_ids, tp)
  write.table(tp_taxa_bin, file = paste0("per_timepoint_v2/", tp, "_binary.txt"), sep = "\t", quote = F, row.names = F)
  write.table(tp_taxa_bin, file = paste0("per_timepoint_v2/", tp, "_binary.noheader.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  tp_taxa_quant <- extract_timepoint_convert_ids(infant_taxa_SGB_all_filt_CLR_quant, tp_ids, tp)
  write.table(tp_taxa_quant, file = paste0("../data/per_timepoint_v2/", tp, "_quant_CLR.txt"), sep = "\t", quote = F, row.names = F)
  write.table(tp_taxa_quant, file = paste0("per_timepoint_v2/", tp, "_quant_CLR.noheader.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  # prepare per timpoint covariates file
  covariates <- meta_tp[tp_ids,c("NEXT_ID", "clean_reads_FQ_1", "dna_conc", "exact_age_months_at_collection")]
  colnames(covariates) <- c("IID", "num_reads", "dna_conc", "age")
  covariates$num_reads <- scale(covariates$num_reads)
  covariates <- cbind(a=0, covariates)
  colnames(covariates)[1] <- "#FID"
  write.table(covariates, file = paste0("per_timepoint_v2/", tp, "_covariates.noheader.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
  
}

baby_meta <- unique(infant_metadata[,c("NEXT_ID", "Type", "FAMILY")])
mother_meta <- unique(metadata[metadata$Type == "mother",c("NEXT_ID", "Type", "FAMILY")])
mother_meta <- mother_meta[-grep("_2", mother_meta$NEXT_ID),]
mother_for_baby <- left_join(baby_meta, mother_meta, by = c("FAMILY" = "FAMILY"))
fut2 <- read.delim("FUT2_secretor.txt", header = T, sep = "\t", as.is= T, check.names = F)
fut2_family <- na.omit(left_join(mother_for_baby, fut2, by = c("NEXT_ID.x" = "NEXT_ID")) %>% left_join(.,fut2, by = c("NEXT_ID.y" = "NEXT_ID")))
