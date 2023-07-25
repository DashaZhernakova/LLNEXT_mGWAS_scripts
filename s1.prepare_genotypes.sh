module load BCFtools
module load PLINK

#
# Convert imputed genotypes from vcf to plink, extract only SNPs with imputation quality > 0.5
#
cd /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/

for chr in `seq 1 22`
do
 plink2 \
 --vcf /groups/umcg-llnext/tmp01/data/imputed/rawdata/LLnextv2_2920merged23122022.vcfs/${chr}.pbwt_reference_impute.vcf.gz \
 --extract-if-info "INFO > 0.5" \
 --make-bed --out ${chr}.info_filtered
 
 awk 'BEGIN {FS=OFS="\t"}; {if ($2 == ".") $2 = $1 ":" $4; print $0}' ${chr}.info_filtered.bim > ${chr}.info_filtered.bim.tmp
 mv ${chr}.info_filtered.bim.tmp ${chr}.info_filtered.bim
 
 awk 'BEGIN {FS=OFS="\t"}; {$2 = $1 ":" $4; print $0}' ${chr}.info_filtered.bim > ${chr}.info_filtered.bim.tmp
 mv ${chr}.info_filtered.bim.tmp ${chr}.info_filtered.bim
 
done


#
# Concatenate chromosomes, excluding multiallelic SNPs
#

rm merge_list.txt
for chr in `seq 1 22`
do
 echo "${chr}.info_filtered" >> merge_list.txt
done

module load PLINK/1.9-beta6-20190617
plink --merge-list merge_list.txt --make-bed --out all_chr.info_filtered

cat all_chr.filtered-merge.missnp > multiallelic_SNPs.txt

for chr in `seq 1 22`
do
 plink --bfile ${chr}.info_filtered \
 --exclude multiallelic_SNPs.txt \
 --make-bed --out ${chr}.info_filtered.nodup
done

rm merge_list.txt
for chr in `seq 1 22`
do
 echo "${chr}.info_filtered.nodup" >> merge_list.txt
done

plink --merge-list merge_list.txt --make-bed --out all_chr.info_filtered

#
# Subset genotypes to include only babies and filter
#

plink \
    --bfile all_chr.info_filtered \
    --keep qc/babies.txt \
    --maf 0.05 --hwe 1e-6 --geno 0.05 \
    --update-sex baby_sex.txt \
    --make-bed \
    --out all_chr.babies.with_rel.flt
    

#
# Make GRMs for fastGWA
#
cd /groups/umcg-llnext/tmp01/umcg-dzhernakova/mGWAS/genotypes
ln -s /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.babies.with_rel.flt.bed
ln -s /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.babies.with_rel.flt.bim
ln -s /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.babies.with_rel.flt.fam

cut -d" " -f1,2,5 all_chr.babies.with_rel.flt.fam > babies_gender.txt

ml PLINK

# make GRMs for GCTA

gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2//tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
grm=GRM/babies_flt_pruned


plink2 --bfile all_chr.babies.with_rel.flt --indep-pairwise 250 50 0.2 --out tmp_pruning
$gcta --bfile all_chr.babies.with_rel.flt  --extract tmp_pruning.prune.in  --make-grm --out $grm
$gcta --grm $grm --make-bK-sparse 0.05 --out ${grm}_sparse 
$gcta --bfile all_chr.babies.with_rel.flt --extract tmp_pruning.prune.in --make-grm-gz --out ${grm}.text

rm tmp.pruning*
