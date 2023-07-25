d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/mGWAS/
cd ${d}

#
# 1. Run GWAS using fastGWA
#
head -1 ${d}/data/M12_CLR.txt | sed "s:\t:\n:g" | tail -n+2 > ${d}/data/all_species.txt

cd ${d}/scripts/

# Submit jobs for CLR-transformed, non-zero CLR-transformed and binary tables
while read line
do
 bac="${line}"
 echo $bac
 sbatch -o logs/run_clr_${bac}.out -e logs/run_clr_${bac}.err -J $bac run_fastGWA.sh $bac clr
 sbatch -o logs/run_quant_${bac}.out -e logs/run_quant_${bac}.err -J $bac run_fastGWA.sh $bac quant
 sbatch -o logs/run_bin_${bac}.out -e logs/run_bin_${bac}.err -J $bac run_fastGWA.sh $bac bin
done < ../data/all_species.txt

#
# 2. Combine results
#

echo -e "Timepoint\tSpecies\tCHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP" > ${d}/results/clr_all_species_and_timepoints.5e-08.txt
echo -e "Timepoint\tSpecies\tCHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP" > ${d}/results/quant_all_species_and_timepoints.5e-08.txt
echo -e "Timepoint\tSpecies\tCHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tT\tSE_T\tP_noSPA\tBETA\tSE\tP\tCONVERGE" > ${d}/results/binary_all_species_and_timepoints.5e-08.txt
sort -m -k12,12g -T $TMPDIR -S 10G --buffer-size=1000  ${d}/results/*/clr/babies_*.fastGWA.5e-08.maf0.1.txt >> ${d}/results/clr_all_species_and_timepoints.5e-08.txt
sort -m -k12,12g -T $TMPDIR -S 10G --buffer-size=1000  ${d}/results/*/quant/babies_*.fastGWA.5e-08.txt >> ${d}/results/quant_all_species_and_timepoints.5e-08.txt
sort -m -k15,15g -T $TMPDIR -S 10G --buffer-size=1000  ${d}/results/*/binary/babies_*.fastGWA.5e-08.txt >> ${d}/results/binary_all_species_and_timepoints.5e-08.txt

echo -e "Timepoint\tSpecies\tCHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP" > ${d}/results/clr_all_species_and_timepoints.1e-05.txt
echo -e "Timepoint\tSpecies\tCHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP" > ${d}/results/quant_all_species_and_timepoints.1e-05.txt
echo -e "Timepoint\tSpecies\tCHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tT\tSE_T\tP_noSPA\tBETA\tSE\tP\tCONVERGE" > ${d}/results/binary_all_species_and_timepoints.1e-05.txt
cat  ${d}/results/*/clr/*.fastGWA.1e-05.maf0.1.txt >> ${d}/results/clr_all_species_and_timepoints.1e-05.txt
cat ${d}/results/*/quant/*.fastGWA.1e-05.txt >> ${d}/results/quant_all_species_and_timepoints.1e-05.txt
cat  ${d}/results/*/binary/*.fastGWA.1e-05.txt >> ${d}/results/binary_all_species_and_timepoints.1e-05.txt


#
# 3. Annotate
#
data_types=(clr quant binary)
for dt in ${data_types[@]}
do
f=${d}/results/${dt}_all_species_and_timepoints.5e-08.txt

# Add rs ids
python3 ${script_dir}/utils/add_rs_by_position.py \
    $f \
    /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/All_20180423.vcf.gz \
    3 \
    > ${f%txt}rsids.txt

# Clump
Rscript ${script_dir}/clump.R ${f%txt}rsids.txt

# Annotate
f=${d}/results/${dt}_all_species_and_timepoints.5e-08.rsids.txt.clumped_0.1.txt

col=6
genes=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/ensembl_b37_genes.bed
genes_prot=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/ensembl_b37_protcoding_genes.bed
gwas_catalog=${d}/gwas_catalog_v1.0.2-associations_e110_r2023-07-20.cut.combined.tsv

module load BEDTools


# Add overlapping genes
awk -v c=${col} '{FS=OFS="\t"}; {split($c, snp, ":"); if (NR != 1) {print snp[1], snp[2] - 1, snp[2], $c}}' $f > ${f}.tmp.bed
echo -e "SNP\tOverlapping_genes" > ${f}.tmp.genes
intersectBed -a ${f}.tmp.bed -b $genes -wa -wb | \
  python2.7 ${script_dir}/utils/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes

 python3 ${script_dir}/utils/add_columns_from_file_v2.py  \
  -i ${f} --header -f ${f}.tmp.genes -i_m 5 -f_m 0 -f_cols 1 \
  > ${f}.tmp.genes1.txt

# Add protein-coding genes within a 250kb window from the SNP
echo -e "SNP\tProtcoding_genes_within_250kb" > ${f}.tmp.genes2
windowBed -a ${f}.tmp.bed -b $genes_prot -w 250000 | \
python2.7 ${script_dir}/utils/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes2

python3 ${script_dir}/utils/add_columns_from_file_v2.py  \
  -i ${f}.tmp.genes1.txt -f ${f}.tmp.genes2 -i_m 5 -f_m 0 -f_cols 1 | \
  python3 ${script_dir}/utils/add_columns_from_file_v2.py \
    -i stdin \
    -f ${gwas_catalog} \
    -i_m 4 -f_m 2 -f_cols 3 --header \
    > ${f}.genes.gwas_catalog.txt


rm ${f}.tmp*
