#!/bin/bash
#SBATCH --job-name=SV
#SBATCH --output=logs/run_GWAS.out
#SBATCH --error=logs/run_GWAS.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

bac=$1
type=$2
timepoints=(W2 M1 M2 M3 M6 M12)

echo "type=$type, bac=$bac"

d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/mGWAS/
gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2//tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1

for tp in ${timepoints[@]}
do 
    mkdir -p ${d}/results/${tp}/
    geno_file=${d}/genotypes/all_chr.babies.with_rel.flt
    gender_file=${d}/genotypes/babies_gender.txt
    grm=${d}/genotypes/GRM/babies_flt_pruned_sparse
    
    covar_file=${d}/data/${tp}_covariates.noheader.txt

    if [ $type == "bin" ]
    then
        gcta_mode="--fastGWA-mlm-binary"
        pheno_file=${d}/data/${tp}_binary.txt
        res_file=${d}/results/${tp}/binary/babies_${tp}_binary_${bac}
        mkdir -p ${d}/results/${tp}/binary/
    elif [ $type == "quant" ]
    then
        gcta_mode="--fastGWA-mlm"
        pheno_file=${d}/data/${tp}_quant_CLR.txt
        res_file=${d}/results/${tp}/quant/babies_${tp}_quant_${bac}
        mkdir -p ${d}/results/${tp}/quant/
    elif [ $type == "clr" ]
    then
        gcta_mode="--fastGWA-mlm"
        pheno_file=${d}/data/${tp}_CLR.txt
        res_file=${d}/results/${tp}/clr/babies_${tp}_clr_${bac}
        mkdir -p ${d}/results/${tp}/clr/
    else
	    echo "Wrong input data type: $type"
	    exit 1
    fi
    echo "type=$type, gcta_mode=$gcta_mode"
    echo "pheno file: $pheno_file"
    
    #get the column number for the SV in the SV table
    col=`head -1 ${pheno_file} | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${bac} | cut -d ":" -f1`


    $gcta \
      --bfile $geno_file \
      --maf 0.05 \
      --grm-sparse ${grm} \
      ${gcta_mode} \
      --pheno ${pheno_file%txt}noheader.txt \
      --mpheno $col \
      --qcovar $covar_file \
      --covar $gender_file \
      --out $res_file
    
    if [ $gcta_mode == "--fastGWA-mlm" ]
    then
        awk -v sp=$bac -v t=$tp 'BEGIN {FS=OFS="\t"}; {if (NR > 1 && $10 < 5e-8 && $6 > 50) print t, sp, $0}' ${res_file}.fastGWA | sort -k12,12g > ${res_file}.fastGWA.5e-08.txt
        awk -v sp=$bac -v t=$tp 'BEGIN {FS=OFS="\t"}; {if (NR > 1 && $10 < 1e-5 && $6 > 50) print t, sp, $0}' ${res_file}.fastGWA > ${res_file}.fastGWA.1e-05.txt
    else
        awk -v sp=$bac -v t=$tp 'BEGIN {FS=OFS="\t"}; {if (NR > 1 && $13 < 5e-8 && $6 > 50) print t, sp, $0}' ${res_file}.fastGWA | sort -k13,13g > ${res_file}.fastGWA.5e-08.txt
        awk -v sp=$bac -v t=$tp 'BEGIN {FS=OFS="\t"}; {if (NR > 1 && $13 < 1e-5 && $6 > 50) print t, sp, $0}' ${res_file}.fastGWA > ${res_file}.fastGWA.1e-05.txt
    fi
    gzip -f ${res_file}.fastGWA 
    if [ $type == "clr" ]
    then
	    awk 'BEGIN {FS=OFS="\t"}; {if ($9 > 0.1) print}' ${res_file}.fastGWA.5e-08.txt > ${res_file}.fastGWA.5e-08.maf0.1.txt
	    awk 'BEGIN {FS=OFS="\t"}; {if ($9 > 0.1) print}' ${res_file}.fastGWA.1e-05.txt > ${res_file}.fastGWA.1e-05.maf0.1.txt
	    rm ${res_file}.fastGWA.5e-08.txt
	    rm ${res_file}.fastGWA.1e-05.txt
    fi
done