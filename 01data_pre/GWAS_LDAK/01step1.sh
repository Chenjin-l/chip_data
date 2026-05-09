#!/bin/bash
#PBS -N growth_gwas
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -l mem=40G
#PBS -t 3-5

# ======================
# 变量
# ======================
n=${PBS_ARRAYID}

GENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/z01data/01phased_imputation/03imputation/02v4-phased-beagle/merge/QC/imp_qc"
PHENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/02pheno/pheno_growth_ldak.rmOutlier.3sd.txt"

COVAR_QUANT="/work/lingcj/03chip_data/02anlyis/01data/03growth/02pheno/covar_quant.txt"
#COVAR_QUANT="/work/lingcj/03chip_data/02anlyis/01data/03growth/02pheno/noage_covar_quant.txt"
COVAR_FACTOR="/work/lingcj/03chip_data/02anlyis/01data/03growth/02pheno/covar_factor.txt"

snp="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/z01data/01phased_imputation/03imputation/02v4-phased-beagle/merge/QC/imp_qc_ld.prune.in"

WORK_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/01step1"

index="/work/lingcj/03chip_data/02anlyis/01data/03growth/02pheno/index_growth.txt"

THREADS=8
pheno_col=$((PBS_ARRAYID-2))

PheName=$(sed -n ${PBS_ARRAYID}p $index)

# ======================
# 工作目录
# ======================
cd $WORK_DIR
mkdir -p ${n}_${PheName}
cd ${n}_${PheName}
rm -rf *

# ======================
# 环境
# ======================
source ~/anaconda3/etc/profile.d/conda.sh  
conda activate ldak_env

# ======================
# Step 1
# ======================
echo "Step1 KVIK..."

ldak6 --kvik-step1 gwas_${PheName}_step1 \
    --bfile $GENO \
    --extract ${snp} \
    --pheno $PHENO \
    --covar $COVAR_QUANT \
    --factors $COVAR_FACTOR \
    --mpheno $pheno_col \
    --max-threads $THREADS

# ======================
# Step 2
# ======================
echo "Step2 KVIK..."

ldak6 --kvik-step2 gwas_${PheName}_step1 \
    --bfile $GENO \
    --pheno $PHENO \
    --covar $COVAR_QUANT \
    --factors $COVAR_FACTOR \
    --mpheno $pheno_col \
    --max-threads $THREADS

# ======================
# 排序p值
# ======================
sort -gk2 gwas_${PheName}_step1.step2.pvalues > sort_min_pvalue.txt
