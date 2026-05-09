#!/bin/bash
#PBS -N ldak_gwas
#PBS -q batch
#PBS -l nodes=1:ppn=16
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -l mem=50G
#PBS -t 3-5

# =========================
# 1. 样本编号（job array）
# =========================
n=${PBS_ARRAYID}

# =========================
# 2. 基础数据（Genotype）
# =========================
GENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/imp_QC"

# SNP过滤列表（LD prune）
snp="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/imp_QC.prune.in"

# =========================
# 3. 表型数据
# =========================
PHENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/pheno_growth_ldak.rmOutlier.3sd.txt"

# =========================
# 4. 协变量（Covariates）
# =========================
COVAR_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar"

COVAR_QUANT="${COVAR_DIR}/covar_quant.txt"
COVAR_FACTOR="${COVAR_DIR}/covar_factor.txt"

# （可选）不含age版本
# COVAR_QUANT="${COVAR_DIR}/noage_covar_quant.txt"

# =========================
# 5. 索引文件 【变量名统一】
# =========================
INDEX_FILE="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/index_growth.txt"

# =========================
# 6. 工作目录
# =========================
WORK_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/02LDAK_gwas/01step1"

THREADS=32
pheno_col=$((PBS_ARRAYID-2))  # 3→1,4→2,5→3 正确

# ======================
# 【修复】这里必须用 INDEX_FILE
# ======================
PheName=$(sed -n "${PBS_ARRAYID}p" "$INDEX_FILE")

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
# Step 2 【修复】必须指定输出名
# ======================
echo "Step2 KVIK..."

ldak6 --kvik-step2 gwas_${PheName}_step2 \
    --bfile $GENO \
    --pheno $PHENO \
    --covar $COVAR_QUANT \
    --factors $COVAR_FACTOR \
    --mpheno $pheno_col \
    --max-threads $THREADS

# ======================
# 排序p值
# ======================
sort -gk2 gwas_${PheName}_step2.pvalues > sort_min_pvalue.txt
