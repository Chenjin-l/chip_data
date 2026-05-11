#!/bin/bash
#PBS -N ldak_mixed
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l mem=30G
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -t 3-5

set -euo pipefail

# =====================================
# 环境
# =====================================

source ~/anaconda3/etc/profile.d/conda.sh
conda activate ldak_env

THREADS=8

# =====================================
# 数据
# =====================================

GENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/newbfile/new"

PHENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/pheno_growth_ldak.rmOutlier.3sd.txt"

COVAR_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar"

COVAR_QUANT="${COVAR_DIR}/covar_quant.txt"

COVAR_FACTOR="${COVAR_DIR}/covar_factor.txt"

INDEX_FILE="${COVAR_DIR}/index_growth.txt"

GRM="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/02LDAK_gwas/GRM/le"

WORK_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/02LDAK_gwas/02_mixed_model"

# =====================================
# phenotype
# =====================================

n=${PBS_ARRAYID}

pheno_col=$((n-2))

PheName=$(sed -n "${n}p" ${INDEX_FILE})

# =====================================
# 工作目录
# =====================================

mkdir -p ${WORK_DIR}/${n}_${PheName}

cd ${WORK_DIR}/${n}_${PheName}

# =====================================
# Mixed model GWAS
# =====================================

echo "Running mixed model GWAS for ${PheName}"

ldak6 \
--linear ${PheName}_mixed \
--bfile ${GENO} \
--pheno ${PHENO} \
--mpheno ${pheno_col} \
--covar ${COVAR_QUANT} \
--factors ${COVAR_FACTOR} \
--grm ${GRM} \
--max-threads ${THREADS}

# =====================================
# 排序P值
# =====================================

sort -gk2 ${PheName}_mixed.pvalues > ${PheName}_mixed.sorted.pvalues
