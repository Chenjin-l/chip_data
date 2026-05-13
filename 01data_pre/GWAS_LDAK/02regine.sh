

GENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/imp_QC"

snp="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/imp_QC.prune.in"

plink2 \
  --bfile ${GENO} \
  --extract ${snp} \
  --make-bed \
  --out fit_pheno_step1 \
  --threads 8
/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/04regeine/01bfile/fit_pheno_step1


#!/bin/bash
#PBS -N regenie_step1
#PBS -q batch
#PBS -l nodes=1:ppn=16
#PBS -l mem=50G
#PBS -o /home/Mzhou/02.F2/log                                                                                                           
#PBS -e /home/Mzhou/02.F2/log 
#PBS -j oe
#PBS -t 3-5

set -euo pipefail

#####################################
# 1. 激活环境
#####################################
source ~/anaconda3/etc/profile.d/conda.sh
conda activate regenie_env


#####################################
# 2. PBS ARRAY ID
#####################################
n=${PBS_ARRAYID}

#####################################
# 3. 文件
#####################################
PHENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/pheno_growth_ldak.rmOutlier.3sd.txt"

COVAR="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/covar_all.txt"

INDEX_FILE="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/index_growth.txt"

GENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/04regeine/01bfile/fit_pheno_step1"

cd /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/04regeine/02regeine/01step1
#####################################
# 4. 提取表型名
#####################################
trait=$(sed -n "${n}p" ${INDEX_FILE})

echo "=================================="
echo "PBS_ARRAYID : ${n}"
echo "Trait       : ${trait}"
echo "=================================="

#####################################
# 5. REGENIE STEP1
#####################################
regenie \
  --step 1 \
  --bed ${GENO} \
  --phenoFile ${PHENO} \
  --phenoCol ${trait} \
  --covarFile ${COVAR} \
  --covarColList age,PC1,PC2,PC3,PC4,PC5 \
  --catCovarList Sex \
  --qt \
  --apply-rint \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp_${trait} \
  --threads 32 \
  --out step1_${trait}

echo "DONE: ${trait}"


