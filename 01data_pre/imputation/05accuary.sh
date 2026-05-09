#!/bin/bash
#PBS -N impute_acc
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l mem=10G
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -t 1-18

set -euo pipefail

chr=${PBS_ARRAYID}

echo "===== chr${chr} ====="

# =========================
# 路径设置
# =========================

BASE_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation"

TRUTH_VCF="${BASE_DIR}/02phased/phase.chr${chr}.vcf.gz"

MASK_SNP="${BASE_DIR}/03mask/chr${chr}.mask.snp"

IMPUTED_VCF="${BASE_DIR}/04imputation/chr${chr}/imputed_chr${chr}.qc.vcf.gz"

OUTDIR="${BASE_DIR}/05accuracy/chr${chr}"

mkdir -p ${OUTDIR}

cd ${OUTDIR}

# =========================
# 1. mask.snp → bed
# ID格式:
# chr_pos_ref_alt
# =========================

awk -F"_" '{
    print $1 "\t" ($2-1) "\t" $2
}' ${MASK_SNP} > mask.bed

# =========================
# 2. 提取 truth genotype
# =========================

bcftools view \
    --threads 4 \
    -R mask.bed \
    ${TRUTH_VCF} \
    -Oz \
    -o truth.vcf.gz

tabix -f -p vcf truth.vcf.gz

# =========================
# 3. 提取 imputed genotype
# =========================

bcftools view \
    --threads 4 \
    -R mask.bed \
    ${IMPUTED_VCF} \
    -Oz \
    -o imp.vcf.gz

tabix -f -p vcf imp.vcf.gz

# =========================
# 4. 转 dosage matrix
# =========================

plink2 \
    --vcf imp.vcf.gz \
    --export A \
    --out imp

plink2 \
    --vcf truth.vcf.gz \
    --export A \
    --out truth

# =========================
# 5. 拷贝R脚本
# =========================

cp ${BASE_DIR}/05accuracy/01accuary.r ./

# =========================
# 6. 运行R

~/anaconda3/envs/R4/bin/Rscript 01accuary.r > 01.log 2>&1 

echo "===== chr${chr} finished ====="
