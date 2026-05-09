


#!/bin/bash
#PBS -N concat
#PBS -q batch
#PBS -l nodes=1:ppn=16
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -l mem=30G

set -euo pipefail

INDIR=/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/04imputation

OUTDIR=/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# =========================
# 构建 + 数值排序
# =========================
find ${INDIR} -name "imputed_chr*.qc.vcf.gz" | sort -V > vcf.list

# =========================
# concat
# =========================

bcftools concat \
    -f vcf.list \
    --threads 16 \
    -Oz \
    -o merged_all_chr.vcf.gz

# =========================
# index
# =========================

tabix -f -p vcf merged_all_chr.vcf.gz

# =========================
# stats
# =========================

bcftools stats merged_all_chr.vcf.gz > stat.txt

echo "Finished merging all chromosomes"
