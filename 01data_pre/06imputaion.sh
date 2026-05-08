#!/bin/bash
#PBS -N beagle_impute
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l mem=60G
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -t 1-18

set -euo pipefail

BEAGLE=/home/Mzhou/02.F2/software/beagle.29Oct24.c8e.jar

REF_PANEL_DIR="/home/Mzhou/02.F2/04temp_data/01txk1000pig"

INPUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/03mask"

OUTPUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/04imputation"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${TMP_DIR}

chr=${PBS_ARRAYID}

echo "===== chr${chr} imputation ====="

IN_VCF=${INPUT_DIR}/chr${chr}.masked.vcf.gz

OUT_PREFIX=${OUTPUT_DIR}/chr${chr}/imputed_chr${chr}

mkdir -p ${OUTPUT_DIR}/chr${chr}

# =========================
# 1. BEAGLE imputation
# =========================
java -Xmx48g -jar ${BEAGLE} \
    ref=${REF_PANEL_DIR}/ffxkt2_chr${chr}.vcf.gz \
    gt=${IN_VCF} \
    chrom=${chr} \
    out=${OUT_PREFIX}

# =========================
# 2. index
# =========================
tabix -f -p vcf ${OUT_PREFIX}.vcf.gz

# =========================
# 3. DR2 filtering
# =========================
bcftools filter \
    -i 'INFO/DR2>=0.85' \
    ${OUT_PREFIX}.vcf.gz \
    -Oz -o ${OUTPUT_DIR}/chr${chr}/imputed_chr${chr}.qc.vcf.gz

# =========================
# 4. index filtered
# =========================
tabix -f -p vcf ${OUTPUT_DIR}/chr${chr}/imputed_chr${chr}.qc.vcf.gz

echo "chr${chr} done"
