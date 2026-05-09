#!/bin/bash
#PBS -N mask_vcf
#PBS -q batch
#PBS -l nodes=1:ppn=2
#PBS -l mem=10G
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -t 1-18

MASK_RATIO=0.05
SEED=222

INPUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/02phased"

OUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/03mask"

mkdir -p ${OUT_DIR}

chr=${PBS_ARRAYID}

echo "===== chr${chr} ====="

VCF=${INPUT_DIR}/phase.chr${chr}.vcf.gz

ALL_SNP=${OUT_DIR}/chr${chr}.all.snp
MASK_SNP=${OUT_DIR}/chr${chr}.mask.snp
OUT_VCF=${OUT_DIR}/chr${chr}.masked.vcf.gz

# =========================
# 1. 提取 SNP ID
# =========================
bcftools query -r ${chr} -f '%ID\n' ${VCF} > ${ALL_SNP}

TOTAL=$(wc -l < ${ALL_SNP})
MASK_N=$(awk -v n=$TOTAL -v r=$MASK_RATIO 'BEGIN{printf "%d", n*r}')

echo "Total=$TOTAL Mask=$MASK_N"

# =========================
# 2. 随机抽 SNP
# =========================
awk -v seed=$((SEED+chr)) '
BEGIN{srand(seed)}
{print rand(), $1}
' ${ALL_SNP} \
| sort -k1,1n \
| head -n ${MASK_N} \
| cut -d" " -f2 \
> ${MASK_SNP}

# =========================
# 3. mask SNP
# =========================
bcftools view \
    -r ${chr} \
    -e 'ID=@'"${MASK_SNP}" \
    -Oz \
    -o ${OUT_VCF} \
    ${VCF}

# =========================
# 4. index
# =========================
bcftools index -f ${OUT_VCF}



