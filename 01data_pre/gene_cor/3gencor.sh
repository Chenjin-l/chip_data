#!/bin/bash
#PBS -N gencor
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l mem=30G
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe

set -euo pipefail

source ~/anaconda3/etc/profile.d/conda.sh
conda activate ldak_env

THREADS=8

# =====================================
# 基础路径
# =====================================

BASE="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/02LDAK_gwas"

TAG="${BASE}/03gencor/01tagg/pig_tag.tagging"

SUMDIR="${BASE}/02_mixed_model"

OUTDIR="${BASE}/03gencor/02result"

mkdir -p ${OUTDIR}

cd ${OUTDIR}

# =====================================
# weight vs backfat
# =====================================

ldak6 \
--sum-cors weight_vs_backfat \
--summary  ${SUMDIR}/3_weight/weight_mixed.summaries \
--summary2 ${SUMDIR}/4_backfat/backfat_mixed.summaries \
--tagfile ${TAG} \
--allow-ambiguous YES \
--check-sums NO \
--max-threads ${THREADS}

# =====================================
# weight vs eye_muscle
# =====================================

ldak6 \
--sum-cors weight_vs_eye_muscle \
--summary  ${SUMDIR}/3_weight/weight_mixed.summaries \
--summary2 ${SUMDIR}/5_eye_muscle/eye_muscle_mixed.summaries \
--tagfile ${TAG} \
--allow-ambiguous YES \
--check-sums NO \
--max-threads ${THREADS}

# =====================================
# backfat vs eye_muscle
# =====================================

ldak6 \
--sum-cors backfat_vs_eye_muscle \
--summary  ${SUMDIR}/4_backfat/backfat_mixed.summaries \
--summary2 ${SUMDIR}/5_eye_muscle/eye_muscle_mixed.summaries \
--tagfile ${TAG} \
--allow-ambiguous YES \
--check-sums NO \
--max-threads ${THREADS}
