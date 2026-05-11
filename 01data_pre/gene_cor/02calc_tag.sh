#!/bin/bash
#PBS -N ldak_tag
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

GENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/newbfile/new"

OUTDIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/02LDAK_gwas/03gencor/01tagg"

mkdir -p ${OUTDIR}

cd ${OUTDIR}

ldak6 \
--calc-tagging pig_tag \
--bfile ${GENO} \
--power -0.25 \
--window-kb 1000 \
--save-matrix YES \
--max-threads ${THREADS}
