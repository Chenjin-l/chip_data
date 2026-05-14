
三列，第一列是FID，第二列IID，第三列是表型数据y，没有行头，空格隔开。
#!/bin/bash
#PBS -N ldak_gwas
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -l mem=30G


GENO="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/imp_QC"
PHENO_COVAR=/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/04gcta_pheno
PHENO=${PHENO_COVAR}/pheno_growth_ldak.rmOutlier.3sd.txt
COVAR_QUANT="${PHENO_COVAR}/covar_quant.txt"
COVAR_FACTOR="${PHENO_COVAR}/covar_factor.txt"
GRM="/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/05gcta/01grm/gcta_grm"


gcta64 \
--bfile ${GENO} \
--make-grm \
--make-grm-alg 1 \
--out gcta_grm \
--threads 30

gcta64 \
--reml-bivar 1 2 \
--grm ${GRM} \
--pheno ${PHENO} \
--qcovar ${COVAR_QUANT} \
--covar ${COVAR_FACTOR} \
--threads 30 \
--out trait12_rg

gcta64 \
--reml-bivar 1 3 \
--grm ${GRM} \
--pheno ${PHENO} \
--qcovar ${COVAR_QUANT} \
--covar ${COVAR_FACTOR} \
--threads 30 \
--out trait13_rg

gcta64 \
--reml-bivar 2 3 \
--grm ${GRM} \
--pheno ${PHENO} \
--qcovar ${COVAR_QUANT} \
--covar ${COVAR_FACTOR} \
--threads 30 \
--out trait23_rg
