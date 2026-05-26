

#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=4
#PBS -l mem=20G
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -t 1-18

set -euo pipefail

echo "Job started at:" `date`

########################################
# 工作目录
########################################
Work_dir=/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/03sexgwas/03gemma/03result/female

cd ${Work_dir}

########################################
# 染色体编号
########################################
chr=${PBS_ARRAYID}

########################################
# 表型编号（手动修改）
########################################
i=3

INDEX="/home/Mzhou/zlingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/index_growth.txt"

name=$(awk -v n=$i 'NR==n{print $1}' ${INDEX})

echo "Trait : ${i} ${name}"
echo "Chr   : ${chr}"

########################################
# 输入文件
########################################
BFILE=/home/Mzhou/zlingcj/03chip_data/02anlyis/01data/03growth/02pheno/04new_adjust/02sex/01pre/sex_split/split_chr/imp_QC_female_chr${chr}

GRM=/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/03sexgwas/03gemma/01GRM/output/female_grm.sXX.txt

PHENO=/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/03sexgwas/03gemma/02pheno/pheno_growth_female_normal.txt


########################################
# 临时目录
########################################
temp=${Work_dir}/${i}_${name}

mkdir -p ${temp}

cd ${temp}

########################################
# GEMMA GWAS
########################################
gemma \
-bfile ${BFILE} \
-k ${GRM} \
-lmm 1 \
-p ${PHENO} \
-n ${i} \
-o chr${chr}_${i}

echo "Finished chr${chr}"
