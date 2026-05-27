

#!/bin/bash
#PBS -N female_gemma_chr
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























#!/bin/bash
#PBS-N gemma_merge
#PBS-q batch
#PBS-l nodes=1:ppn=4
#PBS-l mem=10G
#PBS-o /home/Mzhou/zlingcj/log
#PBS-e /home/Mzhou/zlingcj/log
#PBS-j oe
#PBS-t 3-5
set -euo pipefail

echo "Start: $(date)"

########################################
# 表型编号（手动修改）
########################################
i=${PBS_ARRAYID}

INDEX=/home/Mzhou/zlingcj/03chip_data/02anlyis/01data/03growth/03GWAS/02imputation/02_1kpig_imp_gwas/01data/02pheno_covar/index_growth.txt

name=$(awk -v n=$i 'NR==n{print $1}' ${INDEX})

echo "Trait ID : ${i}"
echo "Trait name: ${name}"

########################################
# 路径
########################################
BASE=/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/03sexgwas/03gemma/03result/male

WORKDIR=${BASE}/${i}_${name}
OUTDIR=${WORKDIR}/output

cd ${WORKDIR}

########################################
# 合并 chr1-18
########################################
echo "Merging..."

cat ${OUTDIR}/chr*_${i}.assoc.txt \
| grep -v "^#" \
| sort -gk 12,12 -S 20% --parallel=30 \
> ${name}_${i}.assoc.txt

sed -i '1,17d' ${name}_${i}.assoc.txt
echo "Done: ${name}_${i}.assoc.txt"
echo "End: $(date)"
02merge.sh (END)







