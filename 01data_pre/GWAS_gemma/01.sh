gemma-0.98.1-linux-static -bfile c -gk 2 -p p.txt 
代码解释：

-bfile 读取plink的二进制文件
-gk 2 标准化的方法计算G矩阵
-p 读取表型数据（这个不能省掉）

gemma-0.98.1-linux-static -bfile c -k output/result.sXX.txt -lmm 1 -p p.txt 
代码解释：
* -bfile plink的二进制文件
* -k 读取G矩阵的文件
* -lmm 1 使用Wald的方法进行SNP检验
* -p 表型数据

######先生成 GRM（只做一次

gemma \
-bfile /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/01bfile/work_pipeline/EE_final \
-gk 2 \
-p /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/02pheno/gemma_pheno.txt \
-n 3 \
-o EE_grm


#######
#!/bin/bash
#PBS -N gemma_gwas
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l mem=50G
#PBS -o /home/Mzhou/02.F2/log
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -t 3-5

set -euo pipefail
cd /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/04result
i=${PBS_ARRAYID}

INDEX=/work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/02pheno/index.txt

name=$(awk -v n=$i 'NR==n{print $1}' ${INDEX})

echo "Running trait: ${i} ${name}"

gemma \
-bfile /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/01bfile/work_pipeline/EE_final \
-k /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/03GRM/output/result.sXX.txt \
-lmm 1 \
-p /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/02pheno/gemma_pheno.txt \
-n ${i} \
-c /work/lingcj/03chip_data/02anlyis/01data/03growth/03GWAS/01chip_gwas/03Gemma_gwas/02pheno/gemma.covar \
-o ${i}_${name}

echo "Finished: ${i}_${name}"
