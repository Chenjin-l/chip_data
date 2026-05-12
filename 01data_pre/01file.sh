#!/bin/bash
set -euo pipefail

############################
# 输入文件
############################
MAP_REF="/work/lingcj/03chip_data/01_data/01check/JXLab_ChinaChip50K_PLUS_SSC11.map"
MAP_OLD="/work/lingcj/03chip_data/01_data/01check/EE.map"
PED_OLD="/work/lingcj/03chip_data/01_data/01check/EE.ped"
KEEP_ID="/work/lingcj/03chip_data/02anlyis/01data/01pre_exchange/ID_compare/growth_keep.txt"

PREFIX="EE"

############################
# 输出目录
############################
WORKDIR="work_pipeline"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

############################
# 1. map转染色体+位置
############################
awk 'NR>1 {print $2"\t"$7"\t"$8}' ../${MAP_REF} > map11.txt

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{
    chr[$1]=$2
    pos[$1]=$3
    next
}
{
    id=$2
    if(id in chr){
        $1=chr[id]
        $4=pos[id]
    } else {
        $4="0"
    }
    print
}' map11.txt ../${MAP_OLD} > ${PREFIX}.SSC11.map

cp ../${PED_OLD} ${PREFIX}.SSC11.ped

############################
# 2. plink生成bfile + 样本过滤
############################
plink --file ${PREFIX}.SSC11 \
      --chr 1-18 \
      --keep ${KEEP_ID} \
      --make-bed \
      --out ${PREFIX}_growth

############################
# 3. QC
############################
plink --bfile ${PREFIX}_growth \
      --maf 0.01 \
      --geno 0.10 \
      --make-bed \
      --out ${PREFIX}_QC

############################
# 4. 转VCF（带压缩）
############################
plink --bfile ${PREFIX}_QC \
      --recode vcf-iid bgz \
      --out ${PREFIX}_QC

tabix -f ${PREFIX}_QC.vcf.gz

############################
# 5. 去重复（bcftools）
############################
bcftools norm -d exact \
    ${PREFIX}_QC.vcf.gz \
    -Oz -o ${PREFIX}_rmdup.vcf.gz

tabix -f -p vcf ${PREFIX}_rmdup.vcf.gz

############################
# 6. 重设ID（标准化）
############################
bcftools annotate \
    -x ID \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    ${PREFIX}_rmdup.vcf.gz \
    -Oz -o ${PREFIX}_final.vcf.gz

tabix -f -p vcf ${PREFIX}_final.vcf.gz

############################
# 7. 回写plink（二次确认一致性）
############################
plink --vcf ${PREFIX}_final.vcf.gz \
      --make-bed \
      --out ${PREFIX}_final

############################
# 8. 清理中间文件（只保留最终结果）
############################
rm -f \
map11.txt \
${PREFIX}.SSC11.map \
${PREFIX}.SSC11.ped \
${PREFIX}_growth.* \
${PREFIX}_QC.* \
${PREFIX}_rmdup.vcf.gz*

echo "DONE: final files kept:"
echo "  VCF: ${PREFIX}_final.vcf.gz + .tbi"
echo "  PLINK: ${PREFIX}_final.bed/.bim/.fam"
