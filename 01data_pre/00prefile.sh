#!/bin/bash
set -euo pipefail

############################################################
# 基础参数
############################################################

MAP_REF="/work/lingcj/03chip_data/01_data/01check/JXLab_ChinaChip50K_PLUS_SSC11.map"

MAP_OLD="/work/lingcj/03chip_data/01_data/01check/EE.map"

PED_OLD="/work/lingcj/03chip_data/01_data/01check/EE.ped"

KEEP_ID="/work/lingcj/03chip_data/02anlyis/01data/01pre_exchange/ID_compare/growth_keep.txt"

REF_PANEL_DIR="/home/Mzhou/zPig_data/01pig_panel"

JAR="/home/Mzhou/02.F2/sofware/conform-gt/conform-gt.24May16.cee.jar"

PREFIX="EE"

############################################################
# 总输出目录
############################################################

OUTDIR="/work/lingcj/03chip_data/02anlyis/02_1KPIG_imputaion/01chipdata/01data"

mkdir -p ${OUTDIR}

cd ${OUTDIR}

############################################################
# 子目录
############################################################

WORKDIR="${OUTDIR}/00work"

SPLIT_DIR="${OUTDIR}/01split_vcf"

CONFORM_DIR="${OUTDIR}/02conform_gt"

MERGE_DIR="${OUTDIR}/03merge"

mkdir -p ${WORKDIR}
mkdir -p ${SPLIT_DIR}
mkdir -p ${CONFORM_DIR}
mkdir -p ${MERGE_DIR}

cd ${WORKDIR}

############################################################
# 1. map 转 SSC11
############################################################

echo "========== STEP1 map convert =========="

awk 'NR>1 {print $2"\t"$7"\t"$8}' ${MAP_REF} > map11.txt

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
}' map11.txt ${MAP_OLD} > ${PREFIX}.SSC11.map

cp ${PED_OLD} ${PREFIX}.SSC11.ped

############################################################
# 2. PED -> PLINK
############################################################

echo "========== STEP2 PED to PLINK =========="

plink --file ${PREFIX}.SSC11 \
      --chr 1-18 \
      --keep ${KEEP_ID} \
      --make-bed \
      --out ${PREFIX}_growth

############################################################
# 3. QC
############################################################

echo "========== STEP3 QC =========="

plink --bfile ${PREFIX}_growth \
      --maf 0.01 \
      --geno 0.10 \
      --make-bed \
      --out ${PREFIX}_QC

############################################################
# 4. 转VCF
############################################################

echo "========== STEP4 PLINK to VCF =========="

plink --bfile ${PREFIX}_QC \
      --recode vcf-iid bgz \
      --out ${PREFIX}_QC

tabix -f -p vcf ${PREFIX}_QC.vcf.gz

############################################################
# 5. 去重复
############################################################

echo "========== STEP5 Remove duplicate =========="

bcftools norm -d exact \
    ${PREFIX}_QC.vcf.gz \
    -Oz \
    -o ${PREFIX}_rmdup.vcf.gz

tabix -f -p vcf ${PREFIX}_rmdup.vcf.gz

############################################################
# 6. 染色体名称统一
############################################################

echo "========== STEP6 Rename chromosome =========="

cat > chr_name.txt <<EOF
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
EOF

bcftools annotate \
    --rename-chrs chr_name.txt \
    ${PREFIX}_rmdup.vcf.gz \
    -Oz \
    -o ${PREFIX}_rmdup_chr.vcf.gz

tabix -f -p vcf ${PREFIX}_rmdup_chr.vcf.gz

############################################################
# 7. 拆分染色体
############################################################

echo "========== STEP7 Split chromosome =========="

for chr in {1..18}; do

    echo "Split chr${chr}"

    bcftools view \
        -r chr${chr} \
        ${PREFIX}_rmdup_chr.vcf.gz \
        -Oz \
        -o ${SPLIT_DIR}/chr${chr}.vcf.gz

    tabix -f -p vcf ${SPLIT_DIR}/chr${chr}.vcf.gz

done

############################################################
# 8. conform-gt
############################################################

echo "========== STEP8 conform-gt =========="

for chr in {1..18}; do

    echo "========== chr${chr} =========="

    mkdir -p ${CONFORM_DIR}/chr${chr}

    cd ${CONFORM_DIR}/chr${chr}

    java -Xmx20g -jar ${JAR} \
        ref=${REF_PANEL_DIR}/chr${chr}.vcf.gz \
        gt=${SPLIT_DIR}/chr${chr}.vcf.gz \
        chrom=chr${chr} \
        match=POS \
        out=mod.chr${chr}

done

############################################################
# 9. 合并染色体
############################################################

echo "========== STEP9 Merge chromosome =========="

cd ${MERGE_DIR}

ls ${CONFORM_DIR}/chr*/mod.chr*.vcf.gz > merge.list

bcftools concat \
    -f merge.list \
    -Oz \
    -o ${PREFIX}_merged_conform_raw.vcf.gz

tabix -f -p vcf ${PREFIX}_merged_conform_raw.vcf.gz

############################################################
# 10. 标准化 variant ID
############################################################

echo "========== STEP10 Reset variant ID =========="

bcftools annotate \
    -x ID \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    ${PREFIX}_merged_conform_raw.vcf.gz \
    -Oz \
    -o ${PREFIX}_merged_conform.vcf.gz

tabix -f -p vcf ${PREFIX}_merged_conform.vcf.gz

############################################################
# 11. 转回 PLINK
############################################################

echo "========== STEP11 VCF to PLINK =========="

plink --vcf ${PREFIX}_merged_conform.vcf.gz \
      --make-bed \
      --out ${PREFIX}_merged_conform

############################################################
# 12. 清理中间文件（可选）
############################################################

echo "========== STEP12 Clean =========="

rm -f \
${WORKDIR}/map11.txt \
${WORKDIR}/${PREFIX}.SSC11.map \
${WORKDIR}/${PREFIX}.SSC11.ped

############################################################
# 完成
############################################################

echo "==========================================="
echo "PIPELINE FINISHED"
echo "==========================================="

echo ""
echo "Final VCF:"
echo "${MERGE_DIR}/${PREFIX}_merged_conform.vcf.gz"

echo ""
echo "Final PLINK:"
echo "${MERGE_DIR}/${PREFIX}_merged_conform.bed"
echo "${MERGE_DIR}/${PREFIX}_merged_conform.bim"
echo "${MERGE_DIR}/${PREFIX}_merged_conform.fam"

echo "==========================================="
