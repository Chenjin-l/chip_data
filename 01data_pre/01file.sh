##转为11.1的位点
###############
awk 'NR>1 {print $2"\t"$7"\t"$8}' JXLab_ChinaChip50K_PLUS_SSC11.map > map11.txt
awk '
BEGIN{
    FS=OFS="\t"
}
NR==FNR{
    pos[$1]=$3
    chr[$1]=$2
    next
}
{
    id=$2
    if(id in pos){
        $1=chr[id]
        $4=pos[id]
    } else {
        $4="NA"
    }
    print
}
' map11.txt EE.map > EE.SSC11.map

###保留常染色体，提取需要保留的样本id
plink --file /work/lingcj/03chip_data/01_data/01check/data/EE.SSC11 \
      --chr 1-18 \
      --keep /work/lingcj/03chip_data/02anlyis/01data/01pre_exchange/ID_compare/growth_keep.txt \
      --make-bed \
      --out EE_growth

# =========================
# 1. 基础QC
# =========================

plink --bfile EE_growth \
  --maf 0.01 \
  --geno 0.10 \
  --hwe 1e-6 \
  --make-bed \
  --out EE_QC_growth

###转为vcf，改id
plink --bfile ../EE_QC_growth --recode vcf-iid --out EE_vcf

bcftools annotate \
  -x ID \
  --set-id '%CHROM\_%POS\_%REF\_%ALT' \
  EE_vcf.vcf -Ov -o EE_vcf.newID.vcf

plink --vcf EE_vcf.newID.vcf --make-bed --out growth_newid

