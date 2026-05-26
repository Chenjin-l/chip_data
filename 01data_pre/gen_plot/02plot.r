

========================================================
# 参数
#========================================================

name=4_backfat

GWAS=/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/02gwas/01LDAKgwas/01step1/4_backfat/gwas_backfat_step1.step2.assoc

bfile=/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/01pre/08merge/01merge/bfile/imp_QC

lead_file=/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/02gwas/01LDAKgwas/01step1/4_backfat/sort_min_pvalue.txt

#========================================================
# lead SNP
#========================================================

top_SNP=$(awk 'NR==2{print $1}' $lead_file)

echo "Lead SNP: $top_SNP"

#========================================================
# LD计算
#========================================================

plink \
  --bfile $bfile \
  --r2 \
  --ld-snp $top_SNP \
  --ld-window-kb 1000 \
  --ld-window 99999 \
  --ld-window-r2 0 \
  --nonfounders \
  --out ${name}.ld

#========================================================
# 提取 LD
#========================================================

awk '
BEGIN{OFS="\t"}
NR>1{
    print $(NF-1), $NF
}' ${name}.ld.ld > ${name}.r2

#========================================================
# GWAS + LD merge
#========================================================

join \
  <(sort -k1,1 ${name}.r2) \
  <(
    awk '
    BEGIN{OFS="\t"}
    NR>1{
        # CHR SNP POS P
        print $2,$1,$3,$7
    }' $GWAS \
    | sort -k1,1
  ) \
| awk '
BEGIN{
    OFS="\t";
    print "snp","chr","pos","p","r2"
}
{
    print $1,$3,$4,$5,$2
}' > ${name}.plotInput
