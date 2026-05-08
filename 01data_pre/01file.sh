
plink --file /work/lingcj/03chip_data/01_data/01check/data/EE.SSC11 \
      --chr 1-18 \
      --keep /work/lingcj/03chip_data/02anlyis/01data/01pre_exchange/ID_compare/growth_keep.txt \
      --make-bed \
      --out EE_growth

# =========================
# 1. 基础QC（芯片数据）
# =========================

plink --bfile EE_growth \
  --maf 0.01 \
  --geno 0.10 \
  --hwe 1e-6 \
  --make-bed \
  --out EE_QC_growth

mkdir -p PCA
plink --bfile EE_QC_growth --pca 20 --out PCA/EE_growth_pca

source ~/anaconda3/etc/profile.d/conda.sh                                                                                               
conda activate R4 
Rscript pca_plot.r > plot.log

plink --bfile ../EE_QC_growth --recode vcf-iid --out EE_vcf

bcftools annotate \
  -x ID \
  --set-id '%CHROM\_%POS\_%REF\_%ALT' \
  EE_vcf.vcf -Ov -o EE_vcf.newID.vcf

plink --vcf EE_vcf.newID.vcf --make-bed --out growth_newid

