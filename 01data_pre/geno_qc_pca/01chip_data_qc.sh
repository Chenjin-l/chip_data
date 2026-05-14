

PLINK v1.90b6.26 64-bit (2 Apr 2022)
Options in effect:
  --bfile /work/lingcj/03chip_data/02anlyis/01data/02bfile/EE
  --keep /work/lingcj/03chip_data/02anlyis/01data/01pre_exchange/ID_compare/FCR_id/FCR_keep.txt
  --make-bed
  --out EE_FCR

plink --bfile EE_growth \
  --maf 0.01 \
  --geno 0.10 \
  --hwe 1e-6 \
  --make-bed \
  --out EE_QC_growth


  
