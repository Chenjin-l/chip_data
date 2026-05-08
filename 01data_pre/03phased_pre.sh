

REF_PANEL_DIR="/home/Mzhou/02.F2/04temp_data/01txk1000pig"
INPUT_VCF_path=""
OUTPUT_DIR=""

# 创建输出目录
mkdir -p ${OUTPUT_DIR}/chr${chr}
cd ${OUTPUT_DIR}/chr${chr}

# =========================
# conform-gt
# =========================

java -jar conform-gt.24May16.cee.jar ref=chr20.1kg.phase3.v5a.vcf.gz gt=chr20.vcf.gz chrom=20 out=mod.chr20 excludesamples=non.eur.excl
