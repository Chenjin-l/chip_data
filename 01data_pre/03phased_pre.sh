plink --bfile qcfile --recode vcf-iid --out EE_QC
bgzip -c EE_QC.vcf > EE_QC.vcf.gz
tabix -f EE_QC.vcf.gz


###指控后的文件转为vcf然后拆分
#!/bin/bash

VCF="EE_QC.vcf.gz"
OUTDIR="split_vcf"

mkdir -p ${OUTDIR}

for chr in {1..18}; do

    echo "Processing chr${chr} ..."

    bcftools view \
        -r ${chr} \
        ${VCF} \
        -Oz \
        -o ${OUTDIR}/phased_chr${chr}.vcf.gz

    tabix -f -p vcf \
        ${OUTDIR}/phased_chr${chr}.vcf.gz

done

###################翻转
#!/bin/bash

REF_PANEL_DIR="/home/Mzhou/02.F2/04temp_data/01txk1000pig"

INPUT_VCF_PATH="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/01data/02vcf/split_vcf"

OUTPUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/01data/03conform-gt"

JAR="/home/Mzhou/02.F2/software/conform-gt/conform-gt.24May16.cee.jar"


for chr in {1..18}; do

    echo "========== chr${chr} =========="

    mkdir -p ${OUTPUT_DIR}/chr${chr}

    cd ${OUTPUT_DIR}/chr${chr}

    java -Xmx20g -jar ${JAR} \
        ref=${REF_PANEL_DIR}/ffxkt2_chr${chr}.vcf.gz \
        gt=${INPUT_VCF_PATH}/phased_chr${chr}.vcf.gz \
        chrom=${chr} \
        match=POS \
        out=mod.chr${chr} 

done
