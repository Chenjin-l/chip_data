

###将指控，去重复的芯片 vcf文件拆分（去除重复）

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

#####匹配参考面板，矫正reference面板
#!/bin/bash                                                                                                                             
#PBS -N lcj_1                                                                                                                      
#PBS -q batch                                                                                                                           
#PBS -l nodes=1:ppn=8                                                                                                                 
#PBS -o /home/Mzhou/02.F2/log                                                                                                           
#PBS -e /home/Mzhou/02.F2/log                                                                                                           
#PBS -j oe                                                                                                                              
#PBS -l mem=30G

REF_PANEL_DIR="/home/Mzhou/zPig_data/01pig_panel"

INPUT_VCF_PATH="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/01data/02vcf/split_vcf"

OUTPUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/01data/03conform-gt"

JAR="/home/Mzhou/02.F2/sofware/conform-gt/conform-gt.24May16.cee.jar"

cd /work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/01data/03conform-gt
for chr in {1..18}; do
    echo "========== chr${chr} =========="

    mkdir -p ${OUTPUT_DIR}/chr${chr}

    cd ${OUTPUT_DIR}/chr${chr}

    java -Xmx20g -jar ${JAR} \
        ref=${REF_PANEL_DIR}/chr${chr}.vcf.gz \
        gt=${INPUT_VCF_PATH}/chr${chr}.vcf.gz \
        chrom=${chr} \
        match=POS \
        out=mod.chr${chr} 

done
