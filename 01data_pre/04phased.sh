
#!/bin/bash
#PBS -N beagle_phase
#PBS -q batch
#PBS -l nodes=1:ppn=10
#PBS -l mem=30G
#PBS -o /home/Mzhou/02.F2/log 
#PBS -e /home/Mzhou/02.F2/log
#PBS -j oe
#PBS -t 1-5

BEAGLE=/home/Mzhou/02.F2/software/beagle.29Oct24.c8e.jar

REF_DIR="/home/Mzhou/02.F2/04temp_data/01txk1000pig"

INPUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/01data/03conform-gt"

OUT_DIR="/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/02phased"

mkdir -p ${OUT_DIR}

chr=${PBS_ARRAYID}

    echo "===== chr${chr} ====="

    java -Xmx24g -jar ${BEAGLE} \
        gt=${INPUT_DIR}/chr${chr}/mod.chr${chr}.vcf.gz \
        ref=${REF_DIR}/ffxkt2_chr${chr}.vcf.gz \
        out=${OUT_DIR}/phase.chr${chr} \
        chrom=${chr} \
        nthreads=8 \
        impute=false

