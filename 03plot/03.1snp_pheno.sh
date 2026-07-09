#!/bin/bash
#PBS -N SNP_Genotype_Plot
#PBS -q batch
#PBS -l nodes=1:ppn=2
#PBS -l mem=4G
#PBS -o /home/Mzhou/zlingcj/log
#PBS -e /home/Mzhou/zlingcj/log
#PBS -j oe

set -e

echo "Job started at: $(date)"

BASE="/home/Mzhou/zlingcj/03chip_data/02anlyis/05_mascript"
VCF="/home/Mzhou/zlingcj/03chip_data/02anlyis/05_mascript/01pheno/01data/03geno_bfile/00vcf/merged_onlyGT.vcf.gz"
PHENO_DIR="${BASE}/01pheno/01data/02pheno" #提供表型文件
OUTPUT_DIR="${BASE}/04snpgenotype_pheno/result"

#SNP_ID="11_7690019_A_T"
SNP_ID="13_15464157_A_G"
mkdir -p ${OUTPUT_DIR}/${SNP_ID}
cd ${OUTPUT_DIR}/${SNP_ID}

SEX_GROUPS=("all" "female" "male")

declare -A PHENO_FILE=(
    ["all"]="${PHENO_DIR}/all_pheno_backfat.txt"
    ["female"]="${PHENO_DIR}/female_pheno.txt"
    ["male"]="${PHENO_DIR}/male_pheno.txt"
)

PHENO_NAME="Backfat thickness"

echo "========================================"
echo "SNP: ${SNP_ID}"
echo "========================================"

# ============================================================
# 1. 提取SNP
# ============================================================
if [[ -f ${SNP_ID}.vcf ]]; then
    echo "SNP VCF already exists, skipping extraction..."
else
    echo "Extracting SNP from VCF..."
    zgrep -w "${SNP_ID}" ${VCF} > ${SNP_ID}.vcf
    if [[ ! -s ${SNP_ID}.vcf ]]; then
        echo "Error: SNP ${SNP_ID} not found in VCF!"
        exit 1
    fi
fi

# ============================================================
# 2. 提取样本ID和基因型
# ============================================================
if [[ -f all_${SNP_ID}.genotype ]]; then
    echo "Genotype file already exists, skipping extraction..."
else
    echo "Extracting sample IDs and genotypes..."
    zgrep "^#CHROM" ${VCF} | cut -f10- | tr '\t' '\n' > all_sample.id
    awk '{for(i=10;i<=NF;i++){split($i,a,":"); print a[1]}}' ${SNP_ID}.vcf > all_genotype.txt
    paste all_sample.id all_genotype.txt > all_${SNP_ID}.genotype
fi

# ============================================================
# 3. 为每个群体筛选基因型
# ============================================================
for SEX in "${SEX_GROUPS[@]}"; do
    echo "----------------------------------------"
    echo "Processing: ${SEX}"
    
    PHENO="${PHENO_FILE[${SEX}]}"
    
    # 提取该群体的样本ID（从表型文件第2列）
    cut -f2 ${PHENO} > ${SEX}_sample.id
    
    # 筛选基因型
    awk 'NR==FNR{a[$1];next} $1 in a' ${SEX}_sample.id all_${SNP_ID}.genotype > ${SNP_ID}_${SEX}.genotype
    
    n_samples=$(wc -l < ${SNP_ID}_${SEX}.genotype)
    echo "${SEX} samples: ${n_samples}"
done

# ============================================================
# 4. 生成R脚本（单独PDF）
# ============================================================
cat > plot_${SNP_ID}.R << 'EOF_R'
library(data.table)
library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly=TRUE)
snp_id <- args[1]
pheno_name <- args[2]

# ============================================================
# 群体配置
# ============================================================
groups <- c("all", "female", "male")
pheno_files <- c(
    "/home/Mzhou/zlingcj/03chip_data/02anlyis/05_mascript/01pheno/01data/02pheno/all_pheno_backfat.txt",
    "/home/Mzhou/zlingcj/03chip_data/02anlyis/05_mascript/01pheno/01data/02pheno/female_pheno.txt",
    "/home/Mzhou/zlingcj/03chip_data/02anlyis/05_mascript/01pheno/01data/02pheno/male_pheno.txt"
)
group_labels <- c("All", "Female", "Male")

all_summary <- data.table()

for (i in seq_along(groups)) {
    group <- groups[i]
    pheno_file <- pheno_files[i]
    
    cat("\n========================================\n")
    cat("Processing:", group, "\n")
    
    # ============================================================
    # 读取基因型
    # ============================================================
    geno_file <- paste0(snp_id, "_", group, ".genotype")
    if (!file.exists(geno_file)) {
        cat("Warning:", geno_file, "not found, skipping\n")
        next
    }
    
    geno <- fread(geno_file, header=FALSE, sep="\t", quote="",
                  colClasses=c("character", "character"))
    setnames(geno, c("IID", "Genotype"))
    
    cat("Genotype samples:", nrow(geno), "\n")
    cat("Genotypes:", paste(unique(geno$Genotype), collapse=", "), "\n")
    
    # ============================================================
    # 读取表型（灵活处理）
    # ============================================================
        if (group == "all") {
        pheno <- fread(pheno_file, header=FALSE, sep=" ", fill=TRUE, quote="")
    } else {
        pheno <- fread(pheno_file, header=FALSE, sep="\t", fill=TRUE, quote="")
    }
    n_cols <- ncol(pheno)
    cat("Phenotype file has", n_cols, "columns\n")
    
    if (n_cols >= 3) {
        pheno <- pheno[, 1:3, with=FALSE]
        setnames(pheno, c("FID", "IID", "Pheno"))
    } else if (n_cols == 2) {
        setnames(pheno, c("IID", "Pheno"))
        pheno[, FID := IID]
    } else {
        cat("Error: Unexpected number of columns in phenotype file:", n_cols, "\n")
        next
    }
    
    cat("Phenotype samples:", nrow(pheno), "\n")
    
    # ============================================================
    # 检查ID匹配
    # ============================================================
    common_iid <- intersect(geno$IID, pheno$IID)
    common_fid <- intersect(geno$IID, pheno$FID)
    cat("Common IDs (geno vs pheno.IID):", length(common_iid), "\n")
    cat("Common IDs (geno vs pheno.FID):", length(common_fid), "\n")
    
    # ============================================================
    # 合并
    # ============================================================
    if (length(common_iid) > 0) {
        dat <- merge(geno, pheno, by="IID")
    } else if (length(common_fid) > 0) {
        dat <- merge(geno, pheno, by.x="IID", by.y="FID")
    } else {
        cat("Warning: No matching IDs for", group, "\n")
        next
    }
    
    cat("Merged samples:", nrow(dat), "\n")
    
    if (nrow(dat) == 0) {
        cat("Warning: Merge produced 0 rows for", group, "\n")
        next
    }
    
    # ============================================================
    # 清理
    # ============================================================
    valid_geno <- c("0|0", "0|1", "1|0", "1|1")
    dat <- dat[Genotype %in% valid_geno, ]
    dat <- dat[!is.na(Pheno), ]
    
    if (nrow(dat) == 0) {
        cat("Warning: No valid genotypes for", group, "\n")
        next
    }
    
    # 合并 0|1 和 1|0
    if ("0|1" %in% dat$Genotype && "1|0" %in% dat$Genotype) {
        dat[Genotype == "1|0", Genotype := "0|1"]
    }
    dat[, Genotype := factor(Genotype, levels=c("0|0", "0|1", "1|1"))]
    dat <- dat[!is.na(Genotype), ]
    
    cat("Final samples:", nrow(dat), "\n")
    print(table(dat$Genotype))
    
    # ============================================================
    # 保存合并数据
    # ============================================================
    write.table(dat, file=paste0(snp_id, "_", group, "_merged.txt"),
                sep="\t", quote=FALSE, row.names=FALSE)
    
    # ============================================================
    # 统计摘要
    # ============================================================
    summary_dt <- dat[, .(N = .N, Mean = round(mean(Pheno), 4), SD = round(sd(Pheno), 4)), by = Genotype]
    summary_dt[, Group := group]
    all_summary <- rbind(all_summary, summary_dt)
    
    # ============================================================
    # 绘图（单独保存）
    # ============================================================
    if (length(unique(dat$Genotype)) >= 2) {
        existing <- unique(dat$Genotype)
        comps <- list()
        if (all(c("0|0","0|1") %in% existing)) comps <- c(comps, list(c("0|0","0|1")))
        if (all(c("0|0","1|1") %in% existing)) comps <- c(comps, list(c("0|0","1|1")))
        if (all(c("0|1","1|1") %in% existing)) comps <- c(comps, list(c("0|1","1|1")))
        
        p <- ggplot(dat, aes(x=Genotype, y=Pheno)) +
            geom_boxplot(outlier.shape=16, outlier.size=2, fill="white") +
            stat_compare_means(comparisons=comps, method="wilcox.test",
                               label="p.signif", size=5) +
            labs(
                title = paste0(snp_id, " (", group_labels[i], ")"),
                x = "Genotype",
                y = pheno_name
            ) +
            theme_bw(base_size=14) +
            theme(
                plot.title = element_text(hjust = 0.5, size=16, face="bold"),
                panel.grid = element_blank()
            )
        
        # 单独保存每个群体的PDF
        pdf_file <- paste0(snp_id, "_", group, "_plot.pdf")
        ggsave(pdf_file, p, width=6, height=5)
        cat("Plot saved:", pdf_file, "\n")
    } else {
        cat("Not enough genotype groups for", group, "\n")
    }
}

# ============================================================
# 保存统计摘要
# ============================================================
if (nrow(all_summary) > 0) {
    write.table(all_summary, file=paste0(snp_id, "_summary.txt"),
                sep="\t", quote=FALSE, row.names=FALSE)
    cat("\nSummary saved:", paste0(snp_id, "_summary.txt"), "\n")
    print(all_summary)
}

cat("\nAll done!\n")
EOF_R

# ============================================================
# 5. 运行R
# ============================================================
source ~/anaconda3/etc/profile.d/conda.sh
conda activate R4

Rscript plot_${SNP_ID}.R ${SNP_ID} "${PHENO_NAME}"

echo ""
echo "========================================"
echo "Done!"
echo "Results saved in: ${OUTPUT_DIR}/${SNP_ID}/"
echo "========================================"
echo "Job finished at: $(date)"
