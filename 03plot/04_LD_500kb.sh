#!/bin/bash
#PBS -N LCJ_LD_Plot
#PBS -q batch
#PBS -l nodes=1:ppn=4
#PBS -l mem=10G
#PBS -o /home/Mzhou/zlingcj/log
#PBS -e /home/Mzhou/zlingcj/log
#PBS -j oe

set -e

echo "Job started at: $(date)"

# ============================================================
# 路径设置
# ============================================================
BASE="/home/Mzhou/zlingcj/03chip_data/02anlyis/05_mascript"
GWAS_DIR="${BASE}/02gwas"
BFILE_DIR="${BASE}/01pheno/01data/03geno_bfile/qc_filtered"
OUTPUT_DIR="/home/Mzhou/zlingcj/03chip_data/02anlyis/05_mascript/06LD_Plot/result"

mkdir -p ${OUTPUT_DIR}

# ============================================================
# 定义三个群体
# ============================================================
declare -A GWAS_FILE=(
    ["all"]="${GWAS_DIR}/all/all_merged.fastGWA"
    ["female"]="${GWAS_DIR}/female/female.fastGWA"
    ["male"]="${GWAS_DIR}/male/male.fastGWA"
)

declare -A BFILE_PREFIX=(
    ["all"]="${BFILE_DIR}/all_merged_qc_LD"
    ["female"]="${BFILE_DIR}/female_qc_LD"
    ["male"]="${BFILE_DIR}/male_qc_LD"
)

# ============================================================
# 测试位点
# ============================================================
#TOP_SNP="11_7690019_A_T"
#TOP_SNP="13_15464157_A_G"
#TOP_SNP="12_14671275_C_T"
#TOP_SNP="12_14650591_A_G"
TOP_SNP="12_14671275_C_T"
WINDOW=500000
Y_MAX=10

echo "========================================"
echo "Testing SNP: ${TOP_SNP}"
echo "Window: ±${WINDOW} bp"
echo "Y-axis max: ${Y_MAX}"
echo "========================================"

# 提取该SNP的位置（从all群体GWAS获取）
GWAS_ALL="${GWAS_FILE[all]}"
SNP_CHR=$(awk -v snp="${TOP_SNP}" '$2==snp{print $1}' ${GWAS_ALL})
SNP_POS=$(awk -v snp="${TOP_SNP}" '$2==snp{print $3}' ${GWAS_ALL})

if [[ -z ${SNP_CHR} || -z ${SNP_POS} ]]; then
    echo "Error: ${TOP_SNP} not found in GWAS"
    exit 1
fi

X_MIN=$((SNP_POS - WINDOW))
X_MAX=$((SNP_POS + WINDOW))

#X_MIN=14650000
#X_MAX=15950000


echo "SNP: ${TOP_SNP}"
echo "Chr: ${SNP_CHR}"
echo "Pos: ${SNP_POS}"
echo "X-axis range: ${X_MIN} - ${X_MAX}"

# ============================================================
# 循环三个群体
# ============================================================
for group in all female male; do
    echo ""
    echo "========================================"
    echo "Processing: ${group}"
    echo "========================================"
    
    GWAS="${GWAS_FILE[$group]}"
    BFILE="${BFILE_PREFIX[$group]}"
    NAME="${group}_${TOP_SNP}"
    
    echo "GWAS: ${GWAS}"
    echo "BFILE: ${BFILE}"
    echo "Output: ${OUTPUT_DIR}/${NAME}"
    
    if [[ ! -f ${GWAS} ]]; then
        echo "Warning: GWAS file not found: ${GWAS}"
        continue
    fi
    
    if [[ ! -f ${BFILE}.bed ]]; then
        echo "Warning: BFILE not found: ${BFILE}"
        continue
    fi
    
    cd ${OUTPUT_DIR}
    
    # ============================================================
    # LD计算
    # ============================================================
    echo "Calculating LD for ${TOP_SNP} in ${group}..."
    
    plink \
        --bfile ${BFILE} \
        --r2 \
        --ld-snp ${TOP_SNP} \
        --ld-window-kb 1000 \
        --ld-window 99999 \
        --ld-window-r2 0 \
        --nonfounders \
        --out ${NAME}
    
    # ============================================================
    # 提取 LD r2
    # ============================================================
    awk '
    BEGIN{OFS="\t"}
    NR>1{
        print $(NF-1), $NF
    }' ${NAME}.ld > ${NAME}.r2
    
    # ============================================================
    # GWAS + LD merge（只保留窗口内SNP）
    # ============================================================
    join \
        <(sort -k1,1 ${NAME}.r2) \
        <(
            awk -v pos=${SNP_POS} -v win=${WINDOW} '
            BEGIN{OFS="\t"}
            NR>1{
                if ($3 >= pos-win && $3 <= pos+win) {
                    print $2, $1, $3, $10
                }
            }' ${GWAS} \
            | sort -k1,1
        ) \
    | awk '
    BEGIN{
        OFS="\t";
        print "snp","chr","pos","p","r2"
    }
    {
        print $1, $3, $4, $5, $2
    }' > ${NAME}.plotInput
    
    # ============================================================
    # R 绘图脚本（固定坐标：x轴固定窗口，y轴固定为8）
    # ============================================================
    cat > ${NAME}.locusplot.R << 'EOF_R'
    library(data.table)
    library(ggplot2)
    
    args <- commandArgs(trailingOnly = TRUE)
     name_prefix <- args[1]
    snp_id <- args[2]
    group_name <- args[3]
    x_min <- as.numeric(args[4])
    x_max <- as.numeric(args[5])
    y_max <- as.numeric(args[6])
    
    dat <- fread(paste0(name_prefix, ".plotInput"))
    
    dat[, pos := as.numeric(pos)]
    dat[, r2  := as.numeric(r2)]
    dat[, p   := as.numeric(p)]
    
    dat <- dat[!is.na(p) & !is.na(r2), ]
    
    # LD 分组
    dat[, LD_group := fifelse(r2 > 0 & r2 < 0.2,  "<0.2",
                       fifelse(r2 < 0.4, "0.2-0.4",
                       fifelse(r2 < 0.6, "0.4-0.6",
                       fifelse(r2 < 0.8, "0.6-0.8",
                       fifelse(r2 < 1,   "0.8-1",
                       fifelse(r2 == 1,  "1", NA_character_))))))]
    
    col_map <- c(
        "1"       = "#A020F0",
        "0.8-1"   = "#FD0000",
        "0.6-0.8" = "#FEA500",
        "0.4-0.6" = "#00FF00",
        "0.2-0.4" = "#82D6FE",
        "<0.2"    = "#00007C"
    )
    
    dat[, LD_group := factor(LD_group,
                             levels = c("1", "0.8-1", "0.6-0.8",
                                        "0.4-0.6", "0.2-0.4", "<0.2"))]
    
    lead <- dat[snp == snp_id]   # 固定位点，紫色菱形

# ============================================================
# 计算阈值线：-log10(1/352479) = log10(352479) ≈ 5.547
# ============================================================
n_tests <- 352479
threshold <- -log10(1 / n_tests)   # 或直接写 log10(n_tests)
threshold_label <- paste0("-log10(1/", n_tests, ") = ", round(threshold, 2))

    
    p <- ggplot(dat, aes(x = pos, y = -log10(p))) +
        geom_point(aes(color = LD_group), size = 3, shape = 16) +
        geom_point(data = lead,
                   aes(x = pos, y = -log10(p)),
                   shape = 18, size = 5, color = "#A020F0") +
    geom_hline(yintercept = threshold,
               linetype = "dashed",
               color = "black",
    name_prefix <- args[1]
    snp_id <- args[2]
    group_name <- args[3]
    x_min <- as.numeric(args[4])
    x_max <- as.numeric(args[5])
    y_max <- as.numeric(args[6])
    
    dat <- fread(paste0(name_prefix, ".plotInput"))
    
    dat[, pos := as.numeric(pos)]
    dat[, r2  := as.numeric(r2)]
    dat[, p   := as.numeric(p)]
    
    dat <- dat[!is.na(p) & !is.na(r2), ]
    
    # LD 分组
    dat[, LD_group := fifelse(r2 > 0 & r2 < 0.2,  "<0.2",
                       fifelse(r2 < 0.4, "0.2-0.4",
                       fifelse(r2 < 0.6, "0.4-0.6",
                       fifelse(r2 < 0.8, "0.6-0.8",
                       fifelse(r2 < 1,   "0.8-1",
                       fifelse(r2 == 1,  "1", NA_character_))))))]
    
    col_map <- c(
        "1"       = "#A020F0",
        "0.8-1"   = "#FD0000",
        "0.6-0.8" = "#FEA500",
        "0.4-0.6" = "#00FF00",
        "0.2-0.4" = "#82D6FE",
        "<0.2"    = "#00007C"
    )
    
    dat[, LD_group := factor(LD_group,
                             levels = c("1", "0.8-1", "0.6-0.8",
                                        "0.4-0.6", "0.2-0.4", "<0.2"))]
    
    lead <- dat[snp == snp_id]   # 固定位点，紫色菱形

# ============================================================
# 计算阈值线：-log10(1/352479) = log10(352479) ≈ 5.547
# ============================================================
n_tests <- 352479
threshold <- -log10(1 / n_tests)   # 或直接写 log10(n_tests)
threshold_label <- paste0("-log10(1/", n_tests, ") = ", round(threshold, 2))

    
    p <- ggplot(dat, aes(x = pos, y = -log10(p))) +
        geom_point(aes(color = LD_group), size = 3, shape = 16) +
        geom_point(data = lead,
                   aes(x = pos, y = -log10(p)),
                   shape = 18, size = 5, color = "#A020F0") +
    geom_hline(yintercept = threshold,
               linetype = "dashed",
               color = "black",
               linewidth = 0.8) +
    annotate("text",
             x = x_min + (x_max - x_min) * 0.02,
             y = threshold + 0.2,
             label = threshold_label,
             color = "red",
             size = 3.5,
             hjust = 0) +
        scale_color_manual(values = col_map, name = expression(LD~r^2)) +
        coord_cartesian(
            xlim = c(x_min, x_max),
            ylim = c(0, y_max),
            expand = FALSE
        ) +
        scale_x_continuous(
            name = "Position (bp)",
            breaks = seq(x_min, x_max, length.out = 6),
            labels = function(x) format(x, big.mark = ","),
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            name = "-log10(P)",
            breaks = seq(0, y_max, by = 2),
            expand = c(0.05, 0.05)
        ) +
        labs(
            title = paste0("LD Manhattan Plot: ", snp_id, " (", group_name, ")")
        ) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", linewidth = 1),
            legend.title = element_text(size = 12),
            legend.text  = element_text(size = 10),
            plot.title   = element_text(hjust = 0.5, size = 14),
            plot.margin = margin(5, 5, 5, 5)
        )
    
    ggsave(paste0(name_prefix, "_manhattan.pdf"),
           p, width = 10, height = 4)
    
    cat("Plot saved:", paste0(name_prefix, "_manhattan.png"), "\n")
EOF_R
    
    # ============================================================
    # 运行 R
    # ============================================================
    source ~/anaconda3/etc/profile.d/conda.sh
    conda activate R4
    
    Rscript ${NAME}.locusplot.R ${NAME} ${TOP_SNP} ${group} ${X_MIN} ${X_MAX} ${Y_MAX}
   
cd ${OUTPUT_DIR}
mkdir -p ${TOP_SNP}

cp *.pdf ${TOP_SNP}/
rm *${TOP_SNP}*.*
 
    echo "Finished: ${group}"
    echo "Output: ${OUTPUT_DIR}/${NAME}_manhattan.png"
    
done

echo ""
echo "========================================"
echo "All done!"
echo "Results saved in: ${OUTPUT_DIR}"
echo "Job finished at: $(date)"
echo "========================================"
