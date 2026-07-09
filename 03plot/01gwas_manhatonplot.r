library(data.table)
library(ggplot2)

# ============================
# 读取GWAS结果
# ============================
gwas <- fread("../female.fastGWA")

# 计算 -log10(P)
gwas[, LOGP := -log10(P)]

# 去除chr前缀（如果有）
gwas[, CHR := as.numeric(gsub("chr", "", CHR))]

# ============================
# 计算累积坐标
# ============================
chr_len <- gwas[, .(len = max(POS, na.rm = TRUE)), by = CHR]
setorder(chr_len, CHR)

chr_len[, offset := cumsum(c(0, head(len, -1)))]

gwas <- merge(
  gwas,
  chr_len[, .(CHR, offset)],
  by = "CHR",
  all.x = TRUE
)

setorder(gwas, CHR, POS)

gwas[, BP_cum := POS + offset]

# ============================
# X轴中心位置
# ============================
axis_df <- gwas[, .(center = mean(BP_cum)), by = CHR]
setorder(axis_df, CHR)

# ============================
# 染色体交替颜色
# ============================
axis_df[, chr_group := rep(c("A", "B"), length.out = .N)]

gwas <- merge(
  gwas,
  axis_df[, .(CHR, chr_group)],
  by = "CHR",
  all.x = TRUE
)

# ============================
# 阈值线
# ============================
n_snp <- 352479

thr_gws <- 1 / n_snp
thr_suggest <- 0.05 / n_snp

line_gws <- -log10(thr_gws)
line_suggest <- -log10(thr_suggest)

cat("Genome-wide threshold :", thr_gws, "\n")
cat("-log10 =", line_gws, "\n\n")

cat("Bonferroni threshold :", thr_suggest, "\n")
cat("-log10 =", line_suggest, "\n")

# ============================
# 自定义Y轴范围
# ============================
#ymax <- ceiling(max(gwas$LOGP, na.rm = TRUE)) + 1

# 例如固定到12：
ymax <- 9

# ============================
# 绘图
# ============================
pdf("femalmanhattan_plot.pdf", width = 12, height = 6)

p <- ggplot(gwas, aes(x = BP_cum, y = LOGP)) +

  geom_point(
    aes(color = chr_group),
    size = 1.5,
    alpha = 0.7,
    shape = 16 
  ) +

  scale_color_manual(
    values = c(
      "A" = "#EA5B5A",
      "B" = "#F3B2AF"
    )
  ) +

  # Genome-wide线 (红色)
  geom_hline(
    yintercept = line_gws,
    colour = "#CBC9E2",
    linetype = "dashed",
    linewidth = 0.8
  ) +

  # Bonferroni线 (蓝色)
  geom_hline(
    yintercept = line_suggest,
    colour = "#542788",
    linetype = "dashed",
    linewidth = 0.8
  ) +

  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = expansion(mult = 0.01)
  ) +

  scale_y_continuous(
    limits = c(0, ymax),
    breaks = seq(0, ymax, by = 2),
    expand = expansion(mult = c(0, 0))
  ) +

  labs(
    x = "Chromosome",
    y = expression(-log[10](P)),
    title = "GWAS Manhattan Plot"
  ) +

  theme_classic() +

  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

print(p)

dev.off()

cat("Manhattan图已保存: manhattan_plot.pdf\n")

ggsave("femalemanhattan_plot.png", p, width = 12, height = 6, dpi = 300)
