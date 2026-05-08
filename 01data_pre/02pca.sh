
mkdir -p PCA
plink --bfile EE_QC_growth --pca 20 --out PCA/EE_growth_pca

cat > pca_plot.r << EOF
# 1. PCA读取
# =========================
pca <- read.table("EE_growth_pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID","IID", paste0("PC", 1:(ncol(pca)-2)))

eig <- scan("EE_growth_pca.eigenval")

var_exp <- eig / sum(eig)

# =========================
# 2. Scree plot
# =========================
library(ggplot2)

df <- data.frame(
  PC = 1:length(var_exp),
  Variance = var_exp
)

pdf("PCA_scree_plot.pdf", width=8, height=6)
ggplot(df, aes(x=PC, y=Variance)) +
  geom_line() +
  geom_point(size=2) +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  theme_classic()

dev.off()

# =========================
# 3. PCA scatter plot
# =========================
pdf("PCA_PC1_PC2.pdf", width=8, height=6)

# =========================
# 坐标范围（正确：±5% buffer）
# =========================
x_min <- min(pca\$PC1) * 1.2
x_max <- max(pca\$PC1) * 1.2

y_min <- min(pca\$PC2) * 1.2
y_max <- max(pca\$PC2) * 1.2

# =========================
# PCA plot
# =========================
plot(pca\$PC1, pca\$PC2,
     pch=16,
     cex=0.8,
     col="steelblue",
     xlab=paste0("PC1 (", round(var_exp[1]*100,2), "%)"),
     ylab=paste0("PC2 (", round(var_exp[2]*100,2), "%)"),
     main="PCA: PC1 vs PC2 (Compact Version)",
     xlim=c(x_min, x_max),
     ylim=c(y_min, y_max),
     xaxs="i",
     yaxs="i"
)

dev.off()

EOF

source ~/anaconda3/etc/profile.d/conda.sh                                                                                               
conda activate R4 
Rscript pca_plot.r > plot.log


