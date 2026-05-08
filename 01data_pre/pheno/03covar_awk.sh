# 读取 PCA
pca <- read.table("/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/01data/PCA/EE_growth_pca.eigenvec", header = FALSE)

# 修改列名
colnames(pca)[1:2] <- c("FID", "IID")

# 把FID改成IID
pca$FID <- pca$IID

# 给PC命名
colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

# <U+2705> 只保留前5个PC
pca_sub <- pca[, c("FID", "IID", paste0("PC", 1:5))]

# 读取协变量
covar <- read.table("covar_growth_ldak.txt", header = TRUE)

# 合并
dat <- merge(covar, pca_sub, by=c("FID","IID"))

# ----------------------------
# 拆分协变量
# ----------------------------

# 1️⃣ 分类变量（factor covariates）
covar_factor <- dat[, c("FID", "IID", "Sex", "farm")]

# 2️⃣ 连续变量（quantitative covariates）
covar_quant <- dat[, c("FID", "IID", "age", paste0("PC", 1:5))]

# ----------------------------
# （可选但强烈建议）标准化 age
# ----------------------------
#covar_quant$age <- scale(covar_quant$age)

# ----------------------------
# 输出文件
# ----------------------------

write.table(covar_factor,
            "covar_factor.txt",
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE)

write.table(covar_quant,
            "covar_quant.txt",
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE)
