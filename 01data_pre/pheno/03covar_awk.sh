#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

library(data.table)

# ===============================
# 1. 读取 PCA
# ===============================
pca <- fread("/work/lingcj/03chip_data/02anlyis/01data/03growth/01bfile/02vcf/02_1k_pigimputation/06merge/bfile/PCA/PCA.eigenvec",
             header = FALSE)

colnames(pca)[1:2] <- c("FID", "IID")
pca$FID <- pca$IID

# PC命名
pc_num <- ncol(pca) - 2
colnames(pca)[3:ncol(pca)] <- paste0("PC", 1:pc_num)

# 只保留前10个PC
pca_sub <- pca[, c("FID", "IID", paste0("PC", 1:min(10, pc_num))), with = FALSE]


# ===============================
# 2. 读取协变量
# ===============================
covar <- fread("covar_growth_ldak.txt", header = TRUE)

# 合并
dat <- merge(covar, pca_sub, by = c("FID", "IID"))


# ===============================
# 3. 自动去除常数列函数
# ===============================
remove_constant_cols <- function(df) {
  keep <- sapply(df, function(x) {
    length(unique(na.omit(x))) > 1
  })
  df[, keep, with = FALSE]
}


# ===============================
# 4. 分类协变量（factor）
# ===============================
covar_factor <- dat[, .(FID, IID, Sex, farm)]

# 去掉单一变量
covar_factor <- remove_constant_cols(covar_factor)


# ===============================
# 5. 连续协变量（quantitative）
# ===============================
covar_quant <- dat[, c("FID", "IID", "age", paste0("PC", 1:min(10, pc_num))), with = FALSE]

# 去掉单一变量
covar_quant <- remove_constant_cols(covar_quant)


# ===============================
# 6. 输出检查（关键！）
# ===============================
cat("\n===== Covariate Summary =====\n")

cat("\nFactor covariates:\n")
print(sapply(covar_factor[, -c("FID","IID"), with=FALSE], function(x) length(unique(x))))

cat("\nQuantitative covariates:\n")
print(sapply(covar_quant[, -c("FID","IID"), with=FALSE], function(x) length(unique(x))))

cat("\nFinal factor cols:\n")
print(colnames(covar_factor))

cat("\nFinal quant cols:\n")
print(colnames(covar_quant))

cat("============================\n")


# ===============================
# 7. 输出文件
# ===============================
fwrite(covar_factor, "covar_factor.txt")
fwrite(covar_quant,  "covar_quant.txt")
