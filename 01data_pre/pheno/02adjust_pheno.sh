rm(list = ls())
options(stringsAsFactors = FALSE)

library(dplyr)

# =========================
# 1. 读取数据
# =========================
phe <- read.table(
  "pheno_growth_ldak.txt",
  header = TRUE,
  check.names = FALSE
)

cat("原始数据维度：", dim(phe), "\n")

# =========================
# 2. 定义异常值处理函数（3SD法）
# =========================
trimPhe <- function(x, threshold = 0.1) {

  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)

  lowlim <- mu - 3 * sdv
  uplim  <- mu + 3 * sdv

  # 异常值判断
  outlier_flag <- x < lowlim | x > uplim
  outlier_count <- sum(outlier_flag, na.rm = TRUE)
  total <- sum(!is.na(x))

  ratio <- outlier_count / total

  # 仅当异常比例 ≤ threshold 才处理
  if (ratio <= threshold) {
    x[outlier_flag] <- NA
  }

  return(list(x = x, outlier = outlier_count, ratio = ratio))
}

# =========================
# 3. 逐性状处理
# =========================
pheno_clean <- phe
report_list <- list()

for (i in 3:ncol(phe)) {  # 前两列FID IID

  res <- trimPhe(phe[, i])

  pheno_clean[, i] <- res$x

  report_list[[colnames(phe)[i]]] <- data.frame(
    Trait = colnames(phe)[i],
    Outliers = res$outlier,
    Outlier_Ratio = res$ratio
  )

  cat("处理完成:", colnames(phe)[i],
      " | outliers:", res$outlier,
      " | ratio:", round(res$ratio, 4), "\n")
}

# =========================
# 4. 合并报告
# =========================
report <- do.call(rbind, report_list)

# =========================
# 5. 输出结果
# =========================
write.table(
  pheno_clean,
  file = "pheno_growth_ldak.rmOutlier.3sd.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  report,
  file = "pheno_outlier_summary_3sd.csv",
  row.names = FALSE
)

cat("\n完成！已输出：\n")
cat("1. pheno_growth_ldak.rmOutlier.3sd.txt（3SD清洗后表型）\n")
cat("2. pheno_outlier_summary_3sd.csv（异常值统计）\n")
