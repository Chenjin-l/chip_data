#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# =========================
# 1. 读取数据
# =========================
phe <- read.table(
  "pheno_growth_ldak.rmOutlier.3sd.txt",
  header=TRUE,
  check.names=FALSE
)

cov1 <- read.table(
  "covar_factor.txt",
  header=TRUE
)

cov2 <- read.table(
  "covar_quant.txt",
  header=TRUE
)

# =========================
# 2. 合并协变量
# =========================
cov <- merge(
  cov1[, c("FID","IID","Sex")],
  cov2[, c("FID","IID","age")],
  by=c("FID","IID"),
  all=FALSE
)

# =========================
# 3. 合并表型
# =========================
dat <- merge(phe, cov, by=c("FID","IID"))

# 类型转换
dat$Sex <- factor(dat$Sex)
dat$age <- as.numeric(dat$age)

# =========================
# 4. 定义性状分组
# =========================

# ✔ 需要校正 age 的性状
traits_age <- c(
  "weight",
  "backfat",
  "eye_muscle"
)

# ✔ 不需要 age 的性状（时间定义性状）
traits_no_age <- c(
  "age_100kg",
  "bf_100kg",
  "age_120kg",
  "bf_120kg"
)

traits <- setdiff(colnames(phe), c("FID","IID"))

# =========================
# 5. residual matrix
# =========================
adjusted <- dat

# =========================
# 6. 主循环
# =========================
for (t in traits) {

  cat("processing:", t, "\n")

  y <- as.numeric(dat[[t]])

  use_age <- t %in% traits_age

  # =========================
  # 模型构建
  # =========================
  if (use_age) {
    df <- data.frame(y=y, Sex=dat$Sex, age=dat$age)
    form <- y ~ Sex + age
  } else {
    df <- data.frame(y=y, Sex=dat$Sex)
    form <- y ~ Sex
  }

  keep <- complete.cases(df)
  df2 <- df[keep, ]

  # Sex 必须有变异
  if (nlevels(df2$Sex) < 2) {
    cat("skip (Sex single level):", t, "\n")
    adjusted[[t]] <- NA
    next
  }

  fit <- lm(form, data=df2)

  res <- rep(NA, nrow(dat))
  res[keep] <- residuals(fit)

  adjusted[[t]] <- res

  cat("done:", t, "\n")
}

# =========================
# 7. 输出
# =========================
write.table(
  adjusted,
  file="pheno_adjust_residual.txt",
  quote=FALSE,
  row.names=FALSE,
  sep="\t"
)

cat("ALL DONE\n")
