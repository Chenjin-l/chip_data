  # =========================
# 1. 读取数据
# =========================
pheno <- read.table("../new_growth.txt_en.txt", header=TRUE, stringsAsFactors=FALSE)

keep_id <- read.table(
  "/work/lingcj/03chip_data/02anlyis/01data/01pre_exchange/ID_compare/growth_geno_fullMatch.txt",
  header=FALSE, stringsAsFactors=FALSE
)

colnames(keep_id) <- c("IID")

# =========================
# 2. 过滤个体
# =========================
pheno_sub <- pheno[pheno$id %in% keep_id$IID, ]

cat("保留个体数:", nrow(pheno_sub), "\n")

# =========================
# 3. 多性状表型（核心）
# =========================
pheno_multi <- data.frame(
  FID = pheno_sub$id,
  IID = pheno_sub$id,
  weight = pheno_sub$weight,
  backfat = pheno_sub$backfat,
  eye_muscle = pheno_sub$eye_muscle,
  age_100kg   = pheno_sub$age_100kg,
  bf_100kg    = pheno_sub$bf_100kg,
  age_120kg   = pheno_sub$age_120kg,
  bf_120kg    = pheno_sub$bf_120kg
)

# 排序（非常重要，避免后面软件报错）
pheno_multi <- pheno_multi[order(pheno_multi$IID), ]

write.table(
  pheno_multi,
  "pheno_growth_ldak.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

# =========================
# 4. 协变量
# =========================
covar <- data.frame(
  FID = pheno_sub$id,
  IID = pheno_sub$id,
  Sex = pheno_sub$sex,
  farm = pheno_sub$farm,
  age = pheno_sub$age
)

covar <- covar[order(covar$IID), ]

write.table(
  covar,
  "covar_growth_ldak.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

# =========================
# 5. 一致性检查
# =========================
cat("\n===== 检查 =====\n")

if(all(pheno_multi$IID == covar$IID)){
  cat("✔ IID 完全一致，可直接用于 LDAK/GCTA\n")
} else {
  cat("<U+274C> IID 不一致（需要检查）\n")
}
