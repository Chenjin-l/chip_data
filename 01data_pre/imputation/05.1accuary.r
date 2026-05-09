#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

library(data.table)

# ===============================
# 1. 解析 SNP ID（统一 key = CHR_POS）
# ===============================
parse_snp <- function(ids) {
  ids <- gsub("^chr", "", ids)
  tmp <- strsplit(ids, "_")

  key <- sapply(tmp, function(x) {
    if (length(x) < 2) return(NA)
    paste0(x[1], "_", x[2])
  })

  ref <- sapply(tmp, function(x) if(length(x) >= 3) x[3] else NA)
  alt <- sapply(tmp, function(x) if(length(x) >= 4) x[4] else NA)

  data.frame(key = key, ref = ref, alt = alt)
}

# ===============================
# 2. flip genotype
# ===============================
flip_gt <- function(x) {
  x <- as.numeric(x)
  ifelse(x == 0, 2,
         ifelse(x == 2, 0, 1))
}

# ===============================
# 3. 读取 PLINK raw
# ===============================
read_raw <- function(file) {
  df <- fread(file, data.table = FALSE, check.names = FALSE)

  fam <- df[,1:6]
  geno <- df[,7:ncol(df), drop = FALSE]

  list(fam = fam, geno = as.matrix(geno))
}

# ===============================
# 4. 主函数
# ===============================
run <- function(true_file, imp_file) {

 cat("Reading files...\n")

  true <- read_raw(true_file)
  imp  <- read_raw(imp_file)

  cat("Parsing SNP IDs...\n")

  true_id <- parse_snp(colnames(true$geno))
  imp_id  <- parse_snp(colnames(imp$geno))

  # 去 NA
  true_id <- true_id[!is.na(true_id$key), ]
  imp_id  <- imp_id[!is.na(imp_id$key), ]

  # ===============================
  # overlap SNP
  # ===============================
  common <- intersect(true_id$key, imp_id$key)

  cat("Shared SNPs:", length(common), "\n")

  if(length(common) < 1) stop("No shared SNPs!")

  true_idx <- match(common, true_id$key)
  imp_idx  <- match(common, imp_id$key)

  # 防止 NA
  ok <- !is.na(true_idx) & !is.na(imp_idx)
  true_idx <- true_idx[ok]
  imp_idx  <- imp_idx[ok]
  common   <- common[ok]

  true_mat <- true$geno[, true_idx, drop=FALSE]
  imp_mat  <- imp$geno[, imp_idx, drop=FALSE]

  # ===============================
  # allele flip
  # ===============================
  flip_flag <- (true_id$ref[true_idx] == imp_id$alt[imp_idx] &
                true_id$alt[true_idx] == imp_id$ref[imp_idx])

  for(i in which(flip_flag)) {
    imp_mat[,i] <- flip_gt(imp_mat[,i])
  }

  # ===============================
  # accuracy
  # ===============================
  acc <- mean(true_mat == imp_mat, na.rm=TRUE)

  snp_acc <- apply(true_mat == imp_mat, 2, mean, na.rm=TRUE)
  sample_acc <- apply(true_mat == imp_mat, 1, mean, na.rm=TRUE)

  r2 <- sapply(1:ncol(true_mat), function(i) {
    cor(true_mat[,i], imp_mat[,i], use="complete.obs")^2
  })

  # ===============================
  # output
  # ===============================
cat("\n========== SUMMARY ==========\n")
  cat("Truth SNPs total   :", ncol(true$geno), "\n")
  cat("Imputed SNPs total :", ncol(imp$geno), "\n")
  cat("Shared SNPs        :", length(common), "\n")
  cat("Shared samples     :", nrow(true_mat), "\n")
  cat("Overall accuracy   :", round(acc,4), "\n")
  cat("Mean sample acc    :", round(mean(sample_acc,na.rm=TRUE),4), "\n")
  cat("Mean SNP acc       :", round(mean(snp_acc,na.rm=TRUE),4), "\n")
  cat("Mean dosage R²     :", round(mean(r2,na.rm=TRUE),4), "\n")
  cat("============================\n")
}

# ===============================
# RUN
# ===============================
run("truth.raw", "imp.raw")
