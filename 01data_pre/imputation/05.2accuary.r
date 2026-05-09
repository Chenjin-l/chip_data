#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

library(data.table)

# ===============================
# 1. 解析 SNP ID（统一 key = CHR_POS）
# ===============================
parse_snp <- function(ids) {

  ids2 <- gsub("^chr", "", ids)

  res <- lapply(ids2, function(x){

    if(grepl("_", x) & !grepl(":", x)){

      tmp <- strsplit(x, "_")[[1]]

      chr <- tmp[1]
      pos <- tmp[2]
      ref <- tmp[3]
      alt <- tmp[4]

    } else if(grepl(":", x)){

      x2 <- gsub("_.*$", "", x)
      tmp <- strsplit(x2, ":")[[1]]

      chr <- tmp[1]
      pos <- tmp[2]
      ref <- tmp[3]
      alt <- tmp[4]

    } else {
      chr <- pos <- ref <- alt <- NA
    }

    data.frame(
      key = paste0(chr, "_", pos),
      ref = ref,
      alt = alt
    )
  })

  do.call(rbind, res)
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
# 4. 主流程
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
  # SNP overlap
  # ===============================
  common <- intersect(true_id$key, imp_id$key)

  cat("Shared SNPs:", length(common), "\n")

  if(length(common) < 1) stop("No shared SNPs!")

  t_idx <- match(common, true_id$key)
  i_idx <- match(common, imp_id$key)

  ok <- !is.na(t_idx) & !is.na(i_idx)

  t_idx <- t_idx[ok]
  i_idx <- i_idx[ok]
  common <- common[ok]

  true_mat <- true$geno[, t_idx, drop=FALSE]
  imp_mat  <- imp$geno[, i_idx, drop=FALSE]

  # ===============================
  # 🔥 allele matching & filtering
  # ===============================
  t_ref <- true_id$ref[t_idx]
  t_alt <- true_id$alt[t_idx]

  i_ref <- imp_id$ref[i_idx]
  i_alt <- imp_id$alt[i_idx]

  same_flag <- (t_ref == i_ref & t_alt == i_alt)
  flip_flag <- (t_ref == i_alt & t_alt == i_ref)

  keep <- same_flag | flip_flag

  cat("Keep SNPs (allele matched):", sum(keep), "\n")

  # 过滤掉不匹配 SNP（非常关键！！！）
  true_mat <- true_mat[, keep, drop=FALSE]
  imp_mat  <- imp_mat[, keep, drop=FALSE]

  flip_idx <- which(flip_flag[keep])

  # ===============================
  # flip genotype
  # ===============================
  if(length(flip_idx) > 0){
    for(i in flip_idx){
      imp_mat[, i] <- flip_gt(imp_mat[, i])
    }
  }

  # ===============================
  # accuracy
  # ===============================
  acc <- mean(true_mat == imp_mat, na.rm=TRUE)

  snp_acc <- apply(true_mat == imp_mat, 2, mean, na.rm=TRUE)
  sample_acc <- apply(true_mat == imp_mat, 1, mean, na.rm=TRUE)

  r2 <- sapply(1:ncol(true_mat), function(i){
    cor(true_mat[,i], imp_mat[,i], use="complete.obs")^2
  })

  # ===============================
  # summary
  # ===============================
  cat("\n========== SUMMARY ==========\n")
  cat("Truth SNPs total   :", ncol(true$geno), "\n")
  cat("Imputed SNPs total :", ncol(imp$geno), "\n")
  cat("Shared SNPs (raw)  :", length(common), "\n")
  cat("Kept SNPs          :", sum(keep), "\n")
  cat("Shared samples     :", nrow(true$geno), "\n")
  cat("Overall accuracy   :", round(acc, 4), "\n")
  cat("Mean sample acc    :", round(mean(sample_acc), 4), "\n")
  cat("Mean SNP acc       :", round(mean(snp_acc), 4), "\n")
  cat("Mean dosage R²     :", round(mean(r2, na.rm=TRUE), 4), "\n")
}

# ===============================
# RUN
# ===============================
run("truth.raw", "imp.raw")
