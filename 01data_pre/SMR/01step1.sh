
#####生成gwas总结文件

library(data.table)

gwas_file <- "/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/02gwas/01LDAKgwas/01step1/4_backfat/gwas_backfat_step1.step2.assoc"
###Chromosome      Predictor       Basepair        A1      A2      Wald_Stat       Wald_P  Effect  SE      MAF     CallRate        MachR2  SPA_Status
###1       1_170_G_T       170     T       G       -1.5056 1.3216e-01      -7.0268e-02     4.6669e-02      0.055556        1.0000  NA      NOT_USED
###1       1_555_C_G       555     G       C       -1.5056 1.3216e-01      -7.0268e-02     4.6669e-02      0.055556        1.0000  NA      NOT_USED
###1       1_689_ACCTG_A   689     A       ACCTG   0.7618  4.4619e-01      2.8019e-02      3.6781e-02      0.092639        1.0000  NA      NOT_USED

gwas <- fread(gwas_file)

#========================================================
# 1. 清理 SNP 引号
#========================================================
gwas[, SNP := gsub('"', '', Predictor)]

#========================================================
# 2. 只保留SMR需要列
#========================================================
gwas_smr <- gwas[, .(
  SNP  = SNP,
  A1   = A1,
  A2   = A2,
  freq = MAF,
  b    = Effect,
  se   = SE,
  p    = Wald_P,
  n    = NA
)]

#========================================================
# 3. 删除任何 NA 行（SMR非常严格）
#========================================================
gwas_smr <- na.omit(gwas_smr)

#========================================================
# 4. 强制写 tab + 禁止quote
#========================================================
fwrite(
  gwas_smr,
  file = "4_backfat.mygwas.ma",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

#生成eqtl总结文件

library(data.table)

eqtl_file <- "/home/Mzhou/zlingcj/03chip_data/02anlyis/04swim/02gwas/01LDAKgwas/02colo/01predata/01eqtl_data_colo/PigGTEx_all_eqtl_for_coloc.txt.gz"

TISSUE <- "Testis"
TARGET_GENE <- "ENSSSCG00000017278"

dat <- fread(eqtl_file)

#========================================================
# filter
#========================================================
sub <- dat[
  tissue == TISSUE &
  gene_id == TARGET_GENE
]

#========================================================
# split snp
# 5_3247554_C_T
#========================================================
tmp <- tstrsplit(sub$snp, "_")

sub[, A1 := tmp[[3]]]
sub[, A2 := tmp[[4]]]

#========================================================
# build ESD
#========================================================
esd <- data.table(
  Chr  = sub$chr,
  SNP  = sub$snp,
  Bp   = sub$pos,
  A1   = sub$A1,
  A2   = sub$A2,
  Freq = sub$maf,
  Beta = sub$beta,
  se   = sub$se,
  p    = sub$p
)

fwrite(
  esd,
  file = paste0(TARGET_GENE, "_", TISSUE, ".esd"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("ESD done\n")


#========================================================
# probe information
#========================================================

probe_chr <- unique(sub$chr)[1]

probe_bp <- round(median(sub$pos))

flist <- data.table(
  Chr = probe_chr,
  ProbeID = TARGET_GENE,
  GeneticDistance = 0,
  ProbeBp = probe_bp,
  Gene = TARGET_GENE,
  Orientation = "+",
  PathOfEsd = paste0(
    getwd(),
    "/",
    TARGET_GENE,
    "_",
    TISSUE,
    ".esd"
  )
)

fwrite(
  flist,
  file = paste0(TARGET_GENE, "_", TISSUE, ".flist"),
  sep = "\t",
  quote = FALSE
)

cat("FLIST done\n")


