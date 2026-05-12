library(data.table)
library(ggplot2)

#========================================================
# input
#========================================================
pat_file <- "/home/Mzhou/02.F2/z02data/14regnie/02_step2/paternal/AntebrachialL/step2_pat_AntebrachialL.regenie"

mat_file <- "/home/Mzhou/02.F2/z02data/14regnie/02_step2/maternal/AntebrachialL/step2_mat_AntebrachialL.regenie"

outdir <- "/home/Mzhou/02.F2/z02data/14regnie/07_datat_muscript/01GWAS_compare"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#========================================================
# significance threshold
#========================================================
sig_cutoff <- -log10(2.35e-9)

#========================================================
# read data
#========================================================
pat <- fread(
  pat_file,
  select = c("ID", "BETA", "LOG10P")
)

mat <- fread(
  mat_file,
  select = c("ID", "BETA", "LOG10P")
)

setnames(
  pat,
  c("BETA", "LOG10P"),
  c("PAT_BETA", "PAT_LOG10P")
)

setnames(
  mat,
  c("BETA", "LOG10P"),
  c("MAT_BETA", "MAT_LOG10P")
)

#========================================================
# merge
#========================================================
 rm(pat, mat)
gc()

#========================================================
# keep SNPs significant in either PAT or MAT
#========================================================
sig_dat <- dat[
  PAT_LOG10P >= sig_cutoff |
    MAT_LOG10P >= sig_cutoff
]

cat("Total significant SNPs:",
    nrow(sig_dat), "\n")

#========================================================
# summary
#========================================================

# PAT only
pat_only <- sum(
  sig_dat$PAT_LOG10P >= sig_cutoff &
    sig_dat$MAT_LOG10P < sig_cutoff
)

# MAT only
mat_only <- sum(
  sig_dat$PAT_LOG10P < sig_cutoff &
    sig_dat$MAT_LOG10P >= sig_cutoff
)

# both
both_sig <- sum(
  sig_dat$PAT_LOG10P >= sig_cutoff &
    sig_dat$MAT_LOG10P >= sig_cutoff
)

# larger effects
pat_larger <- sum(
  abs(sig_dat$PAT_BETA) >
    abs(sig_dat$MAT_BETA)
)

mat_larger <- sum(
  abs(sig_dat$PAT_BETA) <
    abs(sig_dat$MAT_BETA)
)

# opposite direction
opposite_dir <- sum(
:
dat <- merge(pat, mat, by = "ID")

   sign(sig_dat$PAT_BETA) !=
    sign(sig_dat$MAT_BETA)
)

same_dir <- sum(
  sign(sig_dat$PAT_BETA) ==
    sign(sig_dat$MAT_BETA)
)

#========================================================
# correlations
#========================================================
beta_cor <- cor(
  sig_dat$PAT_BETA,
  sig_dat$MAT_BETA,
  use = "complete.obs"
)

p_cor <- cor(
  sig_dat$PAT_LOG10P,
  sig_dat$MAT_LOG10P,
  use = "complete.obs"
)

#========================================================
# paired t-test
#========================================================
ttest <- t.test(
  abs(sig_dat$PAT_BETA),
  abs(sig_dat$MAT_BETA),
  paired = TRUE
)

#========================================================
# write summary
#========================================================
summary_txt <- c(
  paste0("Total significant SNPs: ",
         nrow(sig_dat)),
  "",
  paste0("PAT only significant: ",
         pat_only),
  paste0("MAT only significant: ",
         mat_only),
  paste0("Both significant: ",
         both_sig),
  "",
  paste0("PAT effect larger: ",
         pat_larger),
   paste0("MAT effect larger: ",
         mat_larger),
  "",
  paste0("Same direction SNPs: ",
         same_dir),
  paste0("Opposite direction SNPs: ",
         opposite_dir),
  "",
  paste0("BETA correlation: ",
         round(beta_cor, 4)),
  paste0("LOG10P correlation: ",
         round(p_cor, 4)),
  "",
  paste0("Paired t-test P-value: ",
         signif(ttest$p.value, 4))
)

writeLines(
  summary_txt,
  con = file.path(
    outdir,
    "GWAS_compare_summary_sigSNP.txt"
  )
)

#========================================================
# save significant SNPs
#========================================================
fwrite(
  sig_dat,
  file = file.path(
    outdir,
    "PAT_MAT_sigSNPs.txt.gz"
  ),
  sep = "\t"
)

#========================================================
# plot 1 : LOG10P scatter
#========================================================
p1 <- ggplot(
  sig_dat,
  aes(PAT_LOG10P,
      MAT_LOG10P)
) +
  geom_point(
    alpha = 0.5,
    size = 1
  ) +
:
