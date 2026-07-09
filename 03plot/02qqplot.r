library(data.table)
library(CMplot)

# 读取结果
gwas <- fread("../female.fastGWA")

# 计算lambda
calc_lambda <- function(p){
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  z <- abs(qnorm(p / 2))
  round(median(z^2, na.rm=TRUE) / 0.454, 3)
}

lambda <- calc_lambda(gwas$P)

# 准备CMplot输入
plot_dat <- data.frame(
  SNP = gwas$SNP,
  Chromosome = 1,
  Position = 1:nrow(gwas),
  P = gwas$P
)

# 过滤有效P值
plot_dat <- plot_dat[
  is.finite(plot_dat$P) &
  !is.na(plot_dat$P) &
  plot_dat$P > 0 &
  plot_dat$P <= 1,
]

# QQ图
CMplot(
  plot_dat,
  plot.type = "q",
  conf.int.col = NULL,
  box = TRUE,
  file = "pdf",
  dpi = 500,
  file.output = TRUE,
  file.name = "female",
  main = paste0("QQ plot   λGC = ", lambda),
  verbose = FALSE
)

cat("QQ图已保存: QQ_plot.pdf\n")
cat("λ =", lambda, "\n")
01plot.r (END)
