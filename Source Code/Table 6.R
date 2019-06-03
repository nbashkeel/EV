setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))


dat <- read.table("Data/Outputs/Final Meth Summary.txt")

check_enrich <- function(dat) {
  df <- table(dat$Class, dat$Methylation.Cluster)
  resi <- list("Residuals" = round(chisq.test(df)$stdres,2),
               "P-Value" = chisq.test(df)$p.value)
  resi
}

meth_enrichments <- by(dat, INDICES = dat$Tissue, check_enrich)
meth_enrichments
