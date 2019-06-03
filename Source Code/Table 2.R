setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))

library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)

dat <- read.table("Data/Outputs/Final All 3.txt")
dat$Class <- revalue(dat$Class, c("NV" = "Non-Variable", "Hypervariable" = "Hyper-Variable",
                                  "Hypovariable" = "Hypo-Variable"))


sig <- as.data.frame(table(dat$Class, dat$Sig, dat$Tissue))
sig <- subset(sig, sig$Var2 == "Sig")
colnames(sig) <- c("EV", "Geneset", "Tissue", "Count")
sig$Geneset <- "After Retest"

totals <- as.data.frame(table(dat$Class, dat$Tissue))
colnames(totals) <- c("EV", "Tissue", "Count")
totals$Geneset <- "Before Retest"
totals <- totals[,c(1,4,2,3)]

ggDF <- as.data.frame(rbind(sig, totals))
ggDF$Geneset <- factor(ggDF$Geneset, levels=c("Before Retest", "After Retest"))
ggDF <- ggDF[-c(3,6,9),]


# Calculate NV
dat[which(dat$Sig != "Sig"),"Class"] <- "Non-Variable"
NV_df <- dat %>% filter(Class == "Non-Variable")
sig2 <- as.data.frame(table(NV_df$Class, NV_df$Tissue))
sig2$Geneset <- "After Retest"
sig2 <- sig2[,c(1,4,2,3)]
colnames(sig2) <- c("EV", "Geneset", "Tissue", "Count")
sig2 <- sig2 %>% filter(Count > 0)


ggDF <- rbind(ggDF, sig2)

probetable <- as.data.frame(acast(ggDF, Tissue*EV~Geneset))
probetable$Class <- rownames(probetable)
probetable %>% separate(Class, c("Tissue", "Probe Set"), sep="_") %>%
  select("Probe Set", "Tissue", "Before Retest", "After Retest") %>%
  mutate("% Probes" = format(`After Retest`/`Before Retest`*100, digits=4)) %>%
  arrange(`Probe Set`)

