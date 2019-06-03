setwd(choose.dir(default = "", caption = "Select folder labeled: 'Source Code'"))
library(plyr)
library(plotrix)
library(gplots)
library(illuminaio)

df <- read.table("Data/Outputs/Final All 3.txt")
df <- subset(df, df$Sig == "Sig" & df$Class == "Hypervariable")


####################
# Read Prelim Data #
####################

test <- read.csv("Data/Full Brain Annotations.csv")
colnames(test) <- as.vector(unlist(test[2,]))
test <- test[,-1]
test <- test[c(10,11,12,13),]
test <- as.data.frame(t(test))
colnames(test) <- c("Tissue", "Sex", "Age", "PMI")
samples <- rownames(test)

test <- apply(test, 2, function(dat) {
  dat <- as.character(dat)
  dat <- gsub(".*: ", "", dat)
})

test <- as.data.frame(test)
rownames(test) <- samples
test$Age <- as.numeric(test$Age)
test$PMI <- as.numeric(test$PMI)
test$Sex <- as.factor(test$Sex)

age.cere <- subset(test, test$Tissue == "cerebellum")
cerebellum <- subset(df, df$Tissue == "Cerebellum")
cere.hyper <- read.table("Data/Outputs/Cerebellum Unimodal.txt")
cere.hyper <- cere.hyper[na.omit(match(cerebellum$Probes,rownames(cere.hyper))),]

age.frontal <- subset(test, test$Tissue == "frontal cortex")
frontal <- subset(df, df$Tissue == "Frontal Cortex")
frontal.hyper <- read.table("Data/Outputs/Frontal Unimodal.txt")
frontal.hyper <- frontal.hyper[na.omit(match(frontal$Probes, rownames(frontal.hyper))),]

genelist <- readBGX("Data/HumanHT-12_V3_0_R3_11283641_A.bgx")
genelist <- data.frame(genelist$probes)
genelist <- data.frame(Gene.Symbol=genelist$ILMN_Gene,
                       Probes=genelist$Probe_Id)

LR_functions <- function(dat){
  dat <- t(dat)
  dat <- merge(test, dat, by="row.names")
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  options(contrasts = c("contr.sum","contr.poly"))
  test.DF <- NULL
  for (i in 5:ncol(dat)) {
    gene <- colnames(dat)[i]
    test.lm <- summary(lm(dat[,i]~dat$PMI+dat$Sex+dat$Age))
    test.lm <- as.data.frame(cbind(gene,
                                   test.lm$coefficients['dat$PMI',4],
                                   sign(test.lm$coefficients['dat$PMI',1]),
                                   test.lm$coefficients['dat$Sex1',4],
                                   sign(test.lm$coefficients['dat$Sex1',1]),
                                   test.lm$coefficients['dat$Age',4],
                                   sign(test.lm$coefficients['dat$Age',1])
    ))
    test.DF <- rbind(test.DF, test.lm)}
  colnames(test.DF) <- c("Probe", "PMI.Pr", "PMI.Sign", "Sex.Pr", "Sex.Sign", "Age.Pr", "Age.Sign")
  
  test.DF$PMI.Pr <- as.numeric(as.character(test.DF$PMI.Pr))
  test.DF$PMI.Pr <- p.adjust(test.DF[order(test.DF$PMI.Pr),"PMI.Pr"], method="BH", n=nrow(test.DF))
  
  test.DF$Sex.Pr <- as.numeric(as.character(test.DF$Sex.Pr))
  test.DF$Sex.Pr <- p.adjust(test.DF[order(test.DF$Sex.Pr),"Sex.Pr"], method="BH", n=nrow(test.DF))
  
  test.DF$Age.Pr <- as.numeric(as.character(test.DF$Age.Pr))
  test.DF$Age.Pr <- p.adjust(test.DF[order(test.DF$Age.Pr),"Age.Pr"], method="BH", n=nrow(test.DF))
  
  
  LR.df <- as.data.frame(cbind(nrow(subset(test.DF, test.DF$Sex.Pr < 0.01 & test.DF$Sex.Sign == 1)),
                               nrow(subset(test.DF, test.DF$Sex.Pr < 0.01 & test.DF$Sex.Sign == -1)),
                               nrow(subset(test.DF, test.DF$Sex.Pr < 0.01)),
                               nrow(subset(test.DF, test.DF$PMI.Pr < 0.01 & test.DF$PMI.Sign == 1)),
                               nrow(subset(test.DF, test.DF$PMI.Pr < 0.01 & test.DF$PMI.Sign == -1)),
                               nrow(subset(test.DF, test.DF$PMI.Pr < 0.01)),
                               nrow(subset(test.DF, test.DF$Age.Pr < 0.01 & test.DF$Age.Sign == 1)),
                               nrow(subset(test.DF, test.DF$Age.Pr < 0.01 & test.DF$Age.Sign == -1)),
                               nrow(subset(test.DF, test.DF$Age.Pr < 0.01))))
  colnames(LR.df) <- c( "Sex Up", "Sex Down", "Sex Total",
                        "PMI Up", "PMI Down", "PMI Total",
                        "Age Up", 'Age Down', "Age Total")
  LR.df
}

table_df <- rbind(LR_functions(cere.hyper),
                  LR_functions(frontal.hyper))

table_df