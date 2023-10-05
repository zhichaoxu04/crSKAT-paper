rm(list=ls())
library(tidyverse)
library(cmprsk)
library(ggplot2)
library(nleqslv)
library(CompQuadForm)
library(ICSKAT)
library(bindata)
library(data.table)
library(purrr)
library(rlist)
library(foreach)
require(doSNOW)
library(parallel)
library(iterators)
library(snow)


# source("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/RS/Rscript/Function_092722.R")
# source("S:/Rotation/RS/Rscript/Function_092722.R")
source("/rsrch3/scratch/biostatistics/zxu7/Rotation/RS/Rscript/Function_100922.R")
paste0("---- Packages Loaded ----")
# source("S:/Rotation/RS/Rscript/Function_070522.R")


# input
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])  # Cores
aID2 <- as.numeric(args[2]) # B
aID3 <- as.numeric(args[3]) # Sample Size
aID4 <- as.numeric(args[4]) # Number of Q
aID5 <- as.numeric(args[5]) # Random seed
aID6 <- as.numeric(args[6]) # Effect size par 2
aID7 <- as.numeric(args[7]) # Effect size par 2


numCores <- aID
startB <- 1
endB <- aID2

# -- Parameters for data generation
alpha_1 <- -0.058
alpha_2 <- -0.035
beta_1 <- 0.03
beta_2 <- log(1 - exp(beta_1 / alpha_1)) * alpha_2


SampleSize <- aID3
numQ <- aID4 # 50
rho_test <- 0.1 # Correlation for generating gMat
alpha <- 0.05

paste0("alpha1=", alpha_1, ",alpha2=", alpha_2, ",beta1=", beta_1, ",beta2=", beta_2,
       ",n=", SampleSize, ",Q=", numQ, ",rho=", rho_test)

outputDir <- paste0("/rsrch3/scratch/biostatistics/zxu7/Rotation/RS/Rscript/Sim_091223/Result/", sep="") #### Modify Every Time!

min <- aID6
max <- aID7



# ----------------------
numCausal_test <- seq(2, 14, 1)
pwr_result <- matrix(NA, nrow = length(numCausal_test), ncol = 10)
set.seed(aID5)
effPool <- runif(20, min=min, max=max)

if(aID5 == -1){
  effPool <- rep(c(1,-1), 10) * runif(20, min=min, max=max)
}

for(i in 1:length(numCausal_test)){
  numTemp <- numCausal_test[i]
  effTemp <- effPool[1:numTemp]
  
  tem <- cpICSCAT_pwr(alpha1 = alpha_1, alpha2 = alpha_2, beta1 = beta_1, beta2 = beta_2, 
                      numCausal = numTemp, effectSizes = effTemp,
                      n = SampleSize, q = numQ, rho = rho_test, 
                      init_1 = NULL, quants = NULL, numcores = numCores, B1 = startB, B = endB)
  
  
  tem <- tem %>% mutate_all(as.numeric) %>% slice(1:200)
  ss <- tem %>% filter(complete.cases(.)) %>% nrow(.)
  
  cat("causal SNPs =", numTemp, "\n")
  cat("effect sizes =", effTemp, "\n")
  cat("min effect sizes =", min, " max effect sizes =", max, "\n")
  cat("valid obs =", ss, "\n")
  
  n1 <- sum(!is.na(tem$crICSKAT_pval))
  n2 <- sum(!is.na(tem$crBurden_pval))
  n3 <- sum(!is.na(tem$ICSKAT_pval))
  n4 <- sum(!is.na(tem$ICBurden_pval))
  
  pwr1 <- sum(tem$crICSKAT_pval <= alpha) / n1
  pwr2 <- sum(tem$crBurden_pval <= alpha) / n2
  pwr3 <- sum(tem$ICSKAT_pval <= alpha) / n3
  pwr4 <- sum(tem$ICBurden_pval <= alpha) / n4
  
  # cat("crICSKAT pwr =", pwr1, "\n")
  # cat("crBurden pwr =", pwr2, "\n")
  # cat("ICSKAT pwr =", pwr3, "\n")
  # cat("Burden pwr =", pwr4, "\n")
  paste("Finished at", Sys.time())
  cat("\n")
  
  pwr_result[i, ] <- c(numTemp, pwr1,n1, pwr2,n2, pwr3, n3, pwr4, n4, ss)
  
}

cat(paste0("Finished at ", Sys.time()), "\n")
print(pwr_result)
colnames(pwr_result) <- c("Q", "crSKAT_pwr", "crSKAT_n", "crBurden_pwr", "crBurden_n",
                          "ICSKAT_pwr", "ICSKAT_n", "ICBurden_pwr", "ICBurden_n", "N")
cat("\n")

# write results
setwd(outputDir)
outRoot <- "P2_"
write.table(pwr_result, paste0(outRoot, aID, "_", aID2, "_", aID3, "_", 
                               aID4, "_", aID5, "_", aID6, "_", aID7, ".txt"), append=F, quote=F, row.names=F, col.names=T) 

cat("Job Done! ")
