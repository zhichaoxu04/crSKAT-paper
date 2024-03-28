rm(list=ls())
library(tidyverse)
library(cmprsk)
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
library(tidyr)

source("/YourPathTo/FunctionForSimulations.R")

# input
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])  # Cores
aID2 <- as.numeric(args[2]) # B
aID3 <- as.numeric(args[3]) # Sample Size
aID4 <- as.numeric(args[4]) # Number of Q
aID5 <- as.numeric(args[5]) # Random seed
aID6 <- as.numeric(args[6]) # Rho
aID7 <- as.numeric(args[7]) # Scenarios
aID8 <- as.numeric(args[8]) # WVIC or not
aID9 <- as.numeric(args[9]) # realG or not

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
rho_test <- aID6 # Correlation for generating gMat
alpha <- 0.05

if(aID8==1){
  WVIC_ind=TRUE
  aID8_name="WVIC"
}else{
  WVIC_ind=FALSE
  aID8_name="NOWVIC"
}

if(aID9==1){
  load("/YourPathToRealGeno/RealGeno_raw.rData")

  ## Here we input two times samples for NAN selection
  gMat <- cleanG %>%
    dplyr::slice(1:(SampleSize)) %>%
    dplyr::select(all_of(1:(1+numQ-1))) %>%
    as.matrix(.)
  rm(cleanG)

  aID9_name="RealG"
  
}else{
  gMat=NULL
  aID9_name="SimG"
}

paste0("alpha1=", alpha_1, ",alpha2=", alpha_2, ",beta1=", beta_1, ",beta2=", beta_2,
       ",n=", SampleSize, ",Q=", numQ, ", WVIC=", aID8_name, ", gMat=", aID9_name)

# ----------------------
numCausal_test <- seq(2, 15, 1)
pwr_result <- matrix(NA, nrow = length(numCausal_test), ncol = 12)

if(aID5==1){
  # Define rho_test
  maxMAF_test <- c(0.03, 0.05, 0.15)
  
  # Define unique effect sizes (as it seems you want a 20-length vector for each)
  unique_eff_sizes <- c(0.01, 0.02, 0.03, 0.04, 0.05,
                        0.05, 0.06, 0.07, 0.08, 0.1)
  
  # Initialize a list to hold the combinations
  expanded_list <- list()
  
  # Counter for indexing the list
  counter <- 1
  
  # Loop through each unique effect size
  for (eff_size in unique_eff_sizes) {
    # Create a 20-length vector for the current effect size
    eff_size_vector <- rep(eff_size, 20)
    
    # Loop through each rho value
    for (maxMAF in maxMAF_test) {
      # Store the combination in the list
      expanded_list[[counter]] <- list(eff_size_vector = eff_size_vector, maxMAF = maxMAF)
      counter <- counter + 1
    }
  }
  
  effPool <- expanded_list[[aID7]]$eff_size_vector
  effOut <- as.numeric(effPool[1])
  maxMAF <- expanded_list[[aID7]]$maxMAF
  
}else if(aID5==2){
  set.seed(aID5)
  # Define rho_test
  maxMAF_test <- c(0.03, 0.05, 0.15)
  effPoolSize <- c(0.05, 0.075, 0.1, 0.125, 0.15,
                   1, 1.5, 2, 2.5, 3)
  # Generate effPool_dat with specified ranges of random values
  effPool_dat <- rbind(runif(20, min=-effPoolSize[1], max=effPoolSize[1]), runif(20, min=-effPoolSize[2], max=effPoolSize[2]),
                       runif(20, min=-effPoolSize[3], max=effPoolSize[3]), runif(20, min=-effPoolSize[4], max=effPoolSize[4]), 
                       runif(20, min=-effPoolSize[5], max=effPoolSize[5]), runif(20, min=-effPoolSize[6], max=effPoolSize[6]), 
                       runif(20, min=-effPoolSize[7], max=effPoolSize[7]), runif(20, min=-effPoolSize[8], max=effPoolSize[8]), 
                       runif(20, min=-effPoolSize[9], max=effPoolSize[9]), runif(20, min=-effPoolSize[10], max=effPoolSize[10]))
  
  # Initialize the list to store combinations
  expanded_list <- list()
  
  # Counter to keep track of the list index
  index <- 1
  
  # Loop through each row in effPool_dat
  for (i in 1:nrow(effPool_dat)) {
    # Extract the current row as a vector
    current_effPool <- effPool_dat[i, ]
    
    # Loop through each rho_test value
    for (maxMAF in maxMAF_test) {
      # Add the combination to the list
      expanded_list[[index]] <- list(effPool = current_effPool, maxMAF = maxMAF)
      index <- index + 1
    }
  }
  
  effPool <- expanded_list[[aID7]]$effPool[order(abs(expanded_list[[aID7]]$effPool))]
  effOut <- as.numeric(effPoolSize[((aID7 + 2) %/% 3)])
  maxMAF <- expanded_list[[aID7]]$maxMAF
  
}else if(aID5 == 3){
  maxMAF_test <- c(0.03, 0.05, 0.15)
    effPool_dat <- rbind(rep(c(0.01, -0.01), times=10), rep(c(0.02, -0.02), times=10),
                         rep(c(0.03, -0.03), times=10), rep(c(0.04, -0.04), times=10), rep(c(0.05, -0.05), times=10),
                         rep(c(0.4, -0.4), times=10), rep(c(0.5, -0.5), times=10),
                        rep(c(0.6, -0.6), times=10), rep(c(0.7, -0.7), times=10), rep(c(0.8, -0.8), times=10))
    expanded_list <- list()
    
    index <- 1
    
    # Loop through each row in effPool_dat
    for (i in 1:nrow(effPool_dat)) {
      current_effPool <- effPool_dat[i, ]
      
      for (maxMAF in maxMAF_test) {
        # Add the combination to the list
        expanded_list[[index]] <- list(effPool = current_effPool, maxMAF = maxMAF)
        index <- index + 1
      }
    }
    
    effPool <- expanded_list[[aID7]]$effPool
    effOut <- as.numeric(effPool[1])
    maxMAF <- expanded_list[[aID7]]$maxMAF
}

for(i in 1:length(numCausal_test)){
  numTemp <- numCausal_test[i]
  effTemp <- effPool[1:numTemp]
  
  tem <- cpICSCAT_pwr(alpha1 = alpha_1, alpha2 = alpha_2, beta1 = beta_1, beta2 = beta_2, 
                      numCausal = numTemp, effectSizes = effTemp, gMat=gMat, WVIC=WVIC_ind,
                      n = SampleSize, q = numQ, rho = rho_test, maxMAF=maxMAF,
                      init_1 = NULL, quants = NULL, numcores = numCores, B1 = startB, B = endB)
  
  
  tem <- tem %>% mutate_all(as.numeric) %>% slice(1:200)
  ss <- tem %>% filter(complete.cases(.)) %>% nrow(.)

  cat("causal SNPs =", numTemp, "\n")
  cat("effect sizes =", effTemp, "\n")
  cat("min effect sizes =", min(effTemp), " max effect sizes =", max(effTemp), "\n")
  cat("valid obs =", ss, "\n")

  if(WVIC_ind==TRUE){
    n1 <- sum(!is.na(tem$crICSKAT_pval))
    n2 <- sum(!is.na(tem$crBurden_pval))
    n3 <- sum(!is.na(tem$ICSKAT_pval))
    n4 <- sum(!is.na(tem$ICBurden_pval))
    n5 <- sum(!is.na(tem$pWVIC))

    pwr1 <- sum(tem$crICSKAT_pval <= alpha) / n1
    pwr2 <- sum(tem$crBurden_pval <= alpha) / n2
    pwr3 <- sum(tem$ICSKAT_pval <= alpha) / n3
    pwr4 <- sum(tem$ICBurden_pval <= alpha) / n4
    pwr5 <- sum(tem$pWVIC <= alpha) / n5

  }else{
    n1 <- sum(!is.na(tem$crICSKAT_pval))
    n2 <- sum(!is.na(tem$crBurden_pval))
    n3 <- sum(!is.na(tem$ICSKAT_pval))
    n4 <- sum(!is.na(tem$ICBurden_pval))
    n5 <- NA

    pwr1 <- sum(tem$crICSKAT_pval <= alpha) / n1
    pwr2 <- sum(tem$crBurden_pval <= alpha) / n2
    pwr3 <- sum(tem$ICSKAT_pval <= alpha) / n3
    pwr4 <- sum(tem$ICBurden_pval <= alpha) / n4
    pwr5 <- NA
  }

  paste("Finished at", Sys.time())
  cat("\n")

  pwr_result[i, ] <- c(numTemp, pwr1,n1, pwr2,n2, pwr3, n3, pwr4, n4, pwr5, n5, ss)

}

cat(paste0("Finished at ", Sys.time()), "\n")
colnames(pwr_result) <- c("Q", "crSKAT_pwr", "crSKAT_n", "crBurden_pwr", "crBurden_n",
"ICSKAT_pwr", "ICSKAT_n", "ICBurden_pwr", "ICBurden_n", "WVIC_pwr", "WVIC_n", "N")

pwr_result <- pwr_result %>% 
  as.data.frame() %>% 
  dplyr::mutate(SampleSize=SampleSize, 
                numQ=numQ, 
                rho_test=rho_test, 
                Scen=aID5, 
                EffectSize=effOut, 
                maxMAF=maxMAF)

print(pwr_result)

# write results
outputDir <- paste0("/YourPathToResult/Result/", sep="") #### Modify Every Time!
setwd(outputDir)
outRoot <- "P2_"
write.table(pwr_result, paste0(outRoot, aID, "_", aID2, "_", aID3, "_", 
                               aID4, "_", aID5, "_", aID6, "_", aID7, "_", 
                               aID8_name, "_", aID9_name, ".txt"), append=F, quote=F, row.names=F, col.names=T) 

cat("Job Done! ")
























