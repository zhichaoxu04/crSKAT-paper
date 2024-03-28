rm(list=ls())
library(foreach)
require(doSNOW)
library(parallel)
library(iterators)
library(snow)
library(tidyverse)
library(cmprsk)
library(tidyverse)
library(ggplot2)
library(nleqslv)
library(CompQuadForm)
library(ICSKAT)
library(data.table)
library(purrr)
library(rlist)
library(ALassoSurvIC)
library(nlme)
library(statmod)
source("FunctionForSimulations.R")

# input
args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])  # Cores
aID2 <- as.numeric(args[2]) # B
aID3 <- as.numeric(args[3]) # Sample Size n
aID4 <- as.numeric(args[4]) # Number of Q
aID5 <- as.numeric(args[5]) # Random seed
aID6 <- as.numeric(args[6]) # Output order

# -------------------------
numcores <- aID1
B <- aID2
B1 <- 1+(aID6-1)*B
B2 <- aID6*B
# -------------------------

# ------------------------
n <- aID3
q <- aID4 # q is number of SNPs

alpha1 <- -0.058
alpha2 <- -0.035
beta1 <- 0.03
beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2


# ---- gMat and xMat ----
set.seed(aID5)
load("RealGeno_raw.rData")

gMat <- cleanG %>%
  slice(sample(1:150000, n, replace=F)) %>%
  select(all_of(c(1:q)))

rm(cleanG)

gMat <- as.matrix(gMat)
gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
xMat <- cbind(rnorm(n), rbinom(n=n, size=1, prob=0.5))

# ------ Simulation
# Start parallel computing
cl <- makeCluster(numcores, type = "SOCK", outfile = "sim_012924.txt")
registerDoSNOW(cl)

# ------------
# Foreach: Parallel computing for loop
# ------------

oper <- foreach(sim_it = B1:B2, .inorder = F, .errorhandling = "pass") %dopar% {
  
  library(nleqslv)
  library(CompQuadForm)
  library(stats)
  library(foreach)
  # library(doParallel)
  
  set.seed(sim_it)
  
  # simulate data     
  tempDat <- genData(seed = sim_it, n=n, alpha1=alpha1, alpha2=alpha2, 
                     beta1=beta1, beta2 = beta2)
  
  quants <- stats::quantile(log(tempDat$tempTime), probs=seq(from=0, to=1, length.out=3))
  ej <- (max(quants) - quants[2]) / (max(quants) - min(quants))
  logH1 <- rep(log(-beta1 / alpha1), length(tempDat$tempTime)) + log(1 - exp(alpha1 * tempDat$tempTime))
  logH2 <- rep(log(-beta2 / alpha2), length(tempDat$tempTime)) + log(1 - exp(alpha2 * tempDat$tempTime))
  regDmat <- cbind(1, log(tempDat$tempTime), 
                   pmax(0, (log(tempDat$tempTime) - quants[2])**3) - ej * pmax(0, (log(tempDat$tempTime) - quants[1])**3) -
                     (1 - ej) * pmax(0, (log(tempDat$tempTime) - quants[3])**3))
  init_1 <- c(as.numeric(summary(lm(logH1 ~ cbind(xMat, regDmat) - 1))$coef[, 1]),
              as.numeric(summary(lm(logH2 ~ cbind(xMat, regDmat, gSummed) - 1))$coef[, 1]))
  
  # observation time
  # simulate missing visit
  madeVisit <- matrix(data=rbinom(n=n*7, size=1, prob=0.9), nrow=n, ncol=7)
  visitTime <- sweep(matrix(data=runif(n=n*7, min=-1, max=1), nrow=n, ncol=7),
                     MARGIN=2, STATS=seq(from=4, to=28, by=4), FUN="+")
  # get all visits for each subject
  allVisits <- madeVisit * visitTime
  # make the interval for each subject - USING creatIntNew!
  allInts <- t(mapply(FUN=createIntNew, obsTimes = data.frame(t(allVisits)), eventTime=tempDat$tempTime))
  leftTimes <- allInts[, 1]
  rightTimes <- allInts[, 2]
  # new indicator of whether event was observed
  deltaVecSimple <- ifelse(rightTimes > tempDat$tempTime, 1, 0)
  deltaVec <- deltaVecSimple * tempDat$tempType
  
  # make spline terms
  # for interval censoring, our earliest visit time is 4 and our latest is 28, windows of size 1
  
  dmats <- makeICdmat(xMat=xMat, lt = leftTimes, rt = rightTimes, obsInd = deltaVecSimple, 
                      quant_r=quants, nKnots=1) 
  leftDmat <- dmats$left_dmat
  rightDmat <- dmats$right_dmat
  
  # solve score equations under null with no random effect
  solPartial <- nleqslv(x = init_1, 
                        fn=scoreEqSpline, leftDmat = leftDmat,  rightDmat = rightDmat, leftTimes = leftTimes,
                        deltaVec = deltaVec, gSummed = gSummed, gMat=NULL, estG = FALSE)

  if (solPartial$termcd > 2) {
    next
  }
  
  if(is.numeric(solPartial$x) == FALSE){
    next
  }
  
  # Ugamma with null coefficients
  forUgamma <- scoreEqSpline(x=solPartial$x, leftDmat = leftDmat, 
                             rightDmat = rightDmat, leftTimes = leftTimes, deltaVec = deltaVec,
                             gSummed = gSummed, gMat = gMat, estG = FALSE)
  
  # information
  # partial information, under null
  iMatPartial <- calcInfo(leftDmat = leftDmat, rightDmat = rightDmat, leftTimes = leftTimes, 
                          theta1 = solPartial$x[1:5], theta2 = solPartial$x[6:11], deltaVec = deltaVec,
                          gSummed=gSummed, gMat = NULL)
  
  # Igt
  forIgt <- calcInfo(leftDmat = leftDmat, rightDmat = rightDmat, leftTimes = leftTimes, 
                     theta1 = c(solPartial$x[1:5], rep(0, q)), theta2 = solPartial$x[6:11], deltaVec = deltaVec,
                     gSummed = gSummed, gMat = gMat)
  
  # try SKAT
  skatQ <- t(forUgamma[6:(5+q)]) %*% forUgamma[6:(5+q)]
  burdenQ <- (sum(forUgamma[6:(5+q)]))^2
  Itt <- -iMatPartial
  Igg <- -forIgt[6:(5+q), 6:(5+q)]
  Igt <- -forIgt[c(1:5, (6+q):nrow(forIgt)), 6:(5+q)]
  sig_mat <- Igg - t(Igt) %*% solve(Itt) %*% (Igt)
  lambdaQ <- eigen(sig_mat)$values
  p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq
  B_burden= burdenQ / sum(sig_mat)
  p_burden= 1 - stats::pchisq(B_burden, df = 1)
  
  
  # ------ ICSKAT ------
  # obs_ind <- deltaVecSimple               # Treat cause2 as event
  obs_ind <- ifelse(deltaVec == 1, 1, 0)    # Treat cause2 as censoring
  tpos_ind <- ifelse(leftTimes == 0, 0, 1)
  init_2 <- init_1[1:5]
  nullFit <- ICSKAT::ICSKAT_fit_null(init_beta = init_2, left_dmat = leftDmat, right_dmat = rightDmat,
                             obs_ind = obs_ind, tpos_ind = tpos_ind, lt = leftTimes, rt = rightTimes)
  
  if (is.na(as.numeric(nullFit$beta_fit[1]))){
    pICSKAT <- NA
    pICburden <- NA
    skatQ2 <- NA
    burdenQ2 <- NA
    err <- nullFit$errMsg
  } else{
    # perform the ICSKAT and Burden tests
    icskatOut <- ICSKAT::ICskat(left_dmat = leftDmat, right_dmat = rightDmat, lt = leftTimes, rt = rightTimes,
                        obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, 
                        null_beta = as.numeric(nullFit$beta_fit), Itt = nullFit$Itt)
    pICSKAT <- as.numeric(icskatOut$p_SKAT)
    pICburden <- as.numeric(icskatOut$p_burden)
    skatQ2 <- as.numeric(icskatOut$skatQ)
    burdenQ2 <- as.numeric(icskatOut$burdenQ)
    err <- icskatOut$errMsg
  }
  

  # ----- WVIC ------
  ltWVIC <- ifelse(obs_ind == 1, leftTimes, rightTimes)
  rtWVIC <- ifelse(obs_ind == 1, rightTimes, Inf)
  bound_times <- cbind(ltWVIC, rtWVIC)
  test_results <- WVIC_test(X=as.matrix(gMat),
                            Z=NULL,
                            bound_times=bound_times,
                            trunct=NULL,
                            Gsim=Gsim, Fsim=NULL,
                            r=0, lim=5*10^5, acc=10^(-5),
                            covmat=FALSE)
  pWVIC <- test_results$pval
  rm(Gsim)
  pWVIC <- NA
  
  

  result1 <- as.character(c(skatQ, burdenQ, p_SKAT, p_burden, B_burden, pWVIC,
                            # solPartial$message, solPartial$termcd, init_1[3], solPartial$x, quants,
                            pICSKAT, pICburden, skatQ2, burdenQ2, err, sim_it))
 
  return(result1)
}

# Stop parallel computing
stopCluster(cl)

# ------ Process Result
filtered <- Filter(function(x) length(x) == 12, oper)
skatDF <- as.data.frame(do.call(rbind, filtered))
colnames(skatDF) <- c("skatQ", "burdenQ", "pcrSKAT", "pcrBurden", "B_burden", "pWVIC", 
                    # "message", "termcd", "InitialTheta3", paste0("theta", rep(1:11)), paste0("quant", rep(1:3)),
                    "pICSKAT", "pICburden", "skatQ2", "burdenQ2", "err", "sim_it")

save(list = c("skatDF"),
     file = paste0("/YourPathTo/T1Sim", "_",aID3, "_", aID4, "_", aID6, ".rData"))

