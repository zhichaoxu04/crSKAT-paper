# Apply interval-censored SKAT test to fractures (6151) outcome from questionnaire 
# in the UK Biobank. 
# All chromosomes.
# Weighted Genotype + ICSKAT package functions
library(tidyverse)
library(tidyr)
library(SeqArray)
library(data.table)
library(dplyr)
library(magrittr)
library(Rcpp)
library(GWASTools)
library(SNPRelate)
library(survival)
library(pryr)
library(ICSKAT)
library(nleqslv)
# source("S:/Rotation/RS/Rscript/Function_100922.R")
source("/rsrch3/scratch/biostatistics/zxu7/Rotation/RS/Rscript/Function_100922.R")


# input
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
aID2 <- as.numeric(args[2])
aID3 <- as.numeric(args[3])

# ---------- Some Arguments ------------------
#max_gene_size <- 10^10
max_gene_size <- 10^10
buffer <- 5000
genes_per_run <- aID2
method <- c("Broyden", "Newton")[aID3]
outRoot <- "cricskatFractures_"
# outputDir <- "S:/Rotation/RS/RDA/Out/"
outputDir <- paste0("/rsrch3/scratch/biostatistics/zxu7/Rotation/RS/RDA/Out/091323/Result/", aID5, sep="") #### Modify Every Time!

# Load Data
# workpath <- "S:/Rotation/RS/RDA/"
workpath <- "/rsrch3/scratch/biostatistics/zxu7/Rotation/RS/RDA/"

# Here load the outcomeDat as fullDat
load(paste0(workpath, "fullDat_112922.rData"))
hg19List <- fread(paste0(workpath, "hg19_gene_definitions.txt"))

# ----- Process outcome, gene list and covariates
# longer list
load(file=paste0(workpath, "ensembl_refgene_hg19_20180109.rda"))

all_genes <- ensembl_refgene_hg19_20180109 %>%
  filter(Notes == 0 & Chr <= 22) %>%
  mutate(Diff = txEnd - txStart) %>% 
  arrange(Diff) %>%
  filter(Diff < 10000 | HGNC_name %in% hg19List$V1)

outcomeDat <- fullDat %>%
  # remove the noAns!
  filter(noAns == 0) %>%
  select(eid, Lage1, Rage1, Ind) %>%
  set_colnames(c("ID", "leftDays", "rightDays", "deltaVec")) %>%
  filter(!is.na(ID) & !is.na(leftDays) & !is.na(rightDays) & !is.na(deltaVec)) %>%
  mutate(leftTime = case_when(deltaVec == 0 ~ as.numeric(0), TRUE ~ leftDays),
         rightTime = case_when(deltaVec == 0 ~ as.numeric(leftDays), TRUE ~ rightDays),
         leftTime1 = leftDays,
         rightTime1 = case_when(rightDays == 99999 ~ as.numeric(999), TRUE ~ rightDays),
         deltaVecSimple = case_when(rightTime1 < 999 ~ as.numeric(1), TRUE ~  as.numeric(0)), # obs_ind (0 = right censored)
         leftCensoredInd = case_when(leftTime1 > 0 ~ as.numeric(1), TRUE ~ as.numeric(0))) # tpos_ind (0 = left censored)

ageSexDat <- fread(paste0(workpath, "ukbAgeSex.txt"), header=T, data.table=F) %>%
  set_colnames(c("eid", "Age", "Age2", "Age3", "Age4", "Sex")) %>% 
  # select(eid, Age, Sex) %>% 
  select(eid, Sex) %>% 
  # filter(!is.na(eid) & !is.na(Age) & !is.na(Sex)) %>%
  filter(!is.na(eid) & !is.na(Sex))

heightWeightDat <- fread(paste0(workpath, "UKBheightweight.txt")) %>%
  select(eid, height, weight) %>%
  merge(., ageSexDat, by=c("eid")) %>%
  set_colnames(c("ID", "Height", "Weight", "Sex"))

evecsDat <- fread(paste0(workpath, "evecs_n502524_20200702.txt")) %>%
  set_colnames(c("ID", paste0("evec", 1:40))) %>%
  select(1:11) %>%
  filter(!is.na(ID) & !is.na(evec1))

# merge covariates
allCov <- heightWeightDat %>%
  inner_join(evecsDat, by = "ID") %>% 
  inner_join(outcomeDat, by = "ID") %>%
  mutate(ID = paste0(ID, "_", ID)) 



# slice to genes for this run
startRow <- (aID - 1) * genes_per_run + 1
endRow <- aID * genes_per_run 
gene_info <- all_genes %>%
  slice(startRow:endRow)


# results data frame
resultsDF <- data.frame(gene = gene_info$HGNC_name, chr = gene_info$Chr, start=gene_info$txStart,
                        end=gene_info$txEnd, q=NA, qNB=NA, cleanq=NA, rareq=NA, crskatp=NA, crburdenp=NA, crerr=NA, crInitIter=NA,
                        skatp=NA, burdenp=NA, complex=NA, SKATerr=NA, SKATOp=NA, SKATOr=NA, SKATOdavies=NA, SKATOerr=NA)

# loop through nFilters (make sure this is whole number)
for (gene_it in 1:nrow(gene_info)) {
  cat("Gene ", gene_it, " started:", "\n")
  
  tempChr <- gene_info$Chr[gene_it]
  
  # open GDS file
  setwd("/rsrch3/scratch/biostatistics/zxu7/Rotation/RS/RDA/ukb_gds")
  # setwd("S:/Rotation/RS/RDA/ukb_gds")
  gdsfile <- seqOpen(paste0("ukb_imp_chr", tempChr, "_v3_qc.gds"), allow.duplicate = TRUE)
  
  # get all the SNPs for this run
  snpDF <- data.frame(variant_id = seqGetData(gdsfile, "variant.id"), chromosome = seqGetData(gdsfile, "chromosome"),
                      position = seqGetData(gdsfile, "position"), allele = seqGetData(gdsfile, "allele"),
                      RS = seqGetData(gdsfile, "annotation/id")) %>% 
    filter(chromosome == tempChr)
  
  # pick out the indices
  start_idx <- min(which(snpDF$position >= gene_info$txStart[gene_it] - buffer))
  end_idx <- max(which(snpDF$position <= gene_info$txEnd[gene_it] + buffer))
  
  # stop if only one SNP
  q <- end_idx - start_idx + 1
  resultsDF$q[gene_it] <- q
  if (q <= 1){
    seqResetFilter(gdsfile)
    seqClose(gdsfile)
    cat("q = ", q, ", so skip gene ", gene_it, "\n")
  }
  if (q <= 1){
    next
  }
  
  # q with no buffer
  startIdxNB <- min(which(snpDF$position >= gene_info$txStart[gene_it]))
  endIdxNB <- max(which(snpDF$position <= gene_info$txEnd[gene_it]))
  resultsDF$qNB[gene_it] <- endIdxNB - startIdxNB + 1
  
  # reset the filter
  seqResetFilter(gdsfile) # Will return selected samples and selected variants
  # set the filter
  snpsToFilter <- snpDF$variant_id[start_idx:end_idx]
  seqSetFilter(gdsfile, variant.id=snpsToFilter) # Will return selected variants
  
  # get the sample IDs
  sampleID <- seqGetData(gdsfile, "sample.id")
  
  # get genotypes 
  seqData <- SeqArray::seqGetData(gdsfile, "annotation/format/DS")
  # there is an issue with dplyr where some of the calls don't work right when the matrix has more than 2^31-1 elements
  totalElements <- prod(dim(seqData)) + prod(dim(allCov))
  if (totalElements >= 2^31 - 1) {
    
    cat("TotalElements = ", totalElements, ", so split it.", "\n")
    splitNum <- ceiling(totalElements / (2^31 - 1))
    splitSize <- ceiling(nrow(seqData) / splitNum)
    
    geno <- c()
    for(split_it in 1:splitNum) {
      startE <- (split_it - 1) * splitSize + 1
      endE <- min(split_it * splitSize, nrow(seqData))  # Old Version nrow(seqData$data) gave an error
      tempSplit <- data.frame(seqData[startE:endE, ]) %>%
        mutate(ID = sampleID[startE:endE]) %>%
        merge(., allCov, by="ID") %>%
        drop_na()
      
      # append
      geno <- bind_rows(geno, tempSplit)
      rm(tempSplit)
      cat("Split ", split_it, '\n')
    }
  } else {
    geno <- data.frame(seqData) %>%
      mutate(ID = sampleID) %>%
      inner_join(allCov, by = "ID") %>%
      drop_na()
  }
  # remove seqData for space
  rm(seqData)
  
  # drop NA rows - we filtered SNPs with too much missingness so shouldn't drop too many
  cleanG <- geno %>%
    select(paste0("X", 1:q)) %>%
    as.matrix(.)
  
  # some will have MAF 0
  MAFs <- apply(cleanG, 2, mean) / 2
  MAFpos <- which(MAFs > 0)
  cleanG <- cleanG[, MAFpos]
  MAFs <- MAFs[MAFpos]
  resultsDF$cleanq[gene_it] <- length(MAFs)
  if (length(MAFs) < 2 | length(MAFs) > max_gene_size){
    seqResetFilter(gdsfile)
    seqClose(gdsfile)
    cat("Length of MAFs = ", length(MAFs), ", so skip gene ", gene_it, "\n")
  } 
  
  if (length(MAFs) < 2 | length(MAFs) > max_gene_size){
    next
  }
  
  # flip to minor allele for weighting
  toFlip <- which(MAFs > 0.5)
  if (length(toFlip) > 0) {cleanG[, toFlip] <- apply(as.matrix(cleanG[, toFlip]), 2, flipSNPs)}
  # recalculate MAFs
  MAFs <- apply(cleanG, 2, mean) / 2
  resultsDF$rareq[gene_it] <- length(which(MAFs < 0.05))
  
  # --- weights by beta or NOT weight
  weights <- dbeta(MAFs, 1, 25)
  weightedG <- cleanG %*% diag(weights)
  # weightedG <- cleanG
  
  
  # get the covariates and times
  xMat <- as.matrix(geno %>% select("Sex", paste0("evec", 1:10))) # "Weight", "Height"
  lt <- geno$leftTime
  rt <- geno$rightTime
  lt1 <- geno$leftTime1
  rt1 <- geno$rightTime1
  deltaVec <- geno$deltaVec
  deltaVecSimple <- geno$deltaVecSimple
  gSummed <- matrix(data=apply(weightedG, 1, sum), ncol=1)
  leftCensoredInd <- geno$leftCensoredInd
  
  cat("Sample Size: ", nrow(xMat), "\n")
  # remove cleanG and geno to save space)
  rm(cleanG)
  rm(geno)
  
  # ----- crICSKAT
  # make design matrices
  nKnots <- 1
  dmats <- makeICdmat(xMat=xMat, lt=lt, rt=rt, obsInd=deltaVecSimple, quant_r = NULL, nKnots = nKnots)
  numCov <- ncol(dmats$right_dmat)
  
  cat("Now running crICSKAT", "\n")
  set.seed(0)
  iter <- 0
  solPartial <- NA
  init_1 <- c(0.1, -0.1, -2)
  init_beta <- c(runif(numCov-3, -aID4, aID4), init_1, 
                 runif(numCov-3, -aID4, aID4), init_1, 0.001)
  cat("Initial: ", init_beta, "\n")
  method <- method # Broyden or Newton
  
  # while(iter < 30 & is.na(solPartial)[1]){
  #   solPartial <- tryCatch(nleqslv::nleqslv(x = init_beta, 
  #                                           fn=scoreEqSpline, leftDmat = dmats$left_dmat,  rightDmat = dmats$right_dmat, leftTimes = lt,
  #                                           deltaVec = deltaVec, gSummed = gSummed, gMat=NULL, estG = FALSE, 
  #                                           control=list(allowSingular=TRUE, maxit=300)), method = method, 
  #                          error = function(e){return(as.numeric(NA))})
  #   iter <- iter + 1
  #   
  #   # Add some variation to initial value
  #   set.seed(iter)
  #   init_beta <- init_beta + rnorm(length(init_beta), 0, 0.01)
  # }
  # 
  # 
  # 
  # if(is.na(solPartial[1])){
  #   resultsDF$crskatp[gene_it] <- NA
  #   resultsDF$crburdenp[gene_it] <- NA
  #   resultsDF$crerr[gene_it] <- "nleqslvError"
  #   resultsDF$crInitIter[gene_it] <- NA
  # } else{
  #   resultsDF$crInitIter[gene_it] <- solPartial$iter
  #   # Ugamma with null coefficients
  #   forUgamma <- scoreEqSpline(x=solPartial$x, leftDmat = dmats$left_dmat, rightDmat = dmats$right_dmat, 
  #                              leftTimes=lt, deltaVec = deltaVec,
  #                              gSummed = gSummed, gMat=weightedG, estG = FALSE)
  #   
  #   # information
  #   # partial information, under null
  #   iMatPartial <- calcInfo(leftDmat=dmats$left_dmat, rightDmat=dmats$right_dmat, leftTimes=lt, 
  #                           theta1=solPartial$x[1:numCov], theta2=solPartial$x[(numCov+1):(2*(numCov)+1)], deltaVec=deltaVec,
  #                           gSummed=gSummed, gMat = NULL)
  #   
  #   # Igt
  #   forIgt <- calcInfo(leftDmat=dmats$left_dmat, rightDmat=dmats$right_dmat, leftTimes=lt, 
  #                      theta1=c(solPartial$x[1:numCov], rep(0, q)), theta2=solPartial$x[(numCov+1):(2*(numCov)+1)], deltaVec = deltaVec,
  #                      gSummed=gSummed, gMat=weightedG)
  #   
  #   # try SKAT
  #   skatQ <- t(forUgamma[(numCov+1):(numCov+q)]) %*% forUgamma[(numCov+1):(numCov+q)]
  #   burdenQ <- (sum(forUgamma[(numCov+1):(numCov+q)]))^2
  #   Itt <- -iMatPartial
  #   Igg <- -forIgt[(numCov+1):(numCov+q), (numCov+1):(numCov+q)]
  #   Igt <- -forIgt[c(1:numCov, ((numCov+1)+q):nrow(forIgt)), (numCov+1):(numCov+q)]
  #   sig_mat <- Igg - t(Igt) %*% solve(Itt) %*% (Igt)
  # 
  #   lambdaQ <- tryCatch(eigen(sig_mat)$value, error=function(e) "error")
  #   if(lambdaQ[1] != "error"){
  #     p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq
  #     # as noted in the CompQuadForm documentation, sometimes you need to play with acc or lim parameters
  #     # to get a p-value between 0 and 1
  #     if (!is.na(p_SKAT)) {
  #       if (p_SKAT > 1 | p_SKAT <= 0) {
  #         paramDF <- data.frame(expand.grid(lim = c(10000, 20000, 50000), acc=c(1e-7, 1e-6, 1e-5, 1e-4)))
  #         paramCounter <- 1
  #         while(p_SKAT > 1 | p_SKAT <= 0) {
  #           tempLim <- paramDF$lim[paramCounter]
  #           tempAcc <- paramDF$acc[paramCounter]
  #           p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=tempAcc, lim=tempLim)$Qq
  #           paramCounter <- paramCounter + 1
  #           if (paramCounter > nrow(paramDF)) {break}
  #         }
  #         errCode <- 22
  #         errMsg <- "Had to adjust parameters on CompQuadForm"
  #       }
  #     }
  #     B_burden = burdenQ / sum(sig_mat);
  #     p_burden = 1 - stats::pchisq(B_burden, df = 1)
  #     
  #     
  #     resultsDF$crskatp[gene_it] <- p_SKAT
  #     resultsDF$crburdenp[gene_it] <- p_burden
  #     resultsDF$crerr[gene_it] <- 0
  # 
  #   }else{
  #     resultsDF$crskatp[gene_it] <- NA
  #     resultsDF$crburdenp[gene_it] <- NA
  #     resultsDF$crerr[gene_it] <- "eigenError"
  #     cat("Eigen sig_mat got an error", "\n")
  #   }
  #   
  #   # Save memory
  #   rm(solPartial, forUgamma, iMatPartial, forIgt, Itt, Igg, Igt, sig_mat)
  #   
  #     
  #   }
  
  # check memory	
  cat("Size of weightedG: ", capture.output(pryr::object_size(weightedG)), "\n")
  cat("Memory used: ", capture.output(pryr::mem_used()), "\n") 
  cat("Now running ICSKAT", "\n")
  
  
  
  # ----- ICSKAT
  # obs_ind <- deltaVecSimple
  obs_ind <- ifelse(deltaVec == 1, 1, 0)    # Treat cause2 as censoring
  tpos_ind <- leftCensoredInd
  dmats <- ICSKAT::make_IC_dmat(xMat = xMat, lt = lt1, rt = rt1, obs_ind = obs_ind, tpos_ind = tpos_ind, nKnots = 1)
  null_fit <- ICSKAT::ICSKAT_fit_null(init_beta=c(runif(numCov, 0, 0)), 
                                      lt=lt1, rt=rt1,left_dmat=dmats$left_dmat, right_dmat=dmats$right_dmat, 
                                      obs_ind=obs_ind, tpos_ind=tpos_ind)
  
  cat("ICSKAT:Initial: ", c(runif(numCov, 0, 0)), "\n")
  
  if(is.na(null_fit$beta_fit)[1]){
    resultsDF$skatp[gene_it] <- NA
    resultsDF$burdenp[gene_it] <- NA
    resultsDF$SKATerr[gene_it] <- "ICSKATnullError"
  }else{
    
    # get pvalue
    skat_output <- ICSKAT::ICskat(left_dmat=dmats$left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                                  right_dmat=dmats$right_dmat, gMat=weightedG, lt=lt1, rt=rt1,
                                  null_beta=as.numeric(null_fit$beta_fit), Itt=null_fit$Itt)
    
    resultsDF$skatp[gene_it] <- skat_output$p_SKAT
    resultsDF$burdenp[gene_it] <- skat_output$p_burden
    resultsDF$complex[gene_it] <- skat_output$complex
    resultsDF$SKATerr[gene_it] <- 0
    
    # SKATO
    cat("Size of weightedG: ", capture.output(pryr::object_size(weightedG)), "\n")
    cat("Memory used: ", capture.output(pryr::mem_used()), "\n") 
    cat("Now running SKATO", "\n")
    
    
    skatO <- tryCatch(ICSKAT::ICSKATO(icskatOut = skat_output),
                      error = function(e) NA)
    if(is.na(skatO)[1]){
      resultsDF$SKATOp[gene_it] <- NA
      resultsDF$SKATOr[gene_it] <- NA
      resultsDF$SKATOdavies[gene_it] <- NA
      resultsDF$SKATOerr[gene_it] <- "ICSKATOError"
    }else{
      resultsDF$SKATOp[gene_it] <- skatO$pval
      resultsDF$SKATOr[gene_it] <- skatO$r
      resultsDF$SKATOdavies[gene_it] <- skatO$intDavies
      resultsDF$SKATOerr[gene_it] <- skatO$err
      
    }
    
  }
  
  
  # check memory	
  cat("Size of weightedG: ", capture.output(pryr::object_size(weightedG)), "\n")
  cat("Memory used: ", capture.output(pryr::mem_used()), "\n") 
  
  rm(null_fit, weightedG)
  
  
  # checkpoint
  cat("Done with ", gene_it, " genes.", "\n")
  cat("\n")
  
  # close gds
  seqClose(gdsfile)
}

# write results
setwd(outputDir)
write.table(resultsDF, paste0(outRoot, aID, "_", method, "_", aID4, ".txt"), append=F, quote=F, row.names=F, col.names=T) 

cat("Job ", aID, " Done! ")
