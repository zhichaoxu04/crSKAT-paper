library(foreach)
require(doSNOW)
library(parallel)
library(iterators)
library(snow)
library(tidyverse)
library(data.table)
library(SeqArray)


args <- commandArgs(trailingOnly=TRUE)
aID1 <- as.numeric(args[1])   # DATASET

# Initialize cluster
num_cores <- 80  # Set the number of cores to use
cl <- makeCluster(num_cores)

set.seed(2024)
SampleSplit <- sample(1:300000, aID1, replace = FALSE)

# Define the replicate task function
replicate_task <- function(chr) {
  # Set working directory and open GDS file
  setwd("/rsrch3/scratch/biostatistics/zxu7/Rotation/RS/RDA/ukb_gds")
  gdsfile <- seqOpen(paste0("ukb_imp_chr", chr, "_v3_qc.gds"), allow.duplicate = TRUE)
  variant_id <- seqGetData(gdsfile, "variant.id")
  
  # Randomly select starting point for 100 consecutive SNPs
  start_snp <- sample(1:(max(variant_id) - 99), 1)
  snpsToFilter <- start_snp:(start_snp + 99)
  seqSetFilter(gdsfile, variant.id = snpsToFilter) # Will return selected variants
  
  # Get genotypes
  seqData <- SeqArray::seqGetData(gdsfile, "annotation/format/DS")
  seqClose(gdsfile)
  
  seqDataSub <- seqData[SampleSplit, ]
  
  cleanG <- as.data.frame(seqDataSub) %>%
    drop_na()
  
  rm(seqData, seqDataSub)  # Clean up to save memory
  
  # Calculate correlations with the first SNP and the 2nd, 5th, 10th, ..., 100th SNPs
  correlations <- sapply(c(2, 5, seq(10, 100, by = 5)), function(i) {
    cor(cleanG[, 1], cleanG[, i], use = "pairwise.complete.obs")
  })
  
  return(list(start_snp = start_snp, correlations = correlations))
}

# Export the necessary variables and functions to the cluster
clusterExport(cl, varlist = c("replicate_task", "SampleSplit", "aID1"))

# Load necessary libraries on all nodes
clusterEvalQ(cl, {
  library(dplyr)
  library(SeqArray)
  library(tidyverse)
})

# Run 10 million replications
n_replications <- 1e+07
result <- parLapply(cl, X = rep(9, n_replications), fun = replicate_task)

# Stop and release resources
stopCluster(cl)

# Extract start_snp and correlations
start_snps <- sapply(result, `[[`, "start_snp")
cor_matrix <- do.call(rbind, lapply(result, `[[`, "correlations"))

print(start_snps)

# Calculate quantiles for each correlation set
quantiles_list <- apply(cor_matrix, 2, function(x) {
  quantile(x, probs = seq(0, 1, length = 11), na.rm = TRUE)
})

# Convert the quantiles list to a data frame for easier export
quantiles_df <- as.data.frame(t(quantiles_list))
colnames(quantiles_df) <- paste0("Q", 0:10 * 10, "%")

# Format the quantiles_df to scientific notation with 3 digits after decimal point
quantiles_df[] <- lapply(quantiles_df, function(x) format(x, scientific = TRUE, digits = 3))

# Print the results
print(quantiles_df)

# Write the quantiles to a text file
write.table(quantiles_df, file = paste0("./quantiles_list_N", 
                                        aID1,
                                        ".txt"), 
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)



