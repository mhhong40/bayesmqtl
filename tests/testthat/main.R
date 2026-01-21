## Create test dataset from UCSCXenaTools for mixture model fit test
library(UCSCXenaTools)

clinical <- XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "clinical") %>%
  XenaFilter(filterDatasets = "BRCA")

clinical_download <- XenaQuery(clinical) %>%
  XenaDownload()

clinical <- XenaPrepare(clinical_download)
# saveRDS(clinical, file = "TCGA_BRCA_clinical.RDS")

raw_methylation <- XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "HumanMethylation450") %>%
  XenaFilter(filterDatasets = "BRCA")

raw_methylation_download <- XenaQuery(raw_methylation) %>%
  XenaDownload()

raw_methylation <- XenaPrepare(methylation_download)
# saveRDS(raw_methylation, file = "TCGA_BRCA_raw_methylation.RDS")

## Keep only females with luminal A subtype (ER+, HER2-) who are also PR+
# This turns out to be n_init = 901
patients <- clinical[clinical$Gender_nature2012 == "FEMALE" & clinical$ER_Status_nature2012 == "Positive" & clinical$HER2_Final_Status_nature2012 == "Negative" & clinical$PR_Status_nature2012 == "Positive", ]
patient_IDs <- patients$sampleID
methylation <- raw_methylation[, colnames(raw_methylation) %in% c("sample", patient_IDs)] # n_final = 255

# Remove genes that have any beta values = 0, 1
library(dplyr)
methylation <- methylation %>%
  filter(if_all(everything(), ~ . != 0)) %>%
  filter(if_all(everything(), ~ . != 1))

# Take a random sample of d_final = 20 CpGs for package testing
set.seed(12345)
d_init <- nrow(methylation)
d_final <- 20
final_methylation <- methylation[c(1, sample(2:d_init, d_final, replace = FALSE)), ]
cpgs <- final_methylation$sample
final_methylation <- final_methylation[, -1]

library(minfi)
n_final <- ncol(final_methylation)
mm <- as.matrix(t(final_methylation))
colnames(mm) <- cpgs
# densityPlot(mm, sampGroups = NULL) # just to visualize
saveRDS(mm, file = "TCGA_BRCA_final_methylation_matrix.RDS")

library(maxLik)
library(profvis) # profile code
devtools::load_all()

## Behaves well on toy dataset
Y_1 <- c(rbeta(100, 1, 5), rbeta(150, 10, 2))
Y_2 <- c(rbeta(160, 4, 8), rbeta(90, 20, 3))
Y_new <- as.matrix(cbind(Y_1, Y_2)) # d = 2 for debugging stuff
fit_results <- fit_mixture_model_(Y_new, parameterization = "shape", tol = 0.0001)

mm <- data(mm) # not working for some reason idk

mm <- readRDS("TCGA_BRCA_final_methylation_matrix.RDS")

mm_sample <- mm[, c(1:2)]
fit_mm_sample <- fit_mixture_model_(mm_sample, parameterization = "shape")
labels_time <- system.time(fit_mixture_model_(mm_sample, parameterization = "shape"))

fit_mm_sample_no_labels <- fit_mixture_model_(mm_sample, parameterization = "shape")
no_labels_time <- system.time(fit_mixture_model_(mm_sample, parameterization = "shape"))

mm_sample_2 <- mm[, c(1:5)]
runtime <- system.time(fit_mm_sample_2 <- fit_mixture_model_(mm_sample_2, parameterization = "shape"))
