# Assignment 4 - Microarray Preprocessing Workflow
# Author: Bilal
# Email: bilal.ai.dev@gmail.com
# Dataset: GSE70970

# Step 1: Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("affy", "arrayQualityMetrics", "limma", "GEOquery"))
install.packages(c("ggplot2", "pheatmap"))

# Step 2: Create working directory
dir.create("GSE70970_RAW", showWarnings = FALSE)
setwd("GSE70970_RAW")

# Step 3: Download and extract CEL files
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70970&format=file"
download.file(url, destfile = "GSE70970_RAW.zip", mode = "wb")
unzip("GSE70970_RAW.zip", exdir = "CEL_Files")

# Step 4: Load CEL files
library(affy)
setwd("CEL_Files")
rawData <- ReadAffy()

# Step 5: Quality control before normalization
setwd("..")
dir.create("Results", showWarnings = FALSE)
library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = rawData, outdir = "Results/QC_Before_Normalization", force = TRUE)

# Step 6: Normalize data
normData <- rma(rawData)
data.norm <- exprs(normData)

# Step 7: Quality control after normalization
arrayQualityMetrics(expressionset = normData, outdir = "Results/QC_After_Normalization", force = TRUE)

# Step 8: Filter low-intensity probes
cutoff <- rowMeans(data.norm) > log2(20)
filtered.data <- data.norm[cutoff, ]
cat("Probes before filtering:", nrow(data.norm), "\n")
cat("Probes after filtering:", nrow(filtered.data), "\n")

# Step 9: PCA Plot
pca <- prcomp(t(filtered.data))
pca_data <- as.data.frame(pca$x)
pca_data$Group <- pData(normData)$characteristics_ch1.1
pdf("Results/PCA_plot.pdf")
library(ggplot2)
ggplot(pca_data, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA After Normalization")
dev.off()

# Step 10: Boxplot after normalization
pdf("Results/Boxplot_After_Normalization.pdf")
boxplot(data.norm, main = "Boxplot After Normalization", las = 2, col = "lightblue")
dev.off()

# Step 11: Relabel groups (Normal vs Cancer)
pdata <- pData(normData)
pdata$Group <- ifelse(grepl("normal", pdata$characteristics_ch1.1, ignore.case = TRUE), "Normal", "Cancer")
table(pdata$Group)