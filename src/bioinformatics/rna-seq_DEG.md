## DEG analysis with DESeq2 (input: featureCounts output)

Summary: create a conda R environment with RStudio, install DESeq2, then run a reproducible R script that reads featureCounts output, performs QC, differential expression, shrinkage, and exports results.

Please refer to this tutorial [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)



### Step1 Create conda env and install R / RStudio
Run in a terminal:
```bash
conda install -c conda-forge -c bioconda r-base=4.2 rstudio r-essentials -y
# launch RStudio (optional)
rstudio
```

Note: Bioconductor packages are best installed from within R using BiocManager (below).



### Step2 Install R packages (from inside R or RStudio)
Open R and run:
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2","apeglm","tximport","pheatmap","RColorBrewer","vsn"))
```
- DESeq2: core DEA
- apeglm: LFC shrinkage
- tximport: only if importing transcript-level results (not needed for featureCounts gene counts)
- pheatmap / RColorBrewer / vsn: plotting / normalization helpers



### Step3 Prepare featureCounts input

> Try by yourself
> Finally generate a counts_file (**featureCounts_counts.txt**) for next analysis

### Step4 FeatureCounts merged file

Save as `deseq2_featurecounts_DEA.R` and edit `counts_file` and `coldata` accordingly.

```r
# deseq2_featurecounts_DEA.R
library(DESeq2)
library(apeglm)
library(pheatmap)
library(RColorBrewer)
library(vsn)

# ---- user settings ----
counts_file <- "featureCounts_counts.txt" # merged output
# Define sample table: sample names (matching cols) and condition
coldata <- data.frame(
    sample = c("sampleA_rep1","sampleA_rep2","sampleB_rep1","sampleB_rep2"),
    condition = factor(c("A","A","B","B"))
)
rownames(coldata) <- coldata$sample
# ------------------------

# Read featureCounts merged file (skip comment lines starting with '#')
fc <- read.table(counts_file, header = TRUE, sep = "\t", comment.char="#", stringsAsFactors=FALSE)

# Typical featureCounts layout: Geneid, Chr, Start, End, Strand, Length, sample1, sample2, ...
# Find first counts column (usually where column names match sample names or numeric)
# If counts start at col 7:
counts <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts) <- fc$Geneid
# Optionally subset columns to match coldata
counts <- counts[, rownames(coldata)]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                                        colData = coldata,
                                        design = ~ condition)

# Prefiltering: remove very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq
dds <- DESeq(dds)

# Variance-stabilizing transformation for QC and PCA
vsd <- vst(dds, blind=FALSE)


# Results table
res <- results(dds, alpha=0.05)           # default: condition B vs A depending on level order
res <- lfcShrink(dds, coef=2, type="apeglm")  # adjust coef number to match your design
resOrdered <- res[order(res$padj), ]

# Export results
write.csv(as.data.frame(resOrdered), file="deseq2_results.csv")

# Top genes heatmap
topGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 50)
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)


# Save transformed counts for downstream plots
write.csv(as.data.frame(assay(vsd)), file="vsd_normalized_counts.csv")
```

Notes:
- Adjust `coldata` to reflect your samples and experimental design (use exact column names).
- The `coef` in lfcShrink depends on the ordering of levels in your design; check `resultsNames(dds)` to pick correct coefficient.
- Filtering threshold (here rowSums >= 10) can be changed depending on depth/replicates.




This provides a concise, reproducible pipeline for DESeq2 starting from featureCounts output. Adjust sample metadata and design formula to your experiment.