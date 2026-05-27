# Redirect output and errors to the log file defined in the Snakemake rule
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
if (snakemake@params$aligner_type == "pseudo") {
    library(tximport)
}

# 1. Load Data
counts_path <- snakemake@input$counts
metadata_path <- snakemake@input$runinfo

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)

# Filter metadata to match only the samples present in the counts matrix columns
metadata <- metadata[colnames(counts), , drop = FALSE]

# Ensure 'condition' is treated as a factor baseline
metadata$condition <- as.factor(metadata$condition)

# 2. Handle Multi-threading Setup inside DESeq2
library(BiocParallel)
register(MulticoreParam(snakemake@threads))

# 3. Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ condition)

# 4. Execute Analysis with parameter-driven normalization switch
if (snakemake@params$normz == "FALSE") {
    sizeFactors(dds) <- rep(1, ncol(dds))
    dds <- DESeq(dds, parallel = TRUE)
} else {
    dds <- DESeq(dds, parallel = TRUE)
}

# 5. Extract and Save Results
res <- results(dds, parallel = TRUE)
write.csv(as.data.frame(res), file = snakemake@output$csv)

# 6. Clean up sinks cleanly in reverse order
sink(type = "message")
sink()
close(log)