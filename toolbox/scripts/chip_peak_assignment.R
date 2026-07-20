# Setup logging streams correctly
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

library(ChIPseeker)
library(GenomicFeatures)

# =============================================================================
# 1. Environment and Resource Preparation
# =============================================================================
peak_file  <- snakemake@input$peak_bed
gtf_file   <- snakemake@input$gtf
orgdb_name <- snakemake@params$orgdb

if (!requireNamespace(orgdb_name, quietly = TRUE)) {
    stop(paste("Critical Error: Failed to load required annotation package:", orgdb_name))
}
library(orgdb_name, character.only = TRUE)

message("Building TxDb object from GTF...")
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)

# =============================================================================
# 2. Peak Annotation Execution
# =============================================================================
tss_start <- as.integer(snakemake@params$tss_region_start)
tss_end   <- as.integer(snakemake@params$tss_region_end)

message("Executing coordinate mapping using ChIPseeker...")
peak_annotation <- ChIPseeker::annotatePeak(
    peak      = peak_file,
    TxDb      = txdb,
    tssRegion = c(tss_start, tss_end),
    annoDb    = orgdb_name,
    verbose   = FALSE
)

# =============================================================================
# 3. Export Data and Visualizations
# =============================================================================
annotation_df <- as.data.frame(peak_annotation)
write.csv(annotation_df, file = snakemake@output$annotation_csv, row.names = FALSE)

png(filename = snakemake@output$plot_distribution, width = 800, height = 400, res = 120)
print(ChIPseeker::plotAnnoBar(peak_annotation))
dev.off()

message("ChIP peak assignment completed successfully.")