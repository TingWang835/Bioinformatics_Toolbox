library(pheatmap)

if (!is.null(snakemake@log[[1]])) {
    log_file <- file(snakemake@log[[1]], open = "wt")
    sink(log_file, type = "output")
    sink(log_file, type = "message")
}

# 1. Parse Parameters
batch_col <- snakemake@params$batch
bg_color  <- snakemake@params$bg_color

bg_color_val <- if (bg_color %in% c("transparent", "na", "")) "NA" else bg_color

# 2. Load Data Assets
counts    <- read.csv(snakemake@input$counts, row.names = 1, check.names = FALSE)
run_info  <- read.csv(snakemake@input$runinfo, row.names = 1, check.names = FALSE)
sig_genes <- read.csv(snakemake@input$sig_genes, row.names = 1, check.names = FALSE)

# Extract significant gene identifiers
sig_identifiers <- rownames(sig_genes)

if (length(sig_identifiers) == 0) {
    stop("Aborting: Zero significant genes found in the input matrix threshold file.")
}

# Match sample ordering exactly
counts <- counts[, rownames(run_info), drop = FALSE]

# Filter count matrix for significant genes
mat <- counts[rownames(counts) %in% sig_identifiers, , drop = FALSE]

if (nrow(mat) == 0) {
    stop("Aborting: None of the significant genes could be matched to the raw expression matrix IDs.")
}

# 3. Log-Transform Counts per Million for Scale Normalization
# Use a baseline prior count to prevent taking log of zero
log_mat <- log2(edgeR::cpm(mat, prior.count = 2) + 1)

# Center and scale rows (z-scores) to emphasize patterns across samples
z_mat <- t(scale(t(log_mat)))
z_mat[is.na(z_mat)] <- 0  # Handle zero-variance scale drops safely

# 4. Configure Metadata Annotation Framework
annotation_col <- data.frame(Condition = as.factor(run_info$condition))
rownames(annotation_col) <- rownames(run_info)

if (batch_col != "" && batch_col %in% colnames(run_info) && batch_col != "condition") {
    annotation_col[[batch_col]] <- as.factor(run_info[[batch_col]])
}

# 5. Generate and Export Heatmap Object
pheatmap(
    mat               = z_mat,
    annotation_col    = annotation_col,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = (nrow(z_mat) <= 60),  # Hide labels if feature density blurs grid
    show_colnames     = TRUE,
    scale             = "none",               # Pre-scaled manually above
    clustering_method = "complete",
    color             = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50),
    background_color  = bg_color_val,
    filename          = snakemake@output$png,
    width             = 8,
    height            = 8
)

if (!is.null(snakemake@log[[1]])) {
    sink(type = "message")
    sink(type = "output")
}