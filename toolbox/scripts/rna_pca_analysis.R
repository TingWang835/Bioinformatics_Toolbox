library(ggplot2)
library(DESeq2)
library(edgeR)

if (!is.null(snakemake@log[[1]])) {
    log_file <- file(snakemake@log[[1]], open = "wt")
    sink(log_file, type = "output")
    sink(log_file, type = "message")
}

# 1. Parse Parameters
analyser     <- snakemake@params$analyser
aligner_type <- snakemake@params$aligner_type
batch_col    <- snakemake@params$batch
bg_color     <- snakemake@params$bg_color

# Route to element_blank() if setting matches transparency keywords or is empty
bg_element <- if (bg_color %in% c("transparent", "na", "")) {
    element_blank()
} else {
    element_rect(fill = bg_color, color = NA)
}

# 2. Load Metadata
run_info <- read.csv(snakemake@input$runinfo, row.names = 1, check.names = FALSE)
# Ensure design formula columns are treated as factors
run_info$condition <- as.factor(run_info$condition)
if (batch_col != "" && batch_col %in% colnames(run_info)) {
    run_info[[batch_col]] <- as.factor(run_info[[batch_col]])
}

# 3. Load Count Matrix & Compute Normalized Variance Stabilized Matrix
if (analyser == "deseq2") {
    counts <- read.csv(snakemake@input$counts, row.names = 1, check.names = FALSE)
    counts <- counts[, rownames(run_info)]
    
    dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = run_info, design = ~condition)
    vsd <- rlog(dds, blind = TRUE)
    pca_data <- plotPCA(vsd, intgroup = c("condition", if(batch_col != "") batch_col else "condition"), returnData = TRUE)
    percent_var <- round(100 * attr(pca_data, "percentVar"))
    
} else { # edgeR
    counts <- read.csv(snakemake@input$counts, row.names = 1, check.names = FALSE)
    counts <- counts[, rownames(run_info)]
    
    dge <- DGEList(counts = counts, samples = run_info)
    dge <- calcNormFactors(dge)
    
    log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
    pca_res <- prcomp(t(log_cpm))
    pca_data <- data.frame(
        PC1 = pca_res$x[, 1],
        PC2 = pca_res$x[, 2],
        condition = run_info$condition
    )
    if (batch_col != "" && batch_col %in% colnames(run_info)) {
        pca_data[[batch_col]] <- run_info[[batch_col]]
    }
    
    vars <- pca_res$sdev^2
    percent_var <- round(100 * (vars / sum(vars)))
}

# 4. Generate the PCA Plot
color_aes <- if (batch_col != "" && batch_col %in% colnames(run_info)) aes_string(x = colnames(pca_data)[1], y = colnames(pca_data)[2], color = "condition", shape = batch_col) else aes_string(x = colnames(pca_data)[1], y = colnames(pca_data)[2], color = "condition")

pca_plot <- ggplot(pca_data, color_aes) +
    geom_point(size = 4, alpha = 0.8) +
    labs(
        title = paste("PCA Plot - Aligner Matrix (", toupper(analyser), ")"),
        x = paste0("PC1: ", percent_var[1], "% variance"),
        y = paste0("PC2: ", percent_var[2], "% variance")
    ) +
    theme_minimal() +
    theme(
        panel.background  = bg_element,
        plot.background   = bg_element,
        legend.background = bg_element
    )

# 5. Save Output
ggsave(
    filename = snakemake@output$png,
    plot = pca_plot,
    width = 7,
    height = 6,
    dpi = 300
)

if (!is.null(snakemake@log[[1]])) {
    sink(type = "message")
    sink(type = "output")
}