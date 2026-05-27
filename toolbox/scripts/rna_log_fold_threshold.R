library(ggplot2)

if (!is.null(snakemake@log[[1]])) {
    log_file <- file(snakemake@log[[1]], open = "wt")
    sink(log_file, type = "output")
    sink(log_file, type = "message")
}

# 1. Parse Parameters
fdr_cutoff <- as.numeric(snakemake@params$fdr)
lfc_cutoff <- as.numeric(snakemake@params$lfc)
analyser   <- snakemake@params$analyser
bg_color   <- snakemake@params$bg_color

bg_element <- if (bg_color %in% c("transparent", "na", "")) {
    element_blank()
} else {
    element_rect(fill = bg_color, color = NA)
}

message("Running threshold filtering for: ", analyser)
message("Cutoffs: FDR <= ", fdr_cutoff, " | |Log2FC| >= ", lfc_cutoff)

# 2. Load Data
df <- read.csv(snakemake@input$csv, row.names = 1, check.names = FALSE)

# Standardize column naming across tools
if (analyser == "deseq2") {
    df$p_val_adj <- df$padj
    df$log2FC    <- df$log2FoldChange
} else { # edgeR
    df$p_val_adj <- df$FDR
    df$log2FC    <- df$log2FC
}

# Drop structural NA rows safely
df <- df[!is.na(df$p_val_adj) & !is.na(df$log2FC), ]

# 3. Categorize Significance Matrix
df$Significance <- "Not Significant"
df$Significance[df$p_val_adj <= fdr_cutoff & df$log2FC >= lfc_cutoff]  <- "Upregulated"
df$Significance[df$p_val_adj <= fdr_cutoff & df$log2FC <= -lfc_cutoff] <- "Downregulated"
df$Significance <- factor(df$Significance, levels = c("Not Significant", "Upregulated", "Downregulated"))

# 4. Save Filtered Dataset
sig_df <- df[df$Significance != "Not Significant", ]
write.csv(sig_df, file = snakemake@output$sig_csv, row.names = TRUE)
message("Saved ", nrow(sig_df), " significant genes to: ", snakemake@output$sig_csv)

# 5. Generate Volcano Plot
volcano <- ggplot(df, aes(x = log2FC, y = -log10(p_val_adj), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
        "Not Significant" = "grey70",
        "Upregulated"     = "firebrick",
        "Downregulated"   = "dodgerblue"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
        title = paste("Volcano Plot - Active Analyzer:", toupper(analyser)),
        subtitle = paste0("Thresholds: FDR <= ", fdr_cutoff, ", |Log2FC| >= ", lfc_cutoff),
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(
        legend.position   = "right",
        panel.background  = bg_element,
        plot.background   = bg_element,
        legend.background = bg_element
    )

# 6. Save Plot Asset
ggsave(
    filename = snakemake@output$volcano,
    plot = volcano,
    width = 7,
    height = 6,
    dpi = 300
)

if (!is.null(snakemake@log[[1]])) {
    sink(type = "message")
    sink(type = "output")
}