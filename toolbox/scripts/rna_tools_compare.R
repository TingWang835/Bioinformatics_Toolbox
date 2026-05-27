library(ggplot2)

if (!is.null(snakemake@log[[1]])) {
    log_file <- file(snakemake@log[[1]], open = "wt")
    sink(log_file, type = "output")
    sink(log_file, type = "message")
}

# 1. Parse Parameters
lfc_cutoff <- snakemake@params$lfc
fdr_cutoff <- snakemake@params$fdr
bg_color   <- snakemake@params$bg_color

bg_element <- if (bg_color %in= c("transparent", "na", "")) {
    element_blank()
} else {
    element_rect(fill = bg_color, color = NA)
}

# 2. Load Data
de_deseq2 <- read.csv(snakemake@input$deseq2, row.names = 1, check.names = FALSE)
de_edger  <- read.csv(snakemake@input$edger, check.names = FALSE)

rownames(de_edger) <- de_edger$Geneid

# 3. Merge Data on Common Genes
character_genes <- intersect(rownames(de_deseq2), rownames(de_edger))
merged_df <- data.frame(
    Geneid         = character_genes,
    DESeq2_log2FC  = de_deseq2[character_genes, "log2FoldChange"],
    DESeq2_padj    = de_deseq2[character_genes, "padj"],
    edgeR_log2FC   = de_edger[character_genes, "logFC"],
    edgeR_FDR      = de_edger[character_genes, "FDR"]
)

# Apply Significance Categorization
sig_deseq2_raw <- !is.na(merged_df$DESeq2_padj) & 
                  merged_df$DESeq2_padj <= fdr_cutoff & 
                  !is.na(merged_df$DESeq2_log2FC) & 
                  abs(merged_df$DESeq2_log2FC) >= lfc_cutoff

sig_edger_raw  <- !is.na(merged_df$edgeR_FDR) & 
                  merged_df$edgeR_FDR <= fdr_cutoff & 
                  !is.na(merged_df$edgeR_log2FC) & 
                  abs(merged_df$edgeR_log2FC) >= lfc_cutoff

merged_df$is_sig_deseq2 <- ifelse(is.na(sig_deseq2_raw), FALSE, sig_deseq2_raw)
merged_df$is_sig_edger  <- ifelse(is.na(sig_edger_raw), FALSE, sig_edger_raw)

merged_df$Significance <- "Not Significant"
merged_df$Significance[merged_df$is_sig_deseq2 & !merged_df$is_sig_edger]  <- "DESeq2 Only"
merged_df$Significance[!merged_df$is_sig_deseq2 & merged_df$is_sig_edger]  <- "edgeR Only"
merged_df$Significance[merged_df$is_sig_deseq2 & merged_df$is_sig_edger]    <- "Both"

merged_df$Significance <- factor(merged_df$Significance, levels = c("Not Significant", "DESeq2 Only", "edgeR Only", "Both"))

cat("\n--- Cross-Tool Statistical Signatures Summary ---\n")
print(table(merged_df$Significance))

# 4. Export Combined Dataset
write.csv(merged_df, file = snakemake@output$csv, row.names = FALSE)

# 5. Calculate Pearson Correlation Coefficient
cor_coef <- cor(merged_df$DESeq2_log2FC, merged_df$edgeR_log2FC, method = "pearson", use = "complete.obs")

# 6. Generate Correlation Scatter Plot
plot <- ggplot(merged_df, aes(x = DESeq2_log2FC, y = edgeR_log2FC, color = Significance)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    scale_color_manual(values = c(
        "Not Significant" = "grey75",
        "DESeq2 Only"     = "forestgreen",
        "edgeR Only"      = "dodgerblue",
        "Both"            = "firebrick"
    )) +
    labs(
        title = paste("Log2 Fold Change Correlation (r =", round(cor_coef, 3), ")"),
        subtitle = paste0("Thresholds: FDR <= ", fdr_cutoff, ", |Log2FC| >= ", lfc_cutoff),
        x = "DESeq2 Log2 Fold Change",
        y = "edgeR Log2 Fold Change"
    ) +
    theme_minimal() +
    theme(
        legend.position   = "right",
        panel.background  = bg_element,
        plot.background   = bg_element,
        legend.background = bg_element
    )

ggsave(filename = snakemake@output$venn, plot = plot, width = 7, height = 6, dpi = 300)

if (!is.null(snakemake@log[[1]])) {
    sink(type = "message")
    sink(type = "output")
}