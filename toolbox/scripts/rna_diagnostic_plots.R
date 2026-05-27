library(ggplot2)
library(gridExtra)

if (!is.null(snakemake@log[[1]])) {
    log_file <- file(snakemake@log[[1]], open = "wt")
    sink(log_file, type = "output")
    sink(log_file, type = "message")
}

analyser <- snakemake@params$analyser
bg_color <- snakemake@params$bg_color

# Configure theme background element
bg_element <- if (bg_color %in% c("transparent", "na", "")) {
    element_blank()
} else {
    element_rect(fill = bg_color, color = NA)
}

message("Generating mathematical diagnostics for: ", analyser, " with background color: ", bg_color)

# 1. Load Data
df <- read.csv(snakemake@input$csv, row.names = 1, check.names = FALSE)

# Standardize names across tools for modeling landscape evaluation
if (analyser == "deseq2") {
    df$mean_abundance <- df$baseMean
    df$log2FC         <- df$log2FoldChange
    df$raw_pvalue     <- df$pvalue
} else { # edgeR
    df$mean_abundance <- df$logCPM
    df$log2FC         <- df$log2FC
    df$raw_pvalue     <- df$PValue
}

# Drop any incomplete rows safely
df_clean <- df[!is.na(df$raw_pvalue) & !is.na(df$log2FC) & !is.na(df$mean_abundance), ]

# 2. Build MA Plot
x_aes <- if(analyser == "deseq2") aes(x = log10(mean_abundance + 1), y = log2FC) else aes(x = mean_abundance, y = log2FC)
x_label <- if(analyser == "deseq2") "Log10(Base Mean + 1)" else "Log2 Counts Per Million (LogCPM)"

ma_plot <- ggplot(df_clean, x_aes) +
    geom_point(alpha = 0.4, size = 1, color = "darkslategrey") +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", alpha = 0.6) +
    labs(
        title = paste("MA Plot:", toupper(analyser)),
        x = x_label,
        y = "Log2 Fold Change"
    ) +
    theme_minimal() +
    theme(
        panel.background  = bg_element,
        plot.background   = bg_element,
        legend.background = bg_element
    )

# 3. Build P-Value Histogram Plot
hist_plot <- ggplot(df_clean, aes(x = raw_pvalue)) +
    geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white", boundary = 0) +
    labs(
        title = "P-value Frequency Distribution",
        x = "Raw Unadjusted P-value",
        y = "Frequency Count"
    ) +
    theme_minimal() +
    theme(
        panel.background  = bg_element,
        plot.background   = bg_element,
        legend.background = bg_element
    )

# 4. Save Plots Independently as requested by rule declarations
ggsave(
    filename = snakemake@output$ma,
    plot = ma_plot,
    width = 6,
    height = 5,
    dpi = 300
)

ggsave(
    filename = snakemake@output$hist,
    plot = hist_plot,
    width = 6,
    height = 5,
    dpi = 300
)

message("Diagnostic plots exported successfully.")

if (!is.null(snakemake@log[[1]])) {
    sink(type = "message")
    sink(type = "output")
}