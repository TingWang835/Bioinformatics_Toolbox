library(clusterProfiler)
library(ggplot2)
library(dplyr)

if (!is.null(snakemake@log[[1]])) {
    log_file <- file(snakemake@log[[1]], open = "wt")
    sink(log_file, type = "output")
    sink(log_file, type = "message")
}

# 1. Parse Parameters and Resolve Organism Database
species_name <- snakemake@params$species_latin
bg_color     <- snakemake@params$bg_color

bg_element <- if (bg_color %in% c("transparent", "na", "")) {
    element_blank()
} else {
    element_rect(fill = bg_color, color = NA)
}

message("Starting functional enrichment analysis for species reference string: ", species_name)

# Use strict procedural routing to prevent eager execution of stop()
if (species_name %in% c("saccharomycescerevisiae", "scerevisiae", "yeast")) {
    org_db <- "org.Sc.sgd.db"
} else if (species_name %in% c("homosapiens", "hsapiens", "human")) {
    org_db <- "org.Hs.eg.db"
} else {
    stop(paste("Aborting: Species package mapping not configured for:", species_name))
}

# Dynamic library activation
if (!requireNamespace(org_db, quietly = TRUE)) {
    stop(paste("Aborting: Required annotation package library is missing:", org_db))
}
library(org_db, character.only = TRUE)

# 2. Load and Prepare Gene Signatures
sig_df <- read.csv(snakemake@input$sig_genes, row.names = 1, check.names = FALSE)
gene_list <- unique(rownames(sig_df))

# CLEANING STEP: Strips prefixes ("gene:" or "transcript:") and trailing suffixes ("_mRNA", "-scRNA", etc.)
gene_list <- gsub("^gene:|^transcript:|(_.*|-.*)$", "", gene_list)

if (length(gene_list) == 0) {
    stop("Aborting: Input significant gene list is empty.")
}

# Determine identifier keytype based on string composition
key_type <- if (org_db == "org.Sc.sgd.db") "ORF" else "SYMBOL"
message("Using feature keytype: ", key_type, " with target size: ", length(gene_list))

# 3. Execute Gene Ontology Over-Representation Analysis (ORA)
go_res <- enrichGO(
    gene          = gene_list,
    OrgDb         = get(org_db),
    keyType       = key_type,
    ont           = "ALL",       # Evaluates BP, CC, and MF
    pAdjustMethod = "BH",        # Benjamini-Hochberg FDR correction
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
)

# 4. Handle and Export Diagnostics
dotplotheight <- snakemake@params$dotplot_height

if (is.null(go_res) || nrow(go_res) == 0) {
    message("Warning: No functional terms passed the ORA correction cutoff matrix thresholds.")
    
    write.csv(data.frame(Message = "No enriched terms found"), file = snakemake@output$go_csv, row.names = FALSE)
    
    empty_plot <- ggplot() + 
        annotate("text", x = 1, y = 1, label = "No Enriched Terms Found") + 
        theme_void() +
        theme(
            panel.background  = bg_element,
            plot.background   = bg_element
        )
    ggsave(filename = snakemake@output$go_dot, plot = empty_plot, width = 8, height = 6)
    
} else {
    # Export clean text matrix data
    write.csv(as.data.frame(go_res), file = snakemake@output$go_csv, row.names = FALSE)
    message("Enrichment results saved to: ", snakemake@output$go_csv)
    
    # Render dotplot focusing on the top biological pathways
    dot_plot <- dotplot(go_res, showCategory = 20, split = "ONTOLOGY") +
        facet_grid(ONTOLOGY ~ ., scales = "free") +
        theme_bw() +
        theme(
            panel.grid.major  = element_line(color = "grey92"),
            strip.background  = element_rect(fill = "grey95", color = "grey80"),
            axis.text.y       = element_text(size = 9),
            panel.background  = bg_element,
            plot.background   = bg_element,
            legend.background = bg_element
        )
    
    ggsave(
        filename = snakemake@output$go_dot,
        plot     = dot_plot,
        width    = 9,
        height   = dotplotheight,
        dpi      = 300
    )
    message("Enrichment visualization dotplot saved to: ", snakemake@output$go_dot)
}

if (!is.null(snakemake@log[[1]])) {
    sink(type = "message")
    sink(type = "output")
}