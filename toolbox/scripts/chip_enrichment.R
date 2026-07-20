# Redirect output to the log file specified by Snakemake
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(ggplot2)

# =============================================================================
# 1. Load Data and Standardize IDs
# =============================================================================
message("Loading annotated peaks...")
peaks_df <- read.csv(snakemake@input[["annotated_csv"]], stringsAsFactors = FALSE)

if (!"geneId" %in% colnames(peaks_df)) {
    stop("Error: 'geneId' column not found in the input CSV file.")
}

# Clean prefix from string
peaks_df$geneId <- gsub("^gene:", "", peaks_df$geneId)

# Extract clean unique list
raw_genes <- peaks_df %>%
    filter(!is.na(geneId) & geneId != "") %>%
    pull(geneId) %>%
    unique()

orgdb_name <- snakemake@params[["orgdb"]]
library(orgdb_name, character.only = TRUE)
orgdb_obj <- get(orgdb_name)

# --- ID Translation Layer for Yeast ---
if (orgdb_name == "org.Sc.sgd.db") {
    message("Yeast dataset detected. Translating Systematic ORFs to Entrez IDs...")
    
    # Map ORFs to Entrez (NCBI) IDs using bitr
    gene_conversion <- bitr(
        raw_genes,
        fromType = "ORF",
        toType   = "ENTREZID",
        OrgDb    = orgdb_obj,
        drop     = FALSE
    )
    
    # GO uses the systematic ORFs directly
    go_genes  <- raw_genes
    go_key    <- "ORF"
    
    # KEGG strictly requires the mapped numerical Entrez IDs
    kegg_genes <- gene_conversion %>% 
        filter(!is.na(ENTREZID) & ENTREZID != "") %>% 
        pull(ENTREZID) %>% 
        unique()
    kegg_key   <- "ncbi-geneid"
} else {
    # Default behavior for Human/Other systems where geneId is already Entrez
    go_genes   <- raw_genes
    kegg_genes <- raw_genes
    go_key     <- "ENTREZID"
    kegg_key   <- "ncbi-geneid"
}

message(paste("Genes available for GO analysis:", length(go_genes)))
message(paste("Genes available for KEGG analysis:", length(kegg_genes)))

# =============================================================================
# 2. Enrichment Execution
# =============================================================================
message("Running GO enrichment...")
go_res <- enrichGO(
    gene          = go_genes,
    OrgDb         = orgdb_obj,
    keyType       = go_key,
    ont           = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff  = as.numeric(snakemake@params[["p_cutoff"]]),
    qvalueCutoff  = as.numeric(snakemake@params[["q_cutoff"]])
)

message("Running KEGG enrichment...")
kegg_res <- enrichKEGG(
    gene          = kegg_genes,
    organism      = snakemake@params[["kegg_organism"]],
    keyType       = kegg_key,
    pAdjustMethod = "BH",
    pvalueCutoff  = as.numeric(snakemake@params[["p_cutoff"]]),
    qvalueCutoff  = as.numeric(snakemake@params[["q_cutoff"]])
)

write.csv(as.data.frame(go_res), snakemake@output[["go_results"]], row.names = FALSE)
write.csv(as.data.frame(kegg_res), snakemake@output[["kegg_results"]], row.names = FALSE)

# =============================================================================
# 3. Render Dotplots
# =============================================================================
dotplotheight <- as.numeric(snakemake@params$dotplot_height)
bg_color      <- snakemake@params$bg_color

bg_element <- if (is.null(bg_color) || bg_color %in% c("transparent", "na", "")) {
    element_blank()
} else {
    element_rect(fill = bg_color, color = NA)
}

if (nrow(as.data.frame(go_res)) > 0) {
    p_go <- dotplot(go_res, showCategory = 20, title = paste("GO Enrichment -", snakemake@wildcards[["sample"]])) +
        theme_bw() +
        theme(panel.background = bg_element, plot.background = bg_element, legend.background = bg_element)
    ggsave(snakemake@output[["go_dotplot"]], plot = p_go, width = 8, height = dotplotheight, dpi = 300)
} else {
    file.create(snakemake@output[["go_dotplot"]])
}

if (nrow(as.data.frame(kegg_res)) > 0) {
    p_kegg <- dotplot(kegg_res, showCategory = 20, title = paste("KEGG Enrichment -", snakemake@wildcards[["sample"]])) +
        theme_bw() +
        theme(panel.background = bg_element, plot.background = bg_element, legend.background = bg_element)
    ggsave(snakemake@output[["kegg_dotplot"]], plot = p_kegg, width = 8, height = dotplotheight, dpi = 300)
} else {
    file.create(snakemake@output[["kegg_dotplot"]])
}

message("Enrichment analysis completed successfully.")