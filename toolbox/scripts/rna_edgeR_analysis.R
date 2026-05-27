# Redirect output and errors to the log file defined in the Snakemake rule
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(edgeR)

# 1. Load Data
counts_path <- snakemake@input$counts
metadata_path <- snakemake@input$runinfo

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)

# Filter metadata rows to match count columns precisely
metadata <- metadata[colnames(counts), , drop = FALSE]
group <- factor(metadata$condition)

# 2. Initialize DGEList container
dge <- DGEList(counts = counts, group = group)

# 3. Filter Low Count Genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# 4. Apply Normalization Method from Config
raw_norm <- snakemake@params$normz

if (raw_norm %in% c("TMM", "RLE")) {
    dge <- calcNormFactors(dge, method = raw_norm)
} else {
    dge <- calcNormFactors(dge, method = tolower(raw_norm))
}

# 5. Statistical Modeling Setup
design <- model.matrix(~ group)

if (window_system <- .Platform$OS.type == "unix") {
    options(mc.cores = snakemake@threads)
}

dge <- estimateDisp(dge, design)

# 6. Fit Model and Execute Quasi-Likelihood F-Test
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)

# 7. Format and Export Metrics
res <- topTags(qlf, n = Inf, sort.by = "none")$table
res <- cbind(Geneid = rownames(res), res)
write.csv(res, file = snakemake@output$csv, row.names = FALSE)

# 8. Clean up sinks cleanly in reverse order
sink(type = "message")
sink()
close(log)