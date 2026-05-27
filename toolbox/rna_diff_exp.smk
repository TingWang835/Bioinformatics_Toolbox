
# =============================================================================
# Helper Functions
# =============================================================================
def get_expression_input(wildcards):
    """Selects the correct combined count matrix based on aligner classification."""
    aligner = config.get("ALIGNER", "star").lower()
    if aligner in RNA_SPLICED_ALIGNERS:
        return f"{READS_DIR}/counts/all_samples.{aligner}_counts.csv"
    else:
        return f"{READS_DIR}/counts/all_samples.{aligner}_pseudo_counts.csv"

# =============================================================================
# Rules
# =============================================================================

rule DESeq2:
    """
    Executes differential gene/transcript expression analysis using DESeq2.
    """
    input:
        counts = get_expression_input,
        runinfo = f"{READS_DIR}/{data_src.lower()}_runinfo.csv"
    output:
        csv = ( 
        f"{READS_DIR}/expression/deseq2_expression.{rna_aligner}_counts.csv" 
        if rna_aligner in RNA_SPLICED_ALIGNERS else 
        f"{READS_DIR}/expression/deseq2_expression.{rna_aligner}_pseudo_counts.csv"
        )

    log:
        f"{LOG_DIR}/diff_exp/deseq2/exp_result.{rna_aligner}.log"
    conda:
        "../env/rna_diff_exp.yaml"
    threads: 4
    params:
        normz = config.get("DESEQ2_NORM", "TRUE").upper(),
        aligner_type = "spliced" if rna_aligner in RNA_SPLICED_ALIGNERS else "pseudo"
    script:
        "scripts/rna_deseq2_analysis.R"


rule edgeR:
    """
    Executes differential gene/transcript expression analysis using edgeR.
    """
    input:
        counts = get_expression_input,
        runinfo = lambda w: f"{READS_DIR}/{config.get('DATASOURCE', 'SRA').lower()}_runinfo.csv"
    output:
        csv = ( 
        f"{READS_DIR}/expression/edger_expression.{rna_aligner}_counts.csv" 
        if rna_aligner in RNA_SPLICED_ALIGNERS else 
        f"{READS_DIR}/expression/edger_expression.{rna_aligner}_pseudo_counts.csv"
        )
    log:
        f"{LOG_DIR}/diff_exp/edgeR/exp_result.{config.get('ALIGNER', 'star').lower()}.log"
    conda:
        "../env/rna_diff_exp.yaml"
    threads: 4
    params:
        normz = config.get("EDGER_NORM", "TMM"),
        aligner_type = "spliced" if rna_aligner in RNA_SPLICED_ALIGNERS else "pseudo"
    script:
        "scripts/rna_edgeR_analysis.R"

rule tools_compare:
    """
    Cross-references statistics between DESeq2 and edgeR outputs.
    """
    input:
        deseq = f"{READS_DIR}/expression/deseq2_expression.{rna_aligner}_counts.csv" if rna_aligner in RNA_SPLICED_ALIGNERS else f"{READS_DIR}/expression/deseq2_expression.{rna_aligner}_pseudo_counts.csv",
        edger = f"{READS_DIR}/expression/edger_expression.{rna_aligner}_counts.csv" if rna_aligner in RNA_SPLICED_ALIGNERS else f"{READS_DIR}/expression/edger_expression.{rna_aligner}_pseudo_counts.csv"
    output:
        csv = f"{READS_DIR}/expression/cross_tool_comparison.csv",
        venn = f"{READS_DIR}/expression/plots/diff_exp_venn_overlap.png"
    log:
        f"{LOG_DIR}/diff_exp/tools_compare.log"
    conda:
        "../env/rna_diff_exp.yaml"
    params:
        fdr = config.get("RNA_FDR", 0.05),
        lfc = config.get("RNA_LFC", 1.0),
        bg_color = config.get("BG_COLOR", "transparent").lower()
    script:
        "scripts/rna_tools_compare.R"



rule rna_log_fold_threshold:
    """Filters differential results from the configuration-designated active analyzer."""
    input:
        csv = lambda w: f"{READS_DIR}/expression/{rna_exp_analyser}_expression.{rna_aligner}_counts.csv" if rna_aligner in RNA_SPLICED_ALIGNERS else f"{READS_DIR}/expression/{rna_exp_analyser}_expression.{rna_aligner}_pseudo_counts.csv"
    output:
        sig_csv = f"{READS_DIR}/expression/significant_genes_thresholded.csv",
        volcano = f"{READS_DIR}/expression/plots/volcano_plot.png"
    log:
        f"{LOG_DIR}/diff_exp/threshold_filtering.log"
    conda:
        "../env/rna_diff_exp.yaml"
    params:
        analyser = rna_exp_analyser,
        fdr = config.get("RNA_FDR", 0.05),
        lfc = config.get("RNA_LFC", 1.0),
        bg_color = config.get("BG_COLOR", "transparent").lower()
    script:
        "scripts/rna_log_fold_threshold.R"


rule diagnostic_plots:
    """Generates statistical model diagnostics: MA plots and P-value distribution histograms."""
    input:
        csv = lambda w: f"{READS_DIR}/expression/{rna_exp_analyser}_expression.{rna_aligner}_counts.csv" if rna_aligner in RNA_SPLICED_ALIGNERS else f"{READS_DIR}/expression/{rna_exp_analyser}_expression.{rna_aligner}_pseudo_counts.csv"
    output:
        ma = f"{READS_DIR}/expression/plots/diagnostic_ma_plot.png",
        hist = f"{READS_DIR}/expression/plots/diagnostic_pvalue_histogram.png"
    log:
        f"{LOG_DIR}/diff_exp/diagnostics.log"
    conda:
        "../env/rna_diff_exp.yaml"
    params:
        analyser = rna_exp_analyser,
        bg_color = config.get("BG_COLOR", "transparent").lower()
    script:
        "scripts/rna_diagnostic_plots.R"


rule pca_plot:
    """Generates a principal component analysis plot to evaluate global variance patterns."""
    input:
        counts = get_expression_input,
        runinfo = f"{READS_DIR}/{data_src.lower()}_runinfo.csv"
    output:
        png = f"{READS_DIR}/expression/plots/pca_plot.png"
    log:
        f"{LOG_DIR}/diff_exp/pca.log"
    conda:
        "../env/rna_diff_exp.yaml"
    params:
        analyser = rna_exp_analyser,
        aligner_type = "spliced" if rna_aligner in RNA_SPLICED_ALIGNERS else "pseudo",
        batch = config.get("BATCH_COLUMN", "condition"),
        bg_color = config.get("BG_COLOR", "transparent").lower()
    script:
        "scripts/rna_pca_analysis.R"


rule heatmap:
    """Generates hierarchical clustering heatmaps from significant genes."""
    input:
        counts = get_expression_input,
        runinfo = f"{READS_DIR}/{data_src.lower()}_runinfo.csv",
        sig_genes = f"{READS_DIR}/expression/significant_genes_thresholded.csv"
    output:
        png = f"{READS_DIR}/expression/plots/heatmap.png"
    log:
        f"{LOG_DIR}/diff_exp/heatmap.log"
    conda:
        "../env/rna_diff_exp.yaml"
    params:
        batch = config.get("BATCH_COLUMN", "condition"),
        bg_color = config.get("BG_COLOR", "transparent").lower()
    script:
        "scripts/rna_heatmap.R"


rule functional_enrichment:
    """Performs functional enrichment analysis (GO/KEGG pathways) on thresholded target lists."""
    input:
        sig_genes = f"{READS_DIR}/expression/significant_genes_thresholded.csv"
    output:
        go_csv = f"{READS_DIR}/expression/enrichment_go_results.csv",
        go_dot = f"{READS_DIR}/expression/plots/enrichment_go_dotplot.png"
    log:
        f"{LOG_DIR}/diff_exp/functional_enrichment.log"
    conda:
        "../env/rna_diff_exp.yaml"
    params:
        species_latin = species.lower().replace("_", "").replace(".", "").replace(" ", ""),
        bg_color = config.get("BG_COLOR", "transparent").lower()
    script:
        "scripts/rna_functional_enrichment.R"