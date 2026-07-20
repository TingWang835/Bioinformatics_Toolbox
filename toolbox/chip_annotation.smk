wildcard_constraints:
    aligner = "bowtie2|bwa"
dsource = config.get("DATASOURCE", "SRA").lower()
runinfo_df = pd.read_csv(f"{READS_DIR}/{dsource}_runinfo.csv")
# =============================================================================
# Helper function
# =============================================================================
def get_peak_bed(wildcards):
    """assign path based on peakcaller"""
    if peakcaller == "macs3":
        return {"peak_bed": f"{READS_DIR}/peaks/macs3/{wildcards.sample}_filtered.bed"}
    elif peakcaller == "gem":
        return {"peak_bed": f"{READS_DIR}/peaks/gem/{wildcards.sample}/{wildcards.sample}_GEM_events.bed"}
    else:
        raise ValueError(f"Unsupported CALLER target: {peakcaller}")

def get_chip_matrix_input(wildcards):
    """
    Looks up the treatment sample directly from the wildcard, 
    and finds its corresponding control matching the same group.
    """
    current_sample_row = runinfo_df[runinfo_df["Run"] == wildcards.sample]
    if current_sample_row.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in runinfo.csv")
        
    current_group = current_sample_row["group"].values[0]
    ctrl_cond = config.get("CONTROL_COND", "control_a")
    aligner = config.get("DNA_ALIGNER", "bwa")

    control_match = runinfo_df[
        (runinfo_df["group"] == current_group) & 
        (runinfo_df["condition"] == ctrl_cond)
    ]
    if control_match.empty:
        raise ValueError(f"No matching control found for group {current_group} under condition '{ctrl_cond}'")
        
    control_srr = control_match["Run"].values[0]
    
    return {
        "treatment": f"{READS_DIR}/bam/bigwig/{wildcards.sample}.{aligner}.bw",
        "control": f"{READS_DIR}/bam/bigwig/{control_srr}.{aligner}.bw"
    }

# =============================================================================
# Orgdb mapping
# =============================================================================
ORGDB_MAPPING = {
    "saccharomyces_cerevisiae": "org.Sc.sgd.db",
    "homo_sapiens": "org.Hs.eg.db"
}
chip_orgdb = lambda wildcards: ORGDB_MAPPING.get(config["REF"]["SPECIES"].lower(), "org.Sc.sgd.db")

# =============================================================================
# Riles
# =============================================================================
rule peak_assignment:
    input:
        unpack(get_peak_bed),
        unpack(get_refs)
    output:
        annotation_csv = f"{READS_DIR}/annotation/{{sample}}_annotated_peaks.csv",
        plot_distribution = f"{READS_DIR}/annotation/plots/{{sample}}_genomic_distribution.png"
    params:
        orgdb = chip_orgdb,
        tss_region_start = config.get("TSS_REGION_START", -3000),
        tss_region_end = config.get("TSS_REGION_end", 100)
    log:
        f"{LOG_DIR}/annotation/peak_assignment/{{sample}}.log"
    conda:
        "../env/chip_annotation.yaml"
    script:
        "scripts/chip_peak_assignment.R"



rule tss_intensity_matrix:
    """
    Calculates signal density windows around TSS coordinates using the clean Q30 alignments
    for a specific treatment sample paired directly with its specific group control.
    """
    input:
        unpack(get_chip_matrix_input),
        unpack(get_refs)
    output:
        matrix = f"{READS_DIR}/annotation/{{sample}}_vs_control_tss_matrix.gz"
    params:
        upstream = 3000, 
        downstream = 1000
    log: f"{LOG_DIR}/annotation/matrix/{{sample}}.log"
    threads: 4
    conda:
        "../env/chip_annotation.yaml"
    shell:
        """
        computeMatrix reference-point \
            --referencePoint TSS \
            -b {params.upstream} \
            -a {params.downstream} \
            -R {input.gtf} \
            -S {input.treatment} {input.control} \
            -o {output.matrix} \
            -p {threads} \
            --missingDataAsZero
        """

rule intensity_profile_plot:
    """
    Renders the summary meta-profile line graph based on the spatial matrix calculation,
    overlaying the treatment and its specific matched control.
    """
    input:
        matrix = f"{READS_DIR}/annotation/{{sample}}_vs_control_tss_matrix.gz"
    output:
        profile_plot = f"{READS_DIR}/annotation/plots/{{sample}}_vs_control_tss_profile.png"
    params:
        labels = lambda wildcards: f'"{wildcards.sample}_{config.get("TREAT_COND", "treatment_a")}" "{config.get("CONTROL_COND", "control_a")}"'
    log: f"{LOG_DIR}/annotation/plot/{{sample}}_profileplot.log"
    conda:
        "../env/chip_annotation.yaml"
    shell:
        """
        plotProfile \
            -m {input.matrix} \
            -o {output.profile_plot} \
            --plotType lines \
            --perGroup \
            --samplesLabel {params.labels} 
        """

rule intensity_heatmap:
    """
    Generates a spatial density heatmap showing side-by-side vertical panels
    for both the treatment sample and its corresponding group control.
    """
    input:
        matrix = f"{READS_DIR}/annotation/{{sample}}_vs_control_tss_matrix.gz"
    output:
        heatmap_plot = f"{READS_DIR}/annotation/plots/{{sample}}_vs_control_density_heatmap.png"
    params:
        labels = lambda wildcards: f'"{wildcards.sample}_{config.get("TREAT_COND", "treatment_a")}" "{config.get("CONTROL_COND", "control_a")}"'
    log: f"{LOG_DIR}/annotation/plot/{{sample}}_heatmap.log"
    conda:
        "../env/chip_annotation.yaml"
    shell:
        """
        plotHeatmap \
            -m {input.matrix} \
            -o {output.heatmap_plot} \
            --plotType lines \
            --colorMap Blues \
            --samplesLabel {params.labels} \
            --whatToShow "plot, heatmap and colorbar"
        """

rule chip_enrichment:
    """
    Performs GO and KEGG functional enrichment analysis using clusterProfiler.
    """
    input:
        annotated_csv = f"{READS_DIR}/annotation/{{sample}}_annotated_peaks.csv"
    output:
        go_results = f"{READS_DIR}/annotation/enrichment/{{sample}}_go_enrichment.csv",
        kegg_results = f"{READS_DIR}/annotation/enrichment/{{sample}}_kegg_enrichment.csv",
        go_dotplot = f"{READS_DIR}/annotation/enrichment/{{sample}}_go_dotplot.png",
        kegg_dotplot = f"{READS_DIR}/annotation/enrichment/{{sample}}_kegg_dotplot.png"
    log:
        f"{LOG_DIR}/enrichment/{{sample}}_enrichment.log"
    params:
        orgdb = chip_orgdb,
        kegg_organism = config.get("KEGG_ORG", "hsa"), 
        p_cutoff = config.get("ENRICH_PVAL", 0.05),
        q_cutoff = config.get("ENRICH_QVAL", 0.2),
        bg_color = config.get("BG_COLOR", "transparent").lower(),
        dotplot_height = config.get("DOTPLOT_HEIGHT", 16)
    conda:
        "../env/chip_annotation.yaml"
    script:
        "scripts/chip_enrichment.R"