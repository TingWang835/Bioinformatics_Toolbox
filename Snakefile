import os
import pandas as pd

# =============================================================================
# Setup and Configuration
# =============================================================================
PRJNAME = config.get("PRJNAME") 
if not PRJNAME:
    raise ValueError("ERROR: Specify PRJNAME via --config")

# Load project-specific configuration from the project directory
configfile: f"reads/{PRJNAME}/config.yaml"

# --- Mandatory Configuration Checks---
mandatory_vars = {
    "REF_SOURCE": "choose from ncbi or ensembl",
    "REFNAME": "Reference Genome Name (e.g., Yeast_S288C)",
    "PRJNUMBER": "Project Number or Identifier"
}

for var, description in mandatory_vars.items():
    if not config.get(var):
        raise ValueError(
            f"\n\n[CONFIGURATION ERROR]\n"
            f"Variable '{var}' is missing!\n"
            f"Description: {description}\n"
            f"Please add it to 'config/config.yaml' or provide it via: --config {var}=VALUE\n"
        )

# =============================================================================
# Unified Global Path & Context Variables
# =============================================================================
# Path Roots
ref_src = config.get("REF_SOURCE", "ncbi").lower()
data_src = config.get("DATA_SOURCE", "SRA")
base_src = config.get("DATABASE_SOURCE", "snpeff").lower()

REFS_DIR = f"refs/{config['REFNAME']}/{ref_src}"
READS_DIR = f"reads/{PRJNAME}"
LOG_DIR = f"reads/{PRJNAME}/logs"
DATABASES_DIR = f"database/{base_src}/{config['REFNAME']}"

# aligner variables for path construction
dna_aligner = config.get("DNA_ALIGNER", "bwa").lower()
rna_aligner = config.get("RNA_ALIGNER", "star").lower()

# Tool Wildcard Constraint Configurations
DNA_ALIGNERS = ["bwa", "bowtie2", "minimap2_sr"]
RNA_SPLICED_ALIGNERS = ["star", "hisat2"]
RNA_PSEUDO_ALIGNERS = ["kallisto", "salmon"]

# Ensembl specific variables
species     = config.get("SPECIES_LATIN", "saccharomyces_cerevisiae").lower().replace(" ", "_")
species_cap = species.capitalize()
assembly    = config.get("ASSEMBLY", "R64-1-1")
release     = config.get("RELEASE", 112)

# NCBI specific variables
acc = config.get("ACC", "GCF_000146045.2")


# RNAseq expression analyser
rna_exp_analyser = config.get("EXP_ANALYSER", "deseq2").lower()


# =============================================================================
# Include Modular Rule Files
# =============================================================================
include: "toolbox/getdata.smk" 
include: "toolbox/qc.smk" 
include: "toolbox/dna_aligner.smk" 
include: "toolbox/dna_vcf.smk" 
include: "toolbox/cleanup.smk"
include: "toolbox/dna_rigidity.smk"
include: "toolbox/rna_aligner.smk"
include: "toolbox/rna_diff_exp.smk"

# =============================================================================
# Helper Functions
# =============================================================================

def get_runinfo(wildcards):
    """
    Triggers checkpoints based on DATASOURCE from config.yaml (SRA or LOCAL)
    and resolves the correct runinfo file path dynamically.
    """
    dsource = config.get("DATASOURCE", "SRA").upper() 
    
    if dsource == "SRA":
        cp_output = checkpoints.getdata_sra.get(**wildcards).output.csv 
    elif dsource == "LOCAL":
        cp_output = checkpoints.getdata_local.get(**wildcards).output.csv 
    else:
        raise ValueError(f"Invalid DATASOURCE '{dsource}' in config. Must be SRA or LOCAL.")

    df = pd.read_csv(cp_output) 
    return df['Run'].tolist()

def get_refs(wildcards):
    """Returns a list of reference files needed for alignment, tracking source-specific folder trees."""
    source = config.get("REF_SOURCE", "ncbi").lower()
    
    if source == "ensembl":
        asm = config.get("ASSEMBLY", "R64-1-1")
        rel = config.get("RELEASE", 112)
        return [
            f"{REFS_DIR}/{species_cap}.{asm}.dna.toplevel.fa",
            f"{REFS_DIR}/{species_cap}.{asm}.{rel}.gtf",
            f"{REFS_DIR}/{species_cap}.{asm}.{rel}.gff3"
        ]
    else:
        acc = config.get("ACC", "GCF_000146045.2")
        return [
            f"{REFS_DIR}/{acc}.fa",
            f"{REFS_DIR}/{acc}.gtf",
            f"{REFS_DIR}/{acc}.gff"
        ]

def get_qc(wildcards):
    """Generates QC report targets, merge by multiqc and forces fastq trimming."""
    samples = get_runinfo(wildcards)
    
    fastqc = expand("{rddir}/qc/{s}_{pair}_fastqc.html", 
                    rddir=READS_DIR, s=samples, pair=[1, 2])
    multiqc = [f"{READS_DIR}/qc/multiqc_report.html"]
    
    trim_r1 = expand("{rddir}/qc_trimmed/{s}_1_val_1.fq.gz", rddir=READS_DIR, s=samples)
    trim_r2 = expand("{rddir}/qc_trimmed/{s}_2_val_2.fq.gz", rddir=READS_DIR, s=samples)
  
    return multiqc + trim_r1 + trim_r2

# =============================================================================
# DNA and VCF Helper Functions
# =============================================================================

def get_dna_align(wildcards):
    """Align, generate bam and filtered bam."""
    samples = get_runinfo(wildcards)
    aligner = dna_aligner
    return expand("{rddir}/bam/filtered/{s}.{aln}.filtered.bam",
                  rddir=READS_DIR, s=samples, aln=aligner)

def get_dna_vcf(wildcards):
    """Generates vcf.gz, merged.vcf.gz, ann.vcf."""
    samples = get_runinfo(wildcards)
    aligner = dna_aligner 
    
    annvcf = f"{READS_DIR}/vcf/all_samples.{aligner}.ann.vcf.gz"
    anntbi = f"{READS_DIR}/vcf/all_samples.{aligner}.ann.vcf.gz.tbi"
    annstats = f"{READS_DIR}/vcf/all_samples.{aligner}.snpeff_stats.html"
    consensus = expand("{rddir}/vcf/consensus/{s}.{aln}.consensus.fa", 
                       rddir=READS_DIR, s=samples, aln=aligner)
  
    return [annvcf, anntbi, annstats] + consensus

# =============================================================================
# DNA Ridigity Helper Functions
# =============================================================================
def get_dna_rigid(wildcards): 
    """Generates rigidity profile tracks out of structural data layouts."""
    acc = config.get("ACC", "GCF_000146045.2")
    tsv = f"{READS_DIR}/dna_rigidity/{acc}_rigidity_profile.tsv"
    bg = f"{READS_DIR}/dna_rigidity/{acc}_rigidity_profile.bedGraph"
    return [tsv, bg]

# =============================================================================
# RNAseq Helper Functions
# =============================================================================
def get_rna_align(wildcards):
    """Determines exactly which files to build based on the RNA ALIGNER in config."""
    samples = get_runinfo(wildcards)
    aligner = rna_aligner
    
    if aligner in RNA_PSEUDO_ALIGNERS:
        return expand("{rddir}/counts/{s}.{aln}",
                      rddir=READS_DIR, s=samples, aln=aligner)
    else:
        return expand("{rddir}/bam/{s}.{aln}.bam",
                      rddir=READS_DIR, s=samples, aln=aligner)

def get_rna_bigwig(wildcards):
    """Generates BigWig coverage targets for all samples processed with a spliced aligner."""
    samples = get_runinfo(wildcards)
    aligner = rna_aligner
    
    if aligner in RNA_SPLICED_ALIGNERS:
        return expand("{rddir}/bam/bigwig/{s}.{aln}.bw",
                      rddir=READS_DIR, s=samples, aln=aligner)
    return [] # Returns empty list for pseudo-aligners so it safely skips BigWig generation

def get_rna_counts(wildcards):
    """Compiles the single final merged CSV target path for both alignment styles."""
    aligner = rna_aligner
    
    if aligner in RNA_SPLICED_ALIGNERS:
        return [f"{READS_DIR}/counts/all_samples.{aligner}_counts.csv"]
    else:
        # Automatically routes to your new merge_pseudo matrix for Kallisto/Salmon
        return [f"{READS_DIR}/counts/all_samples.{aligner}_pseudo_counts.csv"]

def get_rna_exp_analyser_output(wildcards):
    """Determines the RNA expression target path based on EXPRESSION_ANALYSER in config."""
    suffix = "counts.csv" if rna_aligner in RNA_SPLICED_ALIGNERS else "pseudo_counts.csv"
    return f"{READS_DIR}/expression/{rna_exp_analyser}_expression.{rna_aligner}_{suffix}"

def get_rna_tools_compare(wildcards):
    """Determines the comprehensive final target path array for the RNA-seq module."""
    return [
        f"{READS_DIR}/expression/cross_tool_comparison.csv",
        f"{READS_DIR}/expression/plots/diff_exp_venn_overlap.png"
    ]

def get_rna_report_outputs(wildcards):
    """Gathers all biological thresholds and quality control diagnostic plots."""
    return [
        f"{READS_DIR}/expression/significant_genes_thresholded.csv",
        f"{READS_DIR}/expression/plots/volcano_plot.png",
        f"{READS_DIR}/expression/plots/diagnostic_ma_plot.png",
        f"{READS_DIR}/expression/plots/diagnostic_pvalue_histogram.png",
        f"{READS_DIR}/expression/plots/pca_plot.png"
    ]

def get_rna_enrich_outputs(wildcards):
    """Generate heatmap and function enrichment study results."""
    return [
        f"{READS_DIR}/expression/plots/heatmap.png",
        f"{READS_DIR}/expression/enrichment_go_results.csv",
        f"{READS_DIR}/expression/plots/enrichment_go_dotplot.png",
    ]


# =============================================================================
# Terminal Rules
# =============================================================================
rule note:
    input:
    run:
        print("\n" + "="*50)
        print("THIS SNAKEFILE SERVES AS TERMINAL FOR BIOINFORMATIC TOOLBOX")
        print("Please specify a target rule:")
        print("Data/Reference Gathering")
        print("  snakemake runinfo           - Get runinfo for samples")
        print("  snakemake download_refs     - Download fa, gff and gtf files")
        print("  snakemake qc                - Download data and run QC")
        print("SNP/Variant Analysis")
        print("  snakemake dna_align         - Run DNA alignments")
        print("  snakemake dna_vcf           - Run DNA VCF")
        print("  snakemake dna_vcf_all       - Run QC, DNA Alignments and VCF")
        print("DNA Rigidity Scoring")
        print("  snakemake dna_rigid         - Run dna rigidity score on fasta")
        print("RNAseq Analysis")
        print("  snakemake rna_align         - Run RNA alignment")
        print("  snakemake rna_exp           - Run RNA expression analysis")
        print("  snakemake exp_tools_compare - (optional) Compare expression analysis tools deseq2 vs edger")
        print("  snakemake rna_report        - Evaluate expression data with pca, MA, p-value histagram and volcano plots,")
        print("                                generate csv for list of significantly different expressed genes")
        print("  snakemake rna_enrich        - Generate heatmap and perform functional enrichment study")
        print("  snakemake rna_all           - Gets samples and refs, run QC, rna align, exp, report and enrich")
        print("Clean up")
        print("  snakemake cleanup           - Clean up dummy R2 files from single end sequencing")
        print("="*50 + "\n")

rule runinfo:
    """Download runinfo.csv or create one from fastq files under reads/project_name"""
    input: get_runinfo

rule download_refs:
    """Utility run-target that safely maps upstream conversion assets."""
    input: get_refs

rule qc:
    """QC for fastq files."""
    input: get_qc

# =============================================================================
# DNA alignment and VCF
# =============================================================================
rule dna_align:
    """Align fastq files with index by chosen aligner in config.yaml"""
    input: get_dna_align

rule dna_vcf:
    """Run vcf call, merge samples, create consensus, build snpeff db, annotate."""
    input: get_dna_vcf

rule dna_vcf_all:
    """A shortcut to run a full course with qc, aligner and vcf."""
    input:  
        get_qc,
        get_dna_align,
        get_dna_vcf

# =============================================================================
# DNA Ridigity Scoring
# =============================================================================
rule dna_rigid: 
    """Run dna rigidity score on fasta."""
    input: get_dna_rigid

# =============================================================================
# RNA Alignment and Expression Analysis
# =============================================================================
rule rna_align:
    """Align RNAseq reads using spliced aligner/ pseudo aligner and generate visualization tracks."""
    input:
        alignments = get_rna_align,
        bigwigs    = get_rna_bigwig,
        counts     = get_rna_counts
        
rule rna_exp:
    """Perform expression analysis using DESeq2 or edgeR."""
    input:
        get_rna_exp_analyser_output

rule exp_tools_compare:
    """Compare statistics between DESeq2 and edgeR outputs."""
    input:
        get_rna_tools_compare

rule rna_report:
    """
    Executes biological threshold filtering,model diagnostics, and global PCA variance plotting.
    """
    input:
        get_rna_report_outputs


rule rna_enrich:
    """
    Generate Heatmap and perform functional enrichment study.
    """
    input:
        get_rna_enrich_outputs

rule rna_all:
    """A shortcut to run a full course of RNAseq with qc, aligner, exp, report and enrich."""
    input:  
        get_rna_align,
        get_rna_bigwig,
        get_rna_counts,
        get_rna_exp_analyser_output,
        get_rna_report_outputs,
        get_rna_enrich_outputs




# =============================================================================
# Clean up
# =============================================================================
rule cleanup:
    """Cleaning up all dummy R2 files created from single end sequencing."""
    input: f"{READS_DIR}/logs/cleanup_complete.done"