import os
import pandas as pd

# 1. Setup and Configuration
PRJNAME = config.get("PRJNAME") 
if not PRJNAME:
    raise ValueError("ERROR: Specify PRJNAME via --config")

# Load project-specific configuration from the project directory
configfile: f"reads/{PRJNAME}/config.yaml"

# 2. Include Modular Rule Files
include: "toolbox/getdata.smk" 
include: "toolbox/qc.smk" 
include: "toolbox/aligner.smk" 
include: "toolbox/vcf.smk" 
include: "toolbox/cleanup.smk"
include: "toolbox/dna_rigidity.smk"
# include: "toolbox/rnaseq.smk" 

# 3. Get Functions
def get_samples(wildcards):
    """
    Trigger checkpoints base on dsource from config.yaml (default: SRA)
    local_runinfo.csv will be generated based on fastq file names under read/[project_name]
    """
    dsource = config.get("DSOURCE", "SRA").upper()
    
    # Map data source to the corresponding checkpoint
    if dsource == "SRA":
        cp_output = checkpoints.getdata_sra.get(**wildcards).output.csv 
    else:
        cp_output = checkpoints.getdata_local.get(**wildcards).output.csv 
        
    df = pd.read_csv(cp_output) 
    return df['Run'].tolist() #make list using "Run" (SRR number) column


def get_qc(wildcards):
    """
    Generates QC report targets, merge by multiqc and forces fastq trimming.
    """
    samples = get_samples(wildcards)
    
    
    # Standard FastQC reports
    fastqc = expand("reads/{prj}/qc/{s}_{pair}_fastqc.html", 
                    prj=PRJNAME, s=samples, pair=[1, 2])
    # Summarise using multiqc
    multiqc = f"reads/{PRJNAME}/qc/multiqc_report.html"
    # Trimmed fastq files (Note: Trim Galore naming uses _1_val_1 and _2_val_2)
    trim_r1 = expand("reads/{prj}/qc_trimmed/{s}_1_val_1.fq.gz", prj=PRJNAME, s=samples)
    trim_r2 = expand("reads/{prj}/qc_trimmed/{s}_2_val_2.fq.gz", prj=PRJNAME, s=samples)
    
    # Return a combined list of all targets
    return [multiqc] + fastqc + trim_r1 + trim_r2


def get_bam(wildcards):
    """
    Align, generate bam and filtered bam.
    """
    samples = get_samples(wildcards)
    aligner = config.get("ALIGNER", "bwa")
    return expand("reads/{prj}/bam/filtered/{s}.{aln}.filtered.bam",
                  prj=PRJNAME, s=samples, aln=aligner)


def get_vcf(wildcards):
    """
    Generates vcf.gz, merged.vcf.gz, ann.vcf 
    """
    samples = get_samples(wildcards)
    aligner = config.get("ALIGNER", "bwa") 
    annvcf =  f"reads/{PRJNAME}/vcf/all_samples.{aligner}.ann.vcf.gz"
    anntbi =  f"reads/{PRJNAME}/vcf/all_samples.{aligner}.ann.vcf.gz.tbi"
    annstats = f"reads/{PRJNAME}/vcf/all_samples.{aligner}.snpeff_stats.html"
    consensus = expand("reads/{prj}/vcf/consensus/{s}.{aln}.consensus.fa", prj=PRJNAME, s=samples, aln=aligner)
    return [annvcf, anntbi, annstats] + consensus
    
def get_rigid(wildcards):
    """
    Generates ACC_rigidity_profile.tsv and ACC_rigidity_profile.bedgraph
    """ 
    tsv = f"reads/{PRJNAME}/dna_rigidity/{config['ACC']}_rigidity_profile.tsv"
    bg = f"reads/{PRJNAME}/dna_rigidity/{config['ACC']}_rigidity_profile.bedGraph"
    return [tsv, bg]

def get_rnaseq(wildcards):
    """
    placeholder.
    """
    samples = get_samples(wildcards)
    return expand("reads/{prj}/counts/{s}.counts.txt", 
                  prj=PRJNAME, s=samples)



# 4. Terminal Rules
rule note:
    input:
    run:
        print("\n" + "="*50)
        print("THIS SNAKEFILE SERVICE AS TERMINAL FOR BIOINFOMATIC TOOLBOX")
        print("Please specify a target rule:")
        print("  snakemake qc      - Run FastQC")
        print("  snakemake bam     - Run Alignments")
        print("  snakemake vcf     - Run VCF")
        print("  snakemake vcf_all - Run FastQC, Alignments, VCF")
        print("  snakemake rigid   - Run dna rigidity score on fasta")
        print("  snakemake rnaseq  - Run RNA-seq analysis (place holder)")
        print("  snakemake cleanup  - clean up dummy R2 files from single end sequencing")
        print("="*50 + "\n")

rule qc:
    """
    QC for fastq files.
    """
    input: get_qc


rule bam:
    """
    Align fastq files with index by chosen aligner in config.yaml
    """
    input: get_bam

rule vcf:
    """
    Run vcf call, merge samples, create consensus, build snpeff db, annotate.
    """
    input: get_vcf

rule vcf_all:
    """
    A short cut to run a full course with qc, aligner and vcf.
    """
    input:  
        lambda wildcards: get_qc(wildcards),
        lambda wildcards: get_bam(wildcards),
        lambda wildcards: get_vcf(wildcards)

rule rigid: 
    """
    run dna rigidity score on fasta
    """
    input: get_rigid

rule rnaseq:
    """
    placeholder.
    """
    input: get_rnaseq



rule cleanup:
    """
    Cleaning up all dummy R2 files created from single end sequencing
    """
    input: f"reads/{PRJNAME}/logs/cleanup_complete.done"




