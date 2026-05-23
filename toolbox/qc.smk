localrules: fastqc, multiqc, trim_reads

# =============================================================================
# Rules
# =============================================================================

rule fastqc:
    input:
        r1 = f"{READS_DIR}/{{sample}}_1.fastq",
        r2 = f"{READS_DIR}/{{sample}}_2.fastq"
    output:
        html1 = f"{READS_DIR}/qc/{{sample}}_1_fastqc.html",
        html2 = f"{READS_DIR}/qc/{{sample}}_2_fastqc.html",
        zip1  = f"{READS_DIR}/qc/{{sample}}_1_fastqc.zip",
        zip2  = f"{READS_DIR}/qc/{{sample}}_2_fastqc.zip"
    log: f"{LOG_DIR}/qc/fastqc/{{sample}}.log"
    conda: "../env/qc.yaml"
    shell:
        """
        # -s checks if the file exists and has size > 0
        if [ -s "{input.r2}" ]; then
            fastqc {input.r1} {input.r2} --outdir {READS_DIR}/qc/ > {log} 2>&1
        else
            # Single-end: only run R1
            fastqc {input.r1} --outdir {READS_DIR}/qc/ > {log} 2>&1
            
            # Create dummy outputs for R2 to satisfy Snakemake
            touch {output.html2} {output.zip2}
            echo "R2 is empty dummy, FastQC skipped R2." >> {log}
        fi
        """

rule multiqc:
    input:
        lambda wildcards: expand("{rddir}/qc/{s}_{pair}_fastqc.html", 
                                 rddir=READS_DIR, s=get_runinfo(wildcards), pair=[1, 2])
    output:
        f"{READS_DIR}/qc/multiqc_report.html"
    params:
        qc_dir = f"{READS_DIR}/qc/"
    conda: "../env/multiqc.yaml"
    shell:
        # --force overwrites old reports so the program won't stop if one exists
        "multiqc {params.qc_dir} -o {params.qc_dir} -n multiqc_report.html --force"

rule trim_reads:
    input:
        r1 = f"{READS_DIR}/{{sample}}_1.fastq",
        r2 = f"{READS_DIR}/{{sample}}_2.fastq"
    output:
        r1_p = f"{READS_DIR}/qc_trimmed/{{sample}}_1_val_1.fq.gz",
        r2_p = f"{READS_DIR}/qc_trimmed/{{sample}}_2_val_2.fq.gz"
    params:
        quality = config.get("QUALITY", "20"),
        length = config.get("LENGTH", "20")
    log: f"{LOG_DIR}/qc/trim_reads/{{sample}}.log"
    conda: "../env/qc.yaml"
    shell:
        """
        if [ -s "{input.r2}" ]; then
            # Real paired-end data
            trim_galore --paired --gzip -q {params.quality} --length {params.length} --output_dir {READS_DIR}/qc_trimmed/ {input.r1} {input.r2} > {log} 2>&1
        else
            # Single-end data
            trim_galore --gzip -q {params.quality} --length {params.length} --output_dir {READS_DIR}/qc_trimmed/ {input.r1} > {log} 2>&1
            
            # 1. Rename SE output.
            # Trim Galore SE output for sample_1.fastq is sample_1_trimmed.fq.gz
            if [ -f "{READS_DIR}/qc_trimmed/{wildcards.sample}_1_trimmed.fq.gz" ]; then
                mv {READS_DIR}/qc_trimmed/{wildcards.sample}_1_trimmed.fq.gz {output.r1_p}
            elif [ -f "{READS_DIR}/qc_trimmed/{wildcards.sample}_1_trimmed.fq" ]; then
                gzip -c {READS_DIR}/qc_trimmed/{wildcards.sample}_1_trimmed.fq > {output.r1_p}
                rm {READS_DIR}/qc_trimmed/{wildcards.sample}_1_trimmed.fq
            fi
            
            # 2. Create a dummy gzipped file for R2
            touch {output.r2_p}
            echo "Single-end mode: Processed R1 and created dummy R2.gz" >> {log}
        fi
        """