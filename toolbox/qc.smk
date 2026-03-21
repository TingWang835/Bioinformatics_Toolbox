rule fastqc:
    input:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq"
    output:
        html1 = f"reads/{PRJNAME}/qc/{{sample}}_1_fastqc.html",
        html2 = f"reads/{PRJNAME}/qc/{{sample}}_2_fastqc.html"
    log: f"reads/{PRJNAME}/logs/fastqc/{{sample}}.log"
    conda: "../env/qc.yaml"
    shell:
        """
        if [ -s "{input.r2}" ]; then
            fastqc {input.r1} {input.r2} --outdir reads/{PRJNAME}/qc/ > {log} 2>&1
        else
            fastqc {input.r1} --outdir reads/{PRJNAME}/qc/ > {log} 2>&1
            touch {output.html2} # Create a dummy output to satisfy Snakemake
        fi
        """

rule trim_reads:
    """
    Handles quality trimming (tail cutting) at Phred 30.
    """
    input:
        r1 = f"reads/{PRJNAME}/{{sample}}_1.fastq",
        r2 = f"reads/{PRJNAME}/{{sample}}_2.fastq"
    output:
        # Trim Galore naming convention for paired: _val_1 and _val_2
        r1_p = f"reads/{PRJNAME}/qc_trimmed/{{sample}}_1_val_1.fq.gz",
        r2_p = f"reads/{PRJNAME}/qc_trimmed/{{sample}}_2_val_2.fq.gz"
    log: f"reads/{PRJNAME}/logs/trimming/{{sample}}.log"
    conda: "../env/qc.yaml"
    shell:
        """
        if [ -s "{input.r2}" ]; then
            trim_galore --paired --gzip -q 30 --output_dir reads/{PRJNAME}/qc_trimmed/ {input.r1} {input.r2} > {log} 2>&1
        else
            trim_galore --gzip -q 30 --output_dir reads/{PRJNAME}/qc_trimmed/ {input.r1} > {log} 2>&1
            touch {output.r2_p}
        fi
        """