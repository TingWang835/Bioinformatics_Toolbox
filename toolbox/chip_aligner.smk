# =============================================================================
# ChIP-Seq Alignment Module
# Note for future development / DAG troubleshooting:
# This module completely reuses the reference index files and the 
# 'get_index_prefix' helper function defined inside 'dna_aligner.smk'.
# =============================================================================

# Constrain the aligner wildcard strictly to prevent string parsing bleed
wildcard_constraints:
    aligner = "bowtie2|bwa"

# =============================================================================
# Alignment & Filtering Rules
# =============================================================================

rule dna_chip_align:
    input:
        unpack(get_align_input),
        idx = get_index_prefix # Reused from dna_aligner.smk
    output:
        bam = temp(f"{READS_DIR}/bam/{{sample}}.{{aligner}}.chip.bam")
    log:
        f"{LOG_DIR}/aligner/align/{{sample}}.{{aligner}}.log"
    conda:
        "../env/dna_aligner.yaml"
    threads: 4
    params:
        cmd = lambda w, threads: (
            f"bwa mem -t {threads} {get_index_prefix(w).replace('.amb', '')}" if w.aligner == "bwa" else
            f"bowtie2 --threads {threads} --local -X 1000 --no-mixed --no-discordant -x {get_index_prefix(w).replace('.1.bt2', '')}"
        ),
        pe = lambda w, input: (
            f"{input.r1} {input.r2}" if w.aligner == "bwa" else
            f"-1 {input.r1} -2 {input.r2}"
        ),
        se = lambda w, input: (
            f"{input.r1}" if w.aligner == "bwa" else
            f"-U {input.r1}"
        )
    shell:
        """
        if [ -s "{input.r2}" ]; then
            {params.cmd} {params.pe} 2> {log} | samtools sort -@ {threads} -o {output.bam}
        else
            {params.cmd} {params.se} 2> {log} | samtools sort -@ {threads} -o {output.bam}
        fi
        """

rule chip_filter_bam:
    input:
        bam = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.chip.bam"
    output:
        filtered = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.chip.filtered.bam",
        bai = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.chip.filtered.bam.bai"
    log:
        f"{LOG_DIR}/aligner/filter_bam/{{sample}}.{{aligner}}.log"
    conda:
        "../env/dna_aligner.yaml"
    threads: 4
    params:
        mapq = config.get("MAPQ", "30")
    shell:
        """
        samtools view -b -q {params.mapq} {input.bam} > {output.filtered} 2> {log}
        samtools index -@ {threads} {output.filtered} >> {log} 2>&1
        """

rule chip_mark_dup:
    input:
        bam = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.chip.filtered.bam"
    output:
        bam = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.chip.ready.bam",
        bai = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.chip.ready.bam.bai",
        metrics = f"{LOG_DIR}/aligner/remove_dup/{{sample}}.{{aligner}}.metrics.txt"
    log:
        f"{LOG_DIR}/aligner/remove_dup/{{sample}}.{{aligner}}.log"
    conda:
        "../env/dna_aligner.yaml"
    threads: 2
    shell:
        """
        # Inject required Read Group tag using tabs (\t) to keep parameters safely bound
        samtools addreplacerg -r "ID:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" -m orphan_only -o - {input.bam} 2> {log} | \
        picard MarkDuplicates \
            I=/dev/stdin \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            ASSUME_SORTED=true \
            READ_NAME_REGEX=null \
            VALIDATION_STRINGENCY=SILENT >> {log} 2>&1

        # Index the final ready BAM file
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """

rule chip_bam_to_bigwig:
    """
    Converts filtered bam into normalized TPM BigWig coverage tracks 
    using deepTools bamCoverage for optimized loading and row-scaling inside IGV.
    """
    input:
        bam = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.chip.filtered.bam",
        bai = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.chip.filtered.bam.bai"
    output:
        bw = f"{READS_DIR}/bam/bigwig/{{sample}}.{{aligner}}.bw"
    log:
        f"{LOG_DIR}/aligner/bigwig/{{sample}}.{{aligner}}.log"
    conda:
        "../env/dna_aligner.yaml"
    threads: 4
    params:
        normalize = "--normalizeUsing BPM"
    shell:
        """
        bamCoverage --bam {input.bam} \
                    --outFileName {output.bw} \
                    --numberOfProcessors {threads} \
                    {params.normalize} > {log} 2>&1
        """