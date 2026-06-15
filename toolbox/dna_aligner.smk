localrules: bwa_index, bowtie2_index, minimap2_index

# =============================================================================
# Helper Functions
# =============================================================================

def get_tag(wildcards):
    """Generates the Read Group tag for BWA/Minimap2 alignment metrics."""
    return r"@RG\tID:{s}\tSM:{s}\tLB:{s}".format(s=wildcards.sample)

def get_index_prefix(wildcards):
    """Retrieves the reference index prefix rooted cleanly inside the reference directory."""
    if wildcards.aligner == "bwa":
        return f"{REFS_DIR}/bwa/index.amb" 
    elif wildcards.aligner == "bowtie2":
        return f"{REFS_DIR}/bowtie2/index.1.bt2"
    elif wildcards.aligner == "minimap2_sr":
        return f"{REFS_DIR}/minimap2.mmi"
    return f"{REFS_DIR}/{wildcards.aligner}/index"

# =============================================================================
# Index Rules
# =============================================================================

rule bwa_index:
    input: 
        unpack(get_refs)
    output: 
        multiext(f"{REFS_DIR}/bwa/index", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log: 
        f"{LOG_DIR}/aligner/bwa_index.log"
    conda: 
        "../env/dna_aligner.yaml"
    params: 
        prefix = f"{REFS_DIR}/bwa/index"
    shell: 
        "bwa index {input.fasta} -p {params.prefix} > {log} 2>&1"

rule bowtie2_index:
    input: 
        unpack(get_refs)
    output: 
        multiext(f"{REFS_DIR}/bowtie2/index", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    log: 
        f"{LOG_DIR}/aligner/bowtie2_index.log"
    conda: 
        "../env/dna_aligner.yaml"
    params: 
        prefix = f"{REFS_DIR}/bowtie2/index"
    shell: 
        "bowtie2 build {input.fasta} {params.prefix} > {log} 2>&1"

rule minimap2_index:
    input: 
        unpack(get_refs)
    output: 
        f"{REFS_DIR}/minimap2.mmi"
    log: 
        f"{LOG_DIR}/aligner/minimap2_index.log"
    conda: 
        "../env/dna_aligner.yaml"
    shell: 
        "minimap2 -d {output} {input.fasta} > {log} 2>&1"

# =============================================================================
# Aligner Rules
# =============================================================================

rule align:
    input:
        unpack(get_align_input),
        idx_check = get_index_prefix
    output:
        bam = temp(f"{READS_DIR}/bam/{{sample}}.{{aligner}}.raw.bam"),
        bai = temp(f"{READS_DIR}/bam/{{sample}}.{{aligner}}.raw.bam.bai")
    log: 
        f"{LOG_DIR}/aligner/align/{{sample}}.{{aligner}}.log"
    conda: 
        "../env/dna_aligner.yaml"
    threads: 4
    params:
        tag = get_tag,
        cmd = lambda w: (
            f"bwa mem -R '{get_tag(w)}' {get_index_prefix(w).replace('.amb', '')}" if w.aligner == "bwa" else
            f"bowtie2 --rg-id '{w.sample}' --rg '{w.sample}' -x {get_index_prefix(w).replace('.1.bt2', '')}" if w.aligner == "bowtie2" else
            f"minimap2 -ax sr -R '{get_tag(w)}' {get_index_prefix(w)}"
        ),
        pe = lambda w, input: (
            f"{input.r1} {input.r2}" if w.aligner in ["bwa", "minimap2_sr"] else
            f"-1 {input.r1} -2 {input.r2}"
        ),
        se = lambda w, input: (
            f"{input.r1}" if w.aligner in ["bwa", "minimap2_sr"] else
            f"-U {input.r1}"
        )
    shell:
        """
        if [ -s "{input.r2}" ]; then
            {params.cmd} {params.pe} | samtools sort -@ {threads} -o {output.bam}
        else
            {params.cmd} {params.se} | samtools sort -@ {threads} -o {output.bam}
        fi
        samtools index {output.bam}
        """

rule filter_bam:
    input:
        bam = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.raw.bam",
        bai = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.raw.bam.bai"
    output:
        filtered = f"{READS_DIR}/bam/filtered/{{sample}}.{{aligner}}.filtered.bam"
    log: 
        f"{LOG_DIR}/aligner/filter_bam/{{sample}}.{{aligner}}.log"
    conda: 
        "../env/dna_aligner.yaml"
    params:
        mapq = config.get("MAPQ", "30"),
        samflag = config.get("SAMFLAG", "3")
    shell:
        """
        samtools view -b -q {params.mapq} -f {params.samflag} {input.bam} > {output.filtered}
        samtools index {output.filtered}
        """