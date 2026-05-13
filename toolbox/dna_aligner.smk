localrules:bwa_index, bowtie2_index, minimap2_index, align, filter_bam
# 1. Input Functions
def get_tag(wildcards):
    """Generates the Read Group tag for BWA/Minimap2."""
    return r"@RG\tID:{s}\tSM:{s}\tLB:{s}".format(s=wildcards.sample)

def get_index_prefix(wildcards):
    """Retrieves the reference index prefix and links rules."""
    prefix = f"refs/{config['REFNAME']}/{wildcards.aligner}"
    if wildcards.aligner == "bwa":
        return f"{prefix}.amb" 
    elif wildcards.aligner == "bowtie2":
        return f"{prefix}.1.bt2"
    #elif wildcards.aligner =="bowtie2_large" #a place holder for genomes larger than 4GB
        #return f"{prefix}.1.bt2l"
    elif wildcards.aligner == "minimap2_sr":
        return f"refs/{config['REFNAME']}/minimap2.mmi"
    return prefix

def get_align_input(wildcards):
    """Points to the output of the trimming rule."""
    return {
        "r1": f"reads/{PRJNAME}/qc_trimmed/{wildcards.sample}_1_val_1.fq.gz",
        "r2": f"reads/{PRJNAME}/qc_trimmed/{wildcards.sample}_2_val_2.fq.gz"
    }

# 2. Rules
rule bwa_index:
    input: ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output: multiext(f"refs/{config['REFNAME']}/bwa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log: f"reads/{PRJNAME}/logs/aligner/bwa_index.log"
    conda: "../env/dna_aligner.yaml"
    params: prefix = f"refs/{config['REFNAME']}/bwa"
    shell: "bwa index {input.ref} -p {params.prefix} > {log} 2>&1"
    

rule bowtie2_index:
    input: ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output: multiext(f"refs/{config['REFNAME']}/bowtie2", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    log: f"reads/{PRJNAME}/logs/aligner/bowtie2_index.log"
    conda: "../env/dna_aligner.yaml"
    params: prefix = f"refs/{config['REFNAME']}/bowtie2"
    shell: "bowtie2 build {input.ref} {params.prefix} > {log} 2>&1"

rule minimap2_index:
    input: ref = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output: f"refs/{config['REFNAME']}/minimap2.mmi"
    log: f"reads/{PRJNAME}/logs/aligner/minimap2_index.log"
    conda: "../env/dna_aligner.yaml"
    shell: "minimap2 -d {output} {input.ref} > {log} 2>&1"

rule align:
    input:
        unpack(get_align_input),
        # This line forces the indexers to run
        idx_check = get_index_prefix
    output:
        bam = temp(f"reads/{PRJNAME}/bam/{{sample}}.{{aligner}}.raw.bam"),
        bai = temp(f"reads/{PRJNAME}/bam/{{sample}}.{{aligner}}.raw.bam.bai")
    log: f"reads/{PRJNAME}/logs/aligner/align/{{sample}}.{{aligner}}.log"
    conda: "../env/dna_aligner.yaml"
    params:
        tag = get_tag,
        # Strip the extension for bwa/bowtie2 so they get the prefix
        cmd = lambda w: (
            f"bwa mem -R '{get_tag(w)}' {get_index_prefix(w).replace('.amb', '')}" if w.aligner == "bwa" else
            f"bowtie2 --rg-id '{w.sample}' --rg '{w.sample}' -x {get_index_prefix(w).replace('.1.bt2', '')}" if w.aligner == "bowtie2" else
            f"minimap2 -ax sr -R '{get_tag(w)}' {get_index_prefix(w)}" # minimap2 uses the full file path
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
            {params.cmd} {params.pe} | samtools sort -o {output.bam}
        else
            {params.cmd} {params.se} | samtools sort -o {output.bam}
        fi
        samtools index {output.bam}
        """
        

rule filter_bam:
    input:
        bam = f"reads/{PRJNAME}/bam/{{sample}}.{{aligner}}.raw.bam",
        bai = f"reads/{PRJNAME}/bam/{{sample}}.{{aligner}}.raw.bam.bai"
    output:
        filtered = f"reads/{PRJNAME}/bam/filtered/{{sample}}.{{aligner}}.filtered.bam"
    log: f"reads/{PRJNAME}/logs/aligner/filter_bam/{{sample}}.{{aligner}}.log"
    conda: "../env/dna_aligner.yaml"
    params:
        mapq = config.get("MAPQ", "30"),
        samflag = config.get("SAMFLAG", "3")
    shell:
        """
        samtools view -b -q {params.mapq} -f {params.samflag} {input.bam} > {output.filtered}
        samtools index {output.filtered}
        """