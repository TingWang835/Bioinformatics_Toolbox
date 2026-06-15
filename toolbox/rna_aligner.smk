import os
wildcard_constraints:
    aligner = "star|hisat2|kallisto|salmon",
    filename = "abundance\.tsv|quant\.sf"
localrules:  # no local rules here due to heavy ram requirement for RNA alignment

# =============================================================================
# Helper Functions
# =============================================================================

def get_rna_index_prefix(wildcards):
    """Retrieves the reference index path based on the aligner type."""
    if wildcards.aligner == "star":
        return f"{REFS_DIR}/star_index"
    elif wildcards.aligner == "hisat2":
        return f"{REFS_DIR}/hisat2/index"
    elif wildcards.aligner == "kallisto":
        return f"{REFS_DIR}/kallisto.idx"
    elif wildcards.aligner == "salmon":
        return f"{REFS_DIR}/salmon_index"
    return f"{REFS_DIR}/{wildcards.aligner}"



# =============================================================================
# Index Rules
# =============================================================================

rule star_index:
    input: 
        unpack(get_refs)
    output: directory(f"{REFS_DIR}/star_index")
    log: f"{LOG_DIR}/aligner/star_index.log"
    conda: "../env/rna_aligner.yaml"
    params: sa = config.get("STAR_SA", "10") 
    threads: 4
    shell:
        """
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --runThreadN {threads} \
             --genomeSAindexNbases {params.sa} > {log} 2>&1
        """

rule hisat2_index:
    input: unpack(get_refs)
    output: multiext(f"{REFS_DIR}/hisat2/index", ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".rev.1.ht2", ".rev.2.ht2")
    log: f"{LOG_DIR}/aligner/hisat2_index.log"
    conda: "../env/rna_aligner.yaml"
    params: prefix = f"{REFS_DIR}/hisat2/index"
    shell: "hisat2-build {input.fasta} {params.prefix} > {log} 2>&1"

rule genome_gtf_to_transcriptome: 
    """Extracts cDNA sequences from the genome using GTF coordinates, used by pseudo aligners."""
    input: unpack(get_refs)
    output:
        fasta = f"{REFS_DIR}/transcriptome.fa"
    conda: "../env/rna_aligner.yaml"
    log: f"{LOG_DIR}/getdata/transcriptome_extract.log"
    shell:
        "gffread {input.gtf} -g {input.fasta} -w {output.fasta} > {log} 2>&1"

rule kallisto_index:
    input: fasta = f"{REFS_DIR}/transcriptome.fa"
    output: idx = f"{REFS_DIR}/kallisto.idx"
    conda: "../env/rna_aligner.yaml"
    log: f"{LOG_DIR}/aligner/kallisto_index.log"
    shell: "kallisto index -i {output.idx} {input.fasta} > {log} 2>&1"

rule salmon_index:
    input: fasta = f"{REFS_DIR}/transcriptome.fa"
    output: idx_dir = directory(f"{REFS_DIR}/salmon_index")
    conda: "../env/rna_aligner.yaml"
    log: f"{LOG_DIR}/aligner/salmon_index.log"
    shell: "salmon index -t {input.fasta} -i {output.idx_dir} > {log} 2>&1"

# =============================================================================
# Alignment Rules
# =============================================================================

rule rna_sequence_align:
    """Standard spliced alignment for STAR and HISAT2 (BAM output)."""
    input:
        unpack(get_align_input),
        idx_check = get_rna_index_prefix
    output:
        bam = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.bam",
        bai = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.bam.bai"
    log: f"{LOG_DIR}/aligner/rna_align/{{sample}}.{{aligner}}.log"
    conda: "../env/rna_aligner.yaml"
    threads: 4
    params:
        cmd = lambda w: (
            f"STAR --genomeDir {get_rna_index_prefix(w)} "
            f"--outTmpDir {READS_DIR}/bam/STARtmp_{w.sample} "
            f"--outFileNamePrefix {LOG_DIR}/aligner/rna_align/{w.sample}_STAR_ "
            f"--outStd SAM "
            f"--readFilesCommand gzip -d -c" if w.aligner == "star" else
            f"hisat2 -x {get_rna_index_prefix(w)} --dta" if w.aligner == "hisat2" else ""
        ),
        pe = lambda w, input: (
            f"--readFilesIn {input.r1} {input.r2}" if w.aligner == "star" else
            f"-1 {input.r1} -2 {input.r2}"
        ),
        se = lambda w, input: (
            f"--readFilesIn {input.r1}" if w.aligner == "star" else
            f"-U {input.r1}"
        )
    shell:
        """
        if [ -s "{input.r2}" ]; then
            # Paired-End Pipeline
            {params.cmd} {params.pe} --runThreadN {threads} | samtools sort -@ {threads} -o {output.bam}
        else
            # Single-End Pipeline
            {params.cmd} {params.se} --runThreadN {threads} | samtools sort -@ {threads} -o {output.bam}
        fi
        samtools index {output.bam}
        """


rule rna_pseudo_align:
    """Fast pseudo-alignment for Kallisto and Salmon (Abundance/Counts output)."""
    input:
        unpack(get_align_input),
        idx_check = get_rna_index_prefix
    output:
        quant_file = f"{READS_DIR}/counts/pseudo_individual/{{sample}}.{{aligner}}/{{filename}}"
    log: f"{LOG_DIR}/aligner/pseudo_align/{{sample}}.{{aligner}}.{{filename}}.log"
    conda: "../env/rna_aligner.yaml"
    params:
        out_dir  = f"{READS_DIR}/counts/pseudo_individual/{{sample}}.{{aligner}}",
        frag_len = config.get("KALLISTO_FRAG_LEN", 200),
        frag_sd  = config.get("KALLISTO_FRAG_SD", 20),
        cmd = lambda w: (
            "kallisto quant" if w.aligner == "kallisto" else
            "salmon quant" if w.aligner == "salmon" else ""
        ),
        pe = lambda w, input: (
            f"-i {get_rna_index_prefix(w)} -o {READS_DIR}/counts/pseudo_individual/{w.sample}.kallisto {input.r1} {input.r2}" if w.aligner == "kallisto" else
            f"-i {get_rna_index_prefix(w)} -l A -1 {input.r1} -2 {input.r2} -o {READS_DIR}/counts/pseudo_individual/{w.sample}.salmon"
        ),
        # Corrected: Removed 'params' argument and accessed config directly
        se = lambda w, input: (
            f"--single -l {config.get('KALLISTO_FRAG_LEN', 200)} -s {config.get('KALLISTO_FRAG_SD', 20)} -i {get_rna_index_prefix(w)} -o {READS_DIR}/counts/pseudo_individual/{w.sample}.kallisto {input.r1}" if w.aligner == "kallisto" else
            f"-i {get_rna_index_prefix(w)} -l A -r {input.r1} -o {READS_DIR}/counts/pseudo_individual/{w.sample}.salmon"
        )
    shell:
        """
        if [ -s "{input.r2}" ]; then
            # Paired-End Execution
            {params.cmd} {params.pe} > {log} 2>&1
        else
            # Single-End Execution
            {params.cmd} {params.se} > {log} 2>&1
        fi
        """
# =============================================================================
# Bam file (Sequence Alignment) after math
# =============================================================================

rule bam_to_bigwig:
    """
    Converts spliced RNA alignments (STAR/HISAT2) into normalized TPM BigWig coverage tracks 
    using deepTools bamCoverage for optimized loading and row-scaling inside IGV.
    """
    input:
        bam = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.bam",
        bai = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.bam.bai"
    output:
        bw = f"{READS_DIR}/bam/bigwig/{{sample}}.{{aligner}}.bw"
    log:
        f"{LOG_DIR}/aligner/bigwig/{{sample}}.{{aligner}}.log"
    conda:
        "../env/rna_aligner.yaml"
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

rule rna_featurecount:
    """
    Quantifies gene-level expression from individual BAM files.
    Outputs raw featureCounts text reports containing genomic annotations.
    """
    input:
        unpack(get_refs),
        bam = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.bam",
        bai = f"{READS_DIR}/bam/{{sample}}.{{aligner}}.bam.bai"
    output:
        txt = f"{READS_DIR}/counts/sequence_individual/{{sample}}.{{aligner}}.txt"
    log:
        f"{LOG_DIR}/aligner/featurecounts/{{sample}}.{{aligner}}.log"
    conda:
        "../env/rna_aligner.yaml"
    threads: 2
    params:
        paired = lambda w: "-p" if os.path.exists(f"{READS_DIR}/{w.sample}_2.fastq") and os.path.getsize(f"{READS_DIR}/{w.sample}_2.fastq") > 0 else "",
        feature_type = config.get("FEATURE_TYPE", "exon"),
        id_attr = config.get("ID_ATTR", "gene_id")
    shell:
        """
        featureCounts {params.paired} \
                      -T {threads} \
                      -a {input.gtf} \
                      -t {params.feature_type} \
                      -g {params.id_attr} \
                      -o {output.txt} \
                      {input.bam} > {log} 2>&1
        """

rule merge_featurecounts:
    """
    Collects individual featureCounts files across all active samples 
    and merges them column-by-column into a single consolidated expression CSV matrix.
    """
    input:
        txts = lambda w: expand(f"{READS_DIR}/counts/sequence_individual/{{sample}}.{w.aligner}.txt", sample=get_runinfo(w))
    output:
        csv = f"{READS_DIR}/counts/sequence/all_samples.{{aligner}}_counts.csv"
    log:
        f"{LOG_DIR}/aligner/featurecounts/merge_all_samples.{{aligner}}.log"
    run:
        import pandas as pd

        # Read the first file to initialize the tracking dataframe with GeneID
        first_file = input.txts[0]
        # featureCounts has a 2-line header; data begins on line 2 (0-indexed)
        master_df = pd.read_csv(first_file, sep="\t", skiprows=1)[["Geneid"]]

        for filepath in input.txts:
            # Extract sample name dynamically out of file path structure
            # e.g., reads/project/counts/individual/SRR11289276.star.txt -> SRR11289276
            sample_name = filepath.split("/")[-1].split(".")[0]
            
            # Read file, keeping only Geneid and the final Count column
            df = pd.read_csv(filepath, sep="\t", skiprows=1)
            count_col = df.columns[-1]
            df = df[["Geneid", count_col]].rename(columns={count_col: sample_name})
            
            # Merge into master matrix
            master_df = pd.merge(master_df, df, on="Geneid", how="outer")

        # Save clean integrated output
        master_df.to_csv(output.csv, index=False)

rule merge_pseudo:
    """
    Collects individual Kallisto abundance files across all active samples 
    and merges their estimated counts column-by-column into a single master matrix.
    """
    input:
        tsvs = lambda w: expand(f"{READS_DIR}/counts/pseudo_individual/{{sample}}.{w.aligner}/abundance.tsv", sample=get_runinfo(w))
    output:
        csv = f"{READS_DIR}/counts/pseudo/all_samples.{{aligner}}_pseudo_counts.csv"
    log:
        f"{LOG_DIR}/aligner/featurecounts/merge_all_samples.{{aligner}}_pseudo.log"
    run:
        import pandas as pd

        # Read the first file to initialize the tracking matrix with transcript 'target_id'
        first_file = input.tsvs[0]
        master_df = pd.read_csv(first_file, sep="\t")[["target_id"]]

        for filepath in input.tsvs:
            # Extract sample name out of file path structure
            # e.g., reads/project/counts/pseudo_individual/SRR11289276.kallisto/abundance.tsv -> SRR11289276
            sample_name = filepath.split("/")[-2].split(".")[0]
            
            # Read file, keeping only target_id and est_counts
            df = pd.read_csv(filepath, sep="\t")[["target_id", "est_counts"]]
            # Rename 'est_counts' column to the current sample name
            df = df.rename(columns={"est_counts": sample_name})
            
            # Merge into master matrix
            master_df = pd.merge(master_df, df, on="target_id", how="outer")

        # Save clean integrated output matrix
        master_df.to_csv(output.csv, index=False)