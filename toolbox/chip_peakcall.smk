import sys
sys.path.append("toolbox/scripts")
from uscs_chromosomes import CHROM_DICTS
import pandas as pd


# =============================================================================
# Peak Caller Variables
# =============================================================================
peakcaller = config.get("PEAKCALL", {}).get("CALLER", "macs3").lower()
g_size = config.get("PEAKCALL", {}).get("G_SIZE", 12100000)
fmt = config.get("PEAKCALL", {}).get("FORMAT", "BAM").upper()
extsize = config.get("PEAKCALL", {}).get("EXTSIZE", 0)
fdr = config.get("PEAKCALL", {}).get("FDR", 0.05)
# pval = config.get("PEAKCALL", {}).get("PVAL", 1e-5)
k_min = config.get("PEAKCALL", {}).get("K_MIN", 6)
k_max = config.get("PEAKCALL", {}).get("K_MAX", 13)
model_val = config.get("PEAKCALL", {}).get("MODEL", None)
model_arg = f"--{model_val}" if model_val else "" 
smooth = config.get("PEAKCALL", {}).get("SMOOTHING", 30)
gem_version = config.get("PEAKCALL", {}).get("GEM_VERSION", "v2.5").lower()
read_dist =  config.get("PEAKCALL", {}).get("READ_DIST", "default.txt")
GEM_REFS_DIR = f"{REFS_DIR}/gem_{gem_version}"

# =============================================================================
# Helper Function
# =============================================================================

runinfo_df = pd.read_csv(f"{READS_DIR}/sra_runinfo.csv")

def get_peakcall_input(wildcards):
    """
    Looks up the treatment sample directly from the wildcard, 
    and finds its corresponding control matching the same group.
    """
    current_sample_row = runinfo_df[runinfo_df["Run"] == wildcards.sample]
    if current_sample_row.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in sra_runinfo.csv")
        
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
        "treatment": f"{READS_DIR}/bam/filtered/{wildcards.sample}.{aligner}.chip.filtered.bam",
        "control": f"{READS_DIR}/bam/filtered/{control_srr}.{aligner}.chip.filtered.bam"
    }

# def get_chip_peakcall_output(wildcards):
#     """
#     Returns the final file path produced by the active peak caller 
#     for a single specific sample wildcard.
#     """
    
#     if peakcaller == "macs3":
#         return {"peaks": f"{READS_DIR}/peaks/macs3/{wildcards.sample}_peaks.narrowPeak"}
#     elif peakcaller == "gem":
#         return {"peaks": f"{READS_DIR}/peaks/gem/{wildcards.sample}_raw.bed"}
#     else:
#         raise ValueError(f"Unsupported CALLER target: {peakcaller}")

# =============================================================================
# Reference Preparation Rules for GEM
# =============================================================================

rule generate_chrom_sizes:
    """generate chrom_sizes file for GEM"""
    input:
        unpack(get_refs)
    output:
        sizes = f"{GEM_REFS_DIR}/chrom.sizes"
    log:
        f"{LOG_DIR}/peakcall/generate_chrom_sizes.log"
    conda:
        "../env/dna_aligner.yaml"  
    shell:
        """
        samtools faidx {input.fasta} > {log} 2>&1
        cut -f1,2 {input.fasta}.fai > {output.sizes} 2>> {log}
        """

rule generate_genome_dir:
    """
    generate genome_dir file for GEM, convert only filename USCS standard, e.g. chrI.fa
    keep clean header as "">I" or ">NC_001133.9" depends on fasta source.
    """
    input:
        unpack(get_refs)
    output:
        chrom_dir = directory(f"{GEM_REFS_DIR}/chromosomes")
    log:
        f"{LOG_DIR}/peakcall/generate_genome_dir.log"
    params:
        out_path = f"{GEM_REFS_DIR}/chromosomes",
        # Pass a map where the key is the RAW header and the value is the UCSC filename
        mapping_str = lambda wildcards: " ".join(
            f"{k},{v}" for k, v in CHROM_DICTS.get(species_low, {}).items()
        )
    shell:
        """
        mkdir -p {params.out_path} > {log} 2>&1
        
        awk -v out="{params.out_path}" -v maps="{params.mapping_str}" '
        BEGIN {{
            n = split(maps, pairs, " ");
            for(i=1; i<=n; i++) {{
                split(pairs[i], kv, ",");
                map[kv[1]] = kv[2];
            }}
        }}
        /^>/ {{
            split($1, a, " ");
            s = substr(a[1], 2); 
            close(f); 
            
            # Find the designated UCSC filename (e.g., "chrI") from the raw header "s"
            # If not found in the dictionary, default back to the original header string
            file_name = (s in map) ? map[s] : s;
            
            # Prepend the directory path and add the physical file extension
            f = out "/" file_name ".fa";
            
            # Write out the raw original header token (e.g., >NC_001133.9) inside the file
            # This ensures perfect matching with your NCBI BAM files
            print ">" s > f;
            next;
        }} 
        {{print > f}}' {input.fasta} 2>> {log}
        """

rule download_and_extract_gem:
    """download and extract gem.jar and Read_Distribution files from MIT GEM page"""
    output:
        jar = f"{GEM_REFS_DIR}/gem.jar",
        exo_dist = f"{GEM_REFS_DIR}/Read_Distribution_ChIP-exo.txt",
        seq_dist = f"{GEM_REFS_DIR}/Read_Distribution_default.txt"
    params:
        target_dir = GEM_REFS_DIR,
        gem_v = gem_version
    conda:"../env/chip_peakcall.yaml" 
    log:
        f"{LOG_DIR}/peakcall/gem/download_and_extract_gem.log"
    shell:
        """
        mkdir -p {params.target_dir} > {log} 2>&1
        
        URL="https://groups.csail.mit.edu/cgs/gem/download/gem.{params.gem_v}.tar.gz"
        wget -O {params.target_dir}/gem_{params.gem_v}.tar.gz "$URL" >> {log} 2>&1
        
        tar -xzvf {params.target_dir}/gem_{params.gem_v}.tar.gz -C {params.target_dir} --strip-components=1 >> {log} 2>&1
        rm {params.target_dir}/gem_{params.gem_v}.tar.gz >> {log} 2>&1
        """

# =============================================================================
# Peak Calling
# =============================================================================

rule peak_calling_macs3:
    """Peak Call with Macs3"""
    input:
        unpack(get_peakcall_input) 
    output:
        peaks = f"{READS_DIR}/peaks/macs3/{{sample}}_peaks.narrowPeak"
    params:
        fmt = fmt,
        g_size = g_size,
        fdr = fdr,
        # pval = pval,
        model = model_arg,
        extsize = extsize
    conda:"../env/chip_peakcall.yaml" 
    log:
        f"{LOG_DIR}/peakcall/macs3/{{sample}}.log"
    shell:
        """
        macs3 callpeak \
            -t {input.treatment} \
            -c {input.control} \
            -f {params.fmt} \
            -g {params.g_size} \
            -q {params.fdr} \
            --extsize {params.extsize} \
            --keep-dup all \
            --outdir {READS_DIR}/peaks/macs3 \
            -n {wildcards.sample} \
            {params.model} > {log} 2>&1
        """

rule peak_calling_gem:
    input:
        unpack(get_peakcall_input),
        chrom_sizes = f"{GEM_REFS_DIR}/chrom.sizes",
        genome_dir = f"{GEM_REFS_DIR}/chromosomes",
        gem_jar = f"{GEM_REFS_DIR}/gem.jar",
        read_dist_path = f"{GEM_REFS_DIR}/Read_Distribution_{read_dist}"
    output:
        events = f"{READS_DIR}/peaks/gem/{{sample}}/{{sample}}_GEM_events.txt",
        pfm    = f"{READS_DIR}/peaks/gem/{{sample}}/{{sample}}_PFM.txt",
        htm    = f"{READS_DIR}/peaks/gem/{{sample}}/{{sample}}_result.htm"
    params:
        fmt = fmt,
        smooth = smooth,
        k_min = k_min,
        k_max = k_max,
        out_prefix = f"{READS_DIR}/peaks/gem/{{sample}}",
        gps_log_dir = f"{LOG_DIR}/peaks/gem"
    conda: "../env/chip_peakcall.yaml" 
    threads: 4
    log:
        f"{LOG_DIR}/peakcall/gem/{{sample}}.log"
    shell:
        """
        mkdir -p {params.gps_log_dir} > {log} 2>&1

        java -Xmx10g -jar {input.gem_jar} \
            --d {input.read_dist_path} \
            --g {input.chrom_sizes} \
            --genome {input.genome_dir} \
            --expt {input.treatment} \
            --ctrl {input.control} \
            --f {params.fmt} \
            --r {params.smooth} \
            --out {params.out_prefix} \
            --min_event_count 20 \
            --k_min {params.k_min} \
            --k_max {params.k_max} \
            --t {threads} >> {log} 2>&1

        # Post-execution cleanup: relocate the stray log file
        if [ -f GPS_Log.txt ]; then
            mv GPS_Log.txt {params.gps_log_dir}/{wildcards.sample}_GPS_Log.txt
        fi
        """

# =============================================================================
# Format convertion for GEM
# =============================================================================
rule gem_to_bed:
    """Converts single-coordinate GPS/GEM events into X basepair BED intervals and sorts them"""
    input:
        events = f"{READS_DIR}/peaks/gem/{{sample}}/{{sample}}_GEM_events.txt"
    output:
        raw_bed = f"{READS_DIR}/peaks/gem/{{sample}}/{{sample}}_GEM_events.bed"
    params:
        padding = 20 // 2
    shell:
        """
        # 1. Skip header (NR>1)
        # 2. split($1, coord, ":") breaks "VII:384892" into coord[1]="VII" and coord[2]=384892
        # 3. Calculate intervals and pull IP (col 2) and Q-score (col 6)
        # 4. Pipe to sort -k1,1 -k2,2n to guarantee chromosome and coordinate order
        
        awk 'NR>1 {{
            split($1, coord, ":");
            print coord[1] "\\t" (coord[2]-{params.padding}) "\\t" (coord[2]+{params.padding}) "\\t" $2 "\\t" $6
        }}' {input.events} | sort -k1,1 -k2,2n > {output.raw_bed}
        """

# =============================================================================
# Peak Filtering for Macs3
# =============================================================================
rule generate_sdust_mask:
    """Runs sdust on the genome fasta to create a low-complexity repeat mask"""
    input:
        unpack(get_refs)
    output:
        sdust_mask = f"{READS_DIR}/peaks/filter_mask_blacklist/sdust_repeats.bed"
    conda: "../env/chip_peakcall.yaml"
    log: f"{LOG_DIR}/peakcall/generate_sdust_mask.log"
    shell:
        """
        sdust {input.fasta} > {output.sdust_mask} 2> {log} 
        """

rule create_raw_blacklist_features:
    """Uses Python to parse the GFF and extract matching raw blacklist features"""
    input:
        unpack(get_refs)
    output:
        mito_raw = temp(f"{READS_DIR}/peaks/filter_mask_blacklist/tmp_mito_raw.bed"),
        rrna_raw = temp(f"{READS_DIR}/peaks/filter_mask_blacklist/tmp_rrna_raw.bed")
    run:
        # Fetch the species translation mapping dictionary
        chrom_map = CHROM_DICTS.get(species_low, {})
        
        with open(input.gff, "r") as f_in, \
             open(output.mito_raw, "w") as f_mito, \
             open(output.rrna_raw, "w") as f_rrna:
             
            for line in f_in:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                
                chrom = parts[0]
                feature = parts[2]
                start = int(parts[3]) - 1
                end = parts[4]
                
                # Check for mitochondrial coordinates using the mapping dictionary
                mapped_chrom = chrom_map.get(chrom, chrom)
                if mapped_chrom == "chrMito":
                    f_mito.write(f"{chrom}\t{start}\t{end}\n")
                
                # Check for ribosomal RNA coordinates
                if feature == "rRNA":
                    f_rrna.write(f"{chrom}\t{start}\t{end}\n")

rule generate_custom_blacklist:
    """Combines, sorts, and merges the raw features into the final blacklist track"""
    input:
        mito_raw = f"{READS_DIR}/peaks/filter_mask_blacklist/tmp_mito_raw.bed",
        rrna_raw = f"{READS_DIR}/peaks/filter_mask_blacklist/tmp_rrna_raw.bed"
    output:
        blacklist = f"{READS_DIR}/peaks/filter_mask_blacklist/blacklist.bed"
    conda: "../env/chip_peakcall.yaml"
    shell:
        """
        # Ensure files exist even if one of the feature categories was empty
        touch {input.mito_raw} {input.rrna_raw}
        
        # Combine, sort coordinates numerically, and merge overlaps
        cat {input.mito_raw} {input.rrna_raw} | sort -k1,1 -k2,2n | bedtools merge > {output.blacklist}
        """


rule chip_peak_filter:
    """Filters either MACS3 or GEM peaks against low-complexity and blacklist tracks"""
    input:
        peaks = f"{READS_DIR}/peaks/macs3/{{sample}}_peaks.narrowPeak",
        blacklist = f"{READS_DIR}/peaks/filter_mask_blacklist/blacklist.bed",
        sdust_mask = f"{READS_DIR}/peaks/filter_mask_blacklist/sdust_repeats.bed"
    output:
        filtered_bed = f"{READS_DIR}/peaks/macs3/{{sample}}_filtered.bed"
    conda: "../env/chip_peakcall.yaml"
    log: f"{LOG_DIR}/peakcall/filter/{{sample}}_bedtools.log"
    shell:
        """
        # Combine masks and clean coordinates
        cat {input.blacklist} {input.sdust_mask} | sort -k1,1 -k2,2n | bedtools merge > tmp_mask_{wildcards.sample}.bed
        
        # Fixed the double stdout redirect bug
        bedtools subtract -a {input.peaks} -b tmp_mask_{wildcards.sample}.bed > {output.filtered_bed} 2> {log}
        
        rm tmp_mask_{wildcards.sample}.bed
        """


# =============================================================================
# Motif analisis for Macs3
# =============================================================================

rule extract_peak_sequences:
    """Extracts fasta sequences from filtered bed coordinates"""
    input:
        unpack(get_refs),
        intervals = f"{READS_DIR}/peaks/macs3/{{sample}}_filtered.bed"
    output:
        result = f"{READS_DIR}/peaks/macs3/{{sample}}_peaks.fasta"
    conda: "../env/chip_peakcall.yaml"
    log: f"{LOG_DIR}/peakcall/motif/{{sample}}_bedtools.log"
    shell:
        """
        bedtools getfasta -fi {input.fasta} -bed {input.intervals} -fo {output.result}
        """

rule motif_discovery_meme_chip:
    """Run MEME-ChIP for comprehensive motif evaluation"""
    input:
        fasta = f"{READS_DIR}/peaks/macs3/{{sample}}_peaks.fasta"
    output:
        out_dir = directory(f"{READS_DIR}/peaks/macs3/motif_analysis/{{sample}}"),
        meme_summary = f"{READS_DIR}/peaks/macs3/motif_analysis/{{sample}}/meme-chip.html"
    conda: "../env/chip_peakcall.yaml"
    params:
        ccut = config.get("CCUT", 100)
    threads: 4
    log:
        f"{LOG_DIR}/peakcall/macs3/motif/{{sample}}_meme-chip.log"
    shell:
        """
        meme-chip \
            -oc {output.out_dir} \
            -meme-p {threads} \
            -dna \
            -ccut {params.ccut} \
            {input.fasta} > {log} 2>&1
        """




