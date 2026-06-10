# this is a pet project inspired by Gemini 3.5 when we chatted about TurboQuant 
# where we found similarity between the this latest spherical compression logic 
# and the mechanism of DNA compression by histone.
# This code chunk is written to celebrate our little discovery on the unification.

localrules: detect_dna_rigidity, generate_igv_track

# =============================================================================
# Rules
# =============================================================================

rule detect_dna_rigidity:
    input:
        unpack(get_refs)
    output:
        tsv = f"{READS_DIR}/dna_rigidity/{acc}_rigidity_profile.tsv"
    conda:
        "../env/dna_rigidity.yaml"
    log: 
        f"{LOG_DIR}/dna_rigidity/{acc}_rigidity_profile.tsv.log"
    params:
        script = "toolbox/scripts/dna_rigidity.py"
    shell:
        "python {params.script} {input.fasta} {output.tsv} > {log} 2>&1"


rule generate_igv_track:
    input:
        tsv = f"{READS_DIR}/dna_rigidity/{acc}_rigidity_profile.tsv"
    output:
        bg = f"{READS_DIR}/dna_rigidity/{acc}_rigidity_profile.bedGraph"
    log: 
        f"{LOG_DIR}/dna_rigidity/{acc}_rigidity_profile.bg.log"
    shell:
        "awk 'NR>1 {{print $1, $2, $2+1, $3}}' {input.tsv} > {output.bg} 2> {log}"