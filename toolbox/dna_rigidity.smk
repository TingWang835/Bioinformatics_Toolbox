localrules: detect_dna_rigidity, generate_igv_track

# this is a pet project inspired by Gemini 3.0 when we chatted about TurboQuant 
# where we found similarity between the this latest spherical compression logic 
# and the mechanism of DNA compression by histone.
# This code chunk is written to celebrate our little discovery on the unification.

rule detect_dna_rigidity:
    input:
        fasta = f"refs/{config['REFNAME']}/{config['ACC']}.fa"
    output:
        tsv = f"reads/{PRJNAME}/dna_rigidity/{config['ACC']}_rigidity_profile.tsv"
    conda:
        "../env/rigidity.yaml"
    params:
        script = "toolbox/scripts/dna_rigidity.py"
    shell:
        "python {params.script} {input.fasta} {output.tsv}"


rule generate_igv_track:
    input:
        tsv = f"reads/{PRJNAME}/dna_rigidity/{config['ACC']}_rigidity_profile.tsv"
    output:
        bg = f"reads/{PRJNAME}/dna_rigidity/{config['ACC']}_rigidity_profile.bedGraph"
    shell:
        "awk 'NR>1 {{print $1 \"\\t\" $2 \"\\t\" $2+1 \"\\t\" $3}}' {input.tsv} > {output.bg}"