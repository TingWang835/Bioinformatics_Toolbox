rule cleanup:
    """
    Deletes dummy files created for single-end data processing.
    Adapts dependency gating automatically based on the active project type.
    """
    input:
        lambda wildcards: (
            f"{READS_DIR}/expression/plots/enrichment_go_dotplot.png"
            if "toolbox/rna_diff_exp.smk" in workflow.included_stack
            else f"{READS_DIR}/vcf/all_samples.{dna_aligner}.snpeff_stats.html"
            if "toolbox/dna_vcf.smk" in workflow.included_stack 
            else []
        )
    output:
        marker = f"{LOG_DIR}/cleanup_complete.done"
    shell:
        """
        echo "Starting clean-up of dummy files for project: {PRJNAME}"

        # 1. Force remove absolute 0-byte structural placeholders if they exist
        find {READS_DIR}/ -name "*_2.fastq" -type f -size 0 -delete
        find {READS_DIR}/qc/ -name "*_2_fastqc.html" -type f -size 0 -delete
        find {READS_DIR}/qc/ -name "*_2_fastqc.zip" -type f -size 0 -delete
        find {LOG_DIR}/ -name "*_2_fastq_trimming_report.txt" -type f -size 0 -delete

        # 2. Handle Trim Galore gzip dummies (escaped for Snakemake via double curly braces)
        find {READS_DIR}/qc_trimmed/ -name "*_2_val_2.fq.gz" -type f -exec sh -c '
            for file; do
                if [ -f "$file" ] && [ $(zcat "$file" | head -n 4 | wc -l) -eq 0 ]; then
                    echo "Deleting empty compressed dummy: $file"
                    rm "$file"
                fi
            done
        ' _ {{}} +

        touch {output.marker}
        echo "Cleanup complete. Marker file generated: {output.marker}"
        """