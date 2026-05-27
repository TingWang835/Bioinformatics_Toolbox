rule final_seal:
    """
    Deletes all dummy files created for single-end data processing.
    Run this only when the analysis is complete and you have your VCFs.
    """
        # Ensures this runs only after the final annotation and any reports are done
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
        echo "Starting cleanup of dummy files for project: {PRJNAME}"

        # 1. Remove 0-byte placeholders from getdata and QC
        find {READS_DIR}/ -name "*_2.fastq" -size 0 -delete
        find {READS_DIR}/qc/ -name "*_2_fastqc.html" -size 0 -delete
        find {READS_DIR}/qc/ -name "*_2_fastqc.zip" -size 0 -delete
        find {LOG_DIR}/ -name "*_2_fastq_trimming_report.txt" -size 0 -delete

        # 2. Handle empty Trim Galore gzip dummies (~20 bytes, verify 0 lines)
        find {READS_DIR}/qc_trimmed/ -name "*_2_val_2.fq.gz" -exec sh -c '
            for file; do
                if [ $(zcat "$file" | wc -l) -eq 0 ]; then
                    rm "$file"
                fi
            done
        ' _ {{}} +

        touch {output.marker}
        echo "Cleanup complete. Dummy files removed safely."
        """