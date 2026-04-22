rule final_seal:
    """
    Deletes all dummy files created for single-end data processing.
    Run this only when the analysis is complete and you have your VCFs.
    """
    input:
        # Ensures this runs only after the final annotation and any reports are done
        expand(f"reads/{{prj}}/vcf/all_samples.{{aligner}}.snpeff_stats.html", 
               prj=PRJNAME, aligner=config.get('ALIGNER', 'bwa'))
    output:
        marker = f"reads/{PRJNAME}/logs/cleanup_complete.done"
    shell:
        """
        echo "Starting cleanup of dummy files for project: {PRJNAME}"

        # 1. Remove 0-byte placeholders from getdata and QC
        # This covers R2 FASTQs, FastQC reports, and Trim Galore logs
        find reads/{PRJNAME}/ -name "*_2.fastq" -size 0 -delete
        find reads/{PRJNAME}/qc/ -name "*_2_fastqc.html" -size 0 -delete
        find reads/{PRJNAME}/qc/ -name "*_2_fastqc.zip" -size 0 -delete
        find reads/{PRJNAME}/logs/ -name "*_2_fastq_trimming_report.txt" -size 0 -delete

        # 2. Specifically handle the empty Trim Galore gzip dummies(which are ~20 bytes, not 0)
        # Verify they contain 0 lines before deleting
        find reads/{PRJNAME}/qc_trimmed/ -name "*_2_val_2.fq.gz" -exec sh -c 'if [ $(zcat {{}} | wc -l) -eq 0 ]; then rm {{}}; fi' \;

        touch {output.marker}
        echo "Cleanup complete. Dummy files removed."
        """