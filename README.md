# Bioinfomatics-Project
## Main contents
This project contains codes for automating bioinformatics processes using management serves from Snakemake and Nextflow;
Automating processes includes:
1. fastq and fasta loading
2. aligning (bam)
3. bam filtering
4. vcf calling
5. vcf filter
6. variant annotation

Snakemake: 
The snakemake file was designed to utilize its ability to streamline commands for a wild variety of bioinfommatics packages (bwa, minimap2, bedtools vcftools etc) from a centeralize Snakemake file under the main directory, while storing project specific variables and parameters under respective directories as config,yaml, which allows managing multiple projects extremely easy.

After Snakemake run analysis the directory will manage files neatly as example shown below

 workstation_repository/
├── workflow/
│    ├── Snakefile
│    └── envs/
│         ├── aligners.yaml
│         ├── multiqc.yaml
│         └── r_plots.yaml
├── reads/
│     └── project_name/
│          ├── read_1.fastq
│          ├── read_2.fastq
│          ├── config.yaml
│          └── bam/
│               └── alignment.bam
└── refs/
│    └── ref_name/ 
│          ├── ref,fa
│          ├── ref,gb
│          ├── ref.gff
│          ├── bwa/
│               └── various_index_files
│          ├── bowtie2/
│          └── minimap2/
