#!/bin/bash
PRJ=$1
shift 

if [ -z "$PRJ" ]; then
    echo "Usage: ./bcfquery.sh [project_name] REGION=... INCLUDE=... VTYPE=..."
    echo "e.g. ./bcfquery.sh ebola_2014 REGION=\"AF086833.2:1-500\""
    exit 1
fi

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
# Define the actual paths the shell will ask Snakemake to create
CSV="reads/${PRJ}/vcf/query/query_${TIMESTAMP}.csv"
VCF="reads/${PRJ}/vcf/query/query_${TIMESTAMP}.vcf.gz"
TBI="reads/${PRJ}/vcf/query/query_${TIMESTAMP}.vcf.gz.tbi"

snakemake \
    "$CSV" "$VCF" "$TBI" \
    --use-conda \
    --cores 4 \
    --printshellcmds \
    --config PRJNAME="$PRJ" "$@"