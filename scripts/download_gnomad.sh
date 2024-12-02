#!/bin/bash

# Parse arguments
gene_file=""
exon_file=""
output_file=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gene) gene_file="$2"; shift ;;
        --exon) exon_file="$2"; shift ;;
        --output) output_file="$2"; shift ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [[ -z "$gene_file" || -z "$output_file" ]]; then
    echo "Usage: $0 --gene <gene_file> [--exon <exon_file>] --output <output_file>"
    exit 1
fi

# Combine and get unique values from gene and exon files
combined_values=$(cat "$gene_file")
if [[ -n "$exon_file" ]]; then
    combined_values+=$'\n'$(cat "$exon_file")
fi
unique_values=$(echo "$combined_values" | sort -u)

# Base path and download URL template
base_path='gnomad_variants/gnomad.joint.v4.1.sites.chr{value}.vcf.bgz'
download_url='https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr{value}.vcf.bgz'

# Process each unique value
while IFS= read -r value; do
    # Substitute the value in the file path
    file_path=${base_path//\{value\}/$value}
    url=${download_url//\{value\}/$value}

    # Check if the file exists
    if [ -f "$file_path" ]; then
        echo "File exists: $file_path"
    else
        echo "File does not exist. Downloading: $url"
        # Download the file
        wget "$url" -O "$file_path"
        wget "$url.tbi" -O "$file_path.tbi"
    fi
done <<< "$unique_values"

# Create the output file
touch "$output_file"
