#!/bin/bash

# This script processes a VCF file to extract sample-level information, calculates allele counts (AC),
# allele numbers (AN), and annotation scores like CADD and spliceAI, and outputs the processed data.

# Input arguments:
#   $1: Input VCF file
#   $2: Output file for processed data
#   $3: Family file linking samples to family IDs

# Extract sample names from the VCF file header
sample_names=$(bcftools query -l "$1")

# Write the header to the output file
echo -e 'pos_id\tAN_cases\tAC_cases\tCADD\tspliceAI' > "$2"

# Extract necessary fields from the VCF using bcftools and process with awk
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\t%FORMAT[\t%GT]\n' "$1" | \
awk -F'\t' -v family_file="$3" -v sample_names="$sample_names" '
BEGIN {
    # Load family information into an array
    while ((getline < family_file) > 0) {
        split($0, fields, "\t");
        sample_to_family[fields[1]] = fields[2];
    }

    # Split sample names into an array
    split(sample_names, samples, "\n");
}
{
    total_alt_alleles = 0;
    total_alleles = 0;
    chr_pos_ref_alt = $1 "-" $2 "-" $3 "-" $4;

    # Initialize family tracking
    delete family_genotypes;

    # Extract the CADD score and spliceAI values from the INFO field
    cadd_score = "unknown";
    splice_AI = "unknown";
    split($5, info_fields, ";");
    for (i in info_fields) {
        if (info_fields[i] ~ /^CSQ=/) {
            split(info_fields[i], csq_values, "=");
            split(csq_values[2], csq_entries, ",");
            split(csq_entries[1], scores, "|");
            cadd_score = scores[1];

            # Extract maximum spliceAI value
            spliceAI_1 = scores[2]; spliceAI_2 = scores[3];
            spliceAI_3 = scores[4]; spliceAI_4 = scores[5];
            splice_AI = spliceAI_1;
            if (spliceAI_2 > splice_AI) splice_AI = spliceAI_2;
            if (spliceAI_3 > splice_AI) splice_AI = spliceAI_3;
            if (spliceAI_4 > splice_AI) splice_AI = spliceAI_4;
            break;
        }
    }

    # Parse FORMAT field to determine positions of genotype-related fields
    split($6, format_fields, ":");
    for (i in format_fields) {
        if (format_fields[i] == "GT") gt_pos = i;
        else if (format_fields[i] == "DP") dp_pos = i;
        else if (format_fields[i] == "AD") ad_pos = i;
        else if (format_fields[i] == "GQ") gq_pos = i;
    }

    # Process genotype information for each sample
    for (i = 7; i <= NF; i++) {
        sample_name = samples[i - 6];
        split($i, genotype_fields, ":");
        split(genotype_fields[gt_pos], alleles, "/");
        split(genotype_fields[ad_pos], ad_values, ",");
        dp = genotype_fields[dp_pos];
        gq = genotype_fields[gq_pos];

        # Skip low-quality genotypes
        if (dp < 10 || gq < 20) continue;
        if (alleles[1] != alleles[2] && (ad_values[1] + ad_values[2]) == 0) continue;
        if (alleles[1] != alleles[2] && (ad_values[1] / (ad_values[1] + ad_values[2])) < 0.2) continue;

        gt_count = 0;
        for (j in alleles) {
            if (alleles[j] != ".") {
                if (alleles[j] != "0") {
                    gt_count++;
                }
            }
        }

        # Update allele counts based on family membership
        if (sample_name in sample_to_family) {
            family_id = sample_to_family[sample_name];
            if (family_id != "") {
                if (family_id in family_genotypes) {
                    family_genotypes[family_id][length(family_genotypes[family_id]) + 1] = gt_count;
                } else {
                    family_genotypes[family_id][1] = gt_count;
                }
            }
        } else {
            total_alleles += 2;
            total_alt_alleles += gt_count;
        }
    }

    # Adjust counts for family genotypes
    for (family_id in family_genotypes) {
        max_gt = 0;
        for (k in family_genotypes[family_id]) {
            if (family_genotypes[family_id][k] > max_gt) {
                max_gt = family_genotypes[family_id][k];
            }
        }
        total_alt_alleles += max_gt;
        total_alleles += 2;
    }

    # Output the results
    print chr_pos_ref_alt, total_alleles, total_alt_alleles, cadd_score, splice_AI;
}' OFS='\t' >> "$2"
