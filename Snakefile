import pandas as pd
import glob
import os
import math


configfile: "scripts/config.yaml"

# Get configfile values
#----------------------------------------------------------
gtf_file = config["gtf_file"]
input_vcf = config["input_vcf"]
fasta = config["fasta"]
sample_family = config["sample_family"]
related_samples = config["relations"]

CADD_threshold = config["CADD_threshold"]
AF_treshold = config["AF_treshold"]
bp_extension_around_gene_pass_to_script = config["bp_extension_around_gene_pass_to_script"]

check_region = config["check_region"]
region = config["region"]

spliceAI_analysis = config["spliceAI_analysis"]
spliceAI_threshold = config["spliceAI_threshold"]

exons_and_genes_pass_to_script = config["exons_and_genes_pass_to_script"]
exons_or_genes_pass_to_script = config["exons_or_genes_pass_to_script"]

all_genes = config["all_genes"]
allgenes_file = config["allgenes_file"]
genes = config["genes"]

cadd_dir_snv = config["vep_cadd_dir_snv"]
cadd_dir_indel = config["vep_cadd_dir_indel"]
spliceai_dir_snv = config["vep_spliceai_dir_snv"]
spliceai_dir_indel = config["vep_spliceai_dir_indel"]
plugin_dir = config["vep_plugin_dir"]
cache_dir = config["vep_cache_dir"]
#----------------------------------------------------------

try:
    os.makedirs("workdir/tmp")
except FileExistsError:
    pass


# Very few steps have to be performed when only regions are to be analyzed
#----------------------------------------------------------
if check_region:

    # Start by 1 to get integer i in the filename from 1 to length
    length = []
    with open(region, 'r') as f:
        for line in f:
            # Remove any leading/trailing whitespace characters
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            # Split the line into columns (assuming tsv / bed)
            columns = line.split('\t')

            # Concatenate the first three columns with underscores
            concatenated_string = '_'.join(columns[:3])
            # Add the concatenated string to the list
            length.append(concatenated_string)


# Steps necessary if exons / genes are to be analyzed
else:


    exons_and_genes = exons_and_genes_pass_to_script if check_region == False else False
    # If exons_and_genes is True, this parameter has to equalize to genes, because exons will be included in the analysis and can be filtered afterwards
    exons_or_genes = "genes" if exons_and_genes == True else exons_or_genes_pass_to_script
    # Extend the region around the gene by this number of base pairs; dont extend if only exons are to be analyzed
    bp_extension_around_gene = bp_extension_around_gene_pass_to_script if exons_or_genes == 'genes' else 0 


    # Get a list of all genes contained in the gtf file
    #----------------------------------------------------------
    if os.path.exists(allgenes_file) == False:

        def extract_gene_name(attribute):
            # Extract gene name or transcript ID from the attribute string
            for attr in attribute.split(';'):
                if 'gene_name' in attr or 'transcript_id' in attr:
                    return attr.split(' ')[-1].replace('"', '')
            return 'unknown'

        # Read GTF file to df
        df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

        # Extract only genes and exclude mitochondrial genes
        genes_df = df[(df['feature'] == 'gene') & (df['seqname'] != 'chrM')]
        # Extract necessary columns for BED format
        bed_df = genes_df[['seqname', 'start', 'end', 'attribute']].copy()
        # Extract gene name
        bed_df['name'] = bed_df['attribute'].apply(lambda x: extract_gene_name(x))
        # Get unique gene names to list
        allgenes_list=bed_df['name'].unique().tolist()

        # Write all gene names to file
        with open(allgenes_file, 'w') as file:
            for gene in allgenes_list:
                file.write(gene + '\n')
    #----------------------------------------------------------


    # Determine on input parameter, which genes should be analyzed
    #----------------------------------------------------------
    # Read allgenes from file to set
    with open(allgenes_file, 'r') as file:
        genes_in_file = set(line.strip() for line in file)

    # If all_genes are to be analyzed, all genes in the gtf file are included
    if all_genes == 1:
        gene_name_list = list(genes_in_file)
    # If only a sublist of genes is to be analyzed, read the gene names from this file and check for their existence in the gtf file
    else:
        gene_name_list = []
        with open(genes, 'r') as file:
            # Read all lines and strip the newline character
            gene_name_list = [line.strip() for line in file]
        # Filter gene_name_list for genes that are in the gtf file
        gene_name_list = [item for item in gene_name_list if item in genes_in_file]
    #----------------------------------------------------------

    length = gene_name_list

    gene_work_file = "workdir/tmp/gene_name_list.txt"
    with open(gene_work_file, "w") as file:
        for gene in gene_name_list:
            file.write(gene + '\n')

#----------------------------------------------------------



# For large number of genes / allgenes, performance is improved, when batching of jobs is done (decreases number of jobs to be submitted)
#----------------------------------------------------------
def group_files(files, batch_size=100):
    return [files[i:i + batch_size] for i in range(0, len(files), batch_size)]

# Construct batches
case_files = [f"workdir/tmp/cases_annotated_{i}.vcf.gz" for i in length]
case_batches = group_files(case_files, batch_size=100)
# Create a dictionary to map batches
cases_batch_map = {i: batch for i, batch in enumerate(case_batches)}
cases_batch_map_length = len(cases_batch_map)

control_files = [f"workdir/tmp/gnomAD_annotated_{i}.vcf.gz" for i in length]
control_batches = group_files(control_files, batch_size=100)
# Create a dictionary to map batches
control_batch_map = {i: batch for i, batch in enumerate(control_batches)}
control_batch_map_length = len(control_batch_map)






## FOR THE MOMENT, chrM IS EXCLUDED! IF THIS SHOULD BE ADAPTED, HARMONIZE_GNOMAD.PY NEEDS TO BE ADAPTED
localrules: download_gnomad


# Define output files
#----------------------------------------------------------
if spliceAI_analysis:
    if exons_and_genes_pass_to_script:
        rule all:
            input:
                os.path.join(f"workdir/results_CADD_{CADD_threshold}.tsv"),
                os.path.join(f"workdir/results_spliceAI_{spliceAI_threshold}.tsv"),
                os.path.join(f"workdir/results_exon_CADD_{CADD_threshold}.tsv"),
                os.path.join(f"workdir/results_exon_spliceAI_{spliceAI_threshold}.tsv")
    else:
        rule all:
                input:
                    os.path.join(f"workdir/results_CADD_{CADD_threshold}.tsv"),
                    os.path.join(f"workdir/results_spliceAI_{spliceAI_threshold}.tsv"),

else:
    if exons_and_genes_pass_to_script:
        rule all:
            input:
                os.path.join(f"workdir/results_CADD_{CADD_threshold}.tsv"),
                os.path.join(f"workdir/results_exon_CADD_{CADD_threshold}.tsv"),
    else:
        rule all:
            input:
                os.path.join(f"workdir/results_CADD_{CADD_threshold}.tsv"),
#----------------------------------------------------------


# Depending on chosen parameters, include additional rules / snakefiles
if spliceAI_analysis:
    include: "scripts/snakefile_spliceAI"
if exons_and_genes:
    include: "scripts/snakefile_exons_and_genes"



# SECTION 0: Prepare the input vcf
#----------------------------------------------------------
rule prepare_cohort_vcf:
    input:
        vcf = input_vcf,
    output:
        os.path.join("workdir/cases.vcf.gz")  
    singularity:
        os.path.join("container/bcftools.sif")
    threads: lambda wildcards, attempt: ( 5 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 40000 * attempt ),
        runtime = lambda wildcards, attempt: ( 200 * attempt ),
        partition = 'epyc-mem',
    params:
        fasta = fasta
    shell:
        """
        bcftools norm --threads {threads} -m-both -f {params.fasta} -Ou {input.vcf} |
        bcftools view --threads {threads} --exclude 'COUNT(GT="0/0" | GT="./.") == N_SAMPLES' -Ou - |
        bcftools annotate --threads {threads} --set-id '%CHROM\:%POS\:%REF\:%ALT' -Ou - | 
        bcftools view --threads {threads} --exclude 'COUNT((GQ >= 20 & DP >= 10 & GT = "het" & (FORMAT/AD[*:1] / (FORMAT/AD[*:1] + FORMAT/AD[*:0])) > 0.2) | (GQ >= 20 & DP >= 10 & GT = "hom")) = 0' -Oz -o {output} 
        
        tabix -p vcf {output}
        """
#----------------------------------------------------------





# SECTION 1: Create input files: xxtract coordinates for regions / genes / exons to be tested; create family ID file
#----------------------------------------------------------
if check_region:

    rule individual_regions:
        input:
            region
        output:
            bed = expand("workdir/tmp/tmp_{i}.bed", i=length),
            chromosomes = "workdir/tmp/chromosomes.txt"
        conda:
            'base'
        threads: lambda wildcards, attempt: ( 1 * attempt ),
        resources: 
            mem_mb = lambda wildcards, attempt: ( 7000 * attempt ),
            runtime = lambda wildcards, attempt: ( 200 * attempt ),
            partition = 'epyc',
        shell:
            "python scripts/individual_regions.py --regions {input} --chromosomes {output.chromosomes}"

else:

    rule extract_exons_or_genes:
        input:
            gtf_file = gtf_file,
        output:
            bed = expand("workdir/tmp/tmp_{i}.bed", i=length),
            chromosomes = "workdir/tmp/chromosomes.txt"
        conda:
            'base'
        threads: lambda wildcards, attempt: ( 1 * attempt ),
        resources: 
            mem_mb = lambda wildcards, attempt: ( 7000 * attempt ),
            runtime = lambda wildcards, attempt: ( 200 * attempt ),
            partition = 'epyc',
        params:
            exons_or_genes = exons_or_genes,
            bp_extension = bp_extension_around_gene,
            gene_work_file = gene_work_file,
        shell:
            "python scripts/extract_exon_to_bed.py --gtf_file {input.gtf_file} --gene_list {params.gene_work_file} --exons_or_genes {params.exons_or_genes} --bp_extension {params.bp_extension} --chromosomes {output.chromosomes}"

if exons_and_genes:
    rule download_gnomad:
        input:
            gene = "workdir/tmp/chromosomes.txt",
            exon = "workdir/tmp/chromosomes_exon.txt"
        output:
            "workdir/tmp/download_complete.txt"
        conda:
            'base'
        threads: lambda wildcards, attempt: ( 1 * attempt ),
        resources: 
            mem_mb = lambda wildcards, attempt: ( 1000 * attempt ),
            runtime = lambda wildcards, attempt: ( 1000 * attempt ),
            partition = 'epyc',
        shell:
            """
            bash scripts/download_gnomad.sh --gene {input.gene} --exon {input.exon} --output {output}
            """
else:
    rule download_gnomad:
        input:
            "workdir/tmp/chromosomes.txt",
        output:
            "workdir/tmp/download_complete.txt"
        conda:
            'base'
        threads: lambda wildcards, attempt: ( 1 * attempt ),
        resources: 
            mem_mb = lambda wildcards, attempt: ( 1000 * attempt ),
            runtime = lambda wildcards, attempt: ( 1000 * attempt ),
            partition = 'epyc',
        shell:
            """
            bash scripts/download_gnomad.sh --gene {input} --output {output}
            """
#----------------------------------------------------------




# SECTION 2: Prepare gnomAD data and annotate with CADD score
#----------------------------------------------------------
rule intersect_gnomAD:
    input:
        bed = os.path.join("workdir/tmp/tmp_{i}.bed"),
        gnomAD = "workdir/tmp/download_complete.txt"
    output:
        os.path.join("workdir/tmp/gnomAD_{i}.vcf.gz")  
    conda:
        'base'
    threads: lambda wildcards, attempt: ( 1 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 2000 * attempt ),
        runtime = lambda wildcards, attempt: ( 400 * attempt ),
        partition = 'epyc',
    params:
        path_to_gnomAD = 'gnomad_variants/gnomad.joint.v4.1.sites.',
        tmp_output = lambda wildcards: os.path.join("workdir/tmp", f"gnomAD_{wildcards.i}.vcf"),
        gnomAD_header = 'reference_files/gnomAD_header.txt'
    shell:
        """
        unique_chromosomes=$(cut -f 1 {input.bed} | sort | uniq)
        cat {params.gnomAD_header} > {params.tmp_output}
        
        for chr in $unique_chromosomes; do
            bed_chr_file="workdir/tmp/tmp_{wildcards.i}_${{chr}}.bed"
            grep -w "$chr" {input.bed} > "$bed_chr_file"


            bedtools intersect -wa -a "{params.path_to_gnomAD}${{chr}}.vcf.bgz" -b ${{bed_chr_file}} >> {params.tmp_output}

            rm "$bed_chr_file"
        done

        gzip {params.tmp_output}
        """


rule annotate_gnomAD:
    input:
        os.path.join("workdir/tmp/gnomAD_{i}.vcf.gz")
    output:
        os.path.join("workdir/tmp/gnomAD_annotated_{i}.vcf.gz") 
    threads: lambda wildcards, attempt: ( 8 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 60000 * attempt ),
        runtime = lambda wildcards, attempt: ( 120 * attempt ),
        partition = 'epyc-mem',
    singularity:
        'container/vep.sif'
    params:
        cadd_dir_snv = cadd_dir_snv,
        cadd_dir_indel = cadd_dir_indel,
        spliceai_dir_snv = spliceai_dir_snv,
        spliceai_dir_indel = spliceai_dir_indel,
        plugin_dir = plugin_dir,
        cache_dir = cache_dir,
        fasta = fasta,
    shell:
        """
        if (( $(zcat {input} | grep -v '^#'  | wc -l) > 0 )); then
            vep \
            -i {input} \
            -o {output} \
            --vcf \
            --fork {threads} \
            --offline \
            --dir_plugins {params.plugin_dir} \
            --plugin CADD,snv={params.cadd_dir_snv},indels={params.cadd_dir_indel} \
            --plugin SpliceAI,snv={params.spliceai_dir_snv},indel={params.spliceai_dir_indel} \
            --no_stats \
            --force_overwrite \
            --fields CADD_PHRED,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL \
            --cache \
            --cache_version 109 \
            --dir_cache {params.cache_dir} \
            --assembly GRCh38  \
            --species homo_sapiens \
            --fasta {params.fasta} \
            --compress_output gzip

            # rm {input}
        else
            mv {input} {output}
        fi
        """


rule harmonize_gnomAD:
    input:
        lambda wildcards: control_batch_map[int(wildcards.i)]
    output:
        "workdir/tmp/ctl_batch_{i}.done"
    conda:
        'base'
    threads: lambda wildcards, attempt: ( 1 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 7000 * attempt ),
        runtime = lambda wildcards, attempt: ( 200 * attempt ),
        partition = 'epyc',
    params:
        exons_and_genes = exons_and_genes,
    shell:
        """
        for in_file in {input}; do
            out_file=$(echo "$in_file" | sed 's/gnomAD_annotated_/gnomAD_df_/' | sed 's/.vcf.gz/.tsv/')
            if (( $(zcat "$in_file" | grep -v '^#'  | wc -l) > 0 )); then
                python scripts/harmonize_gnomAD.py --vcf_input "$in_file" --output "$out_file"
            else
                echo -e 'pos_id\tAN\tAC\tAN_nfe\tAC_nfe\tCADD\tspliceAI' > "$out_file"
            fi
        done

        if [ "{params.exons_and_genes}" = "True" ]; then
            echo ''
        else
            echo ''
            # rm {input}
        fi

        touch {output}
        """
#----------------------------------------------------------



# SECTION 3: Prepare cases_data and combine AN&AC with CADD scores
#----------------------------------------------------------
rule intersect_cases:
    input:
        bed = os.path.join("workdir/tmp/tmp_{i}.bed"),
        vcf = os.path.join("workdir/cases.vcf.gz"),
        gnomAD = "workdir/tmp/download_complete.txt"
    output:
        os.path.join("workdir/tmp/cases_{i}.vcf")  
    singularity:
        'container/bedtools.sif'
    threads: lambda wildcards, attempt: ( 1 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 2000 * attempt ),
        runtime = lambda wildcards, attempt: ( 200 * attempt ),
        partition = 'epyc',
    shell:
        """
        bedtools intersect -wa -header -a {input.vcf} -b {input.bed} > {output}
        """

rule annotate_cases:
    input:
        os.path.join("workdir/tmp/cases_{i}.vcf")
    output:
        os.path.join("workdir/tmp/cases_annotated_{i}.vcf.gz") 
    threads: lambda wildcards, attempt: ( 8 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 60000 * attempt ),
        runtime = lambda wildcards, attempt: ( 60 * attempt ),
        partition = 'epyc-mem',
    singularity:
        'container/vep.sif'
    params:
        cadd_dir_snv = cadd_dir_snv,
        cadd_dir_indel = cadd_dir_indel,
        spliceai_dir_snv = spliceai_dir_snv,
        spliceai_dir_indel = spliceai_dir_indel,
        plugin_dir = plugin_dir,
        cache_dir = cache_dir,
        fasta = fasta,
    shell:
        """
        if (( $(cat {input} | grep -v '^#'  | wc -l) > 0 )); then
            vep \
            -i {input} \
            -o {output} \
            --vcf \
            --fork {threads} \
            --offline \
            --dir_plugins {params.plugin_dir} \
            --plugin CADD,snv={params.cadd_dir_snv},indels={params.cadd_dir_indel} \
            --plugin SpliceAI,snv={params.spliceai_dir_snv},indel={params.spliceai_dir_indel} \
            --no_stats \
            --force_overwrite \
            --fields CADD_PHRED,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL \
            --cache \
            --cache_version 109 \
            --dir_cache {params.cache_dir} \
            --assembly GRCh38  \
            --species homo_sapiens \
            --fasta {params.fasta} \
            --compress_output gzip
        else
            cat {input} | gzip > {output}
        fi
        # rm {input}
        """

rule get_AC_AN_cases:
    input:
        vcf = lambda wildcards: cases_batch_map[int(wildcards.i)],
    output:
        "workdir/tmp/cases_batch_{i}.done"
    conda:
        'base'
    threads: lambda wildcards, attempt: ( 1 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 7000 * attempt ),
        runtime = lambda wildcards, attempt: ( 200 * attempt ),
        partition = 'epyc',
    params:
        exons_and_genes = exons_and_genes,
        family = "related_samples.txt",
    shell:
        """
        if [ ! -e {params.family} ]; then
            touch {params.family}
        fi
        
        for in_file in {input.vcf}; do
            out_file=$(echo "$in_file" | sed 's/cases_annotated_/cases_df_/' | sed 's/.vcf.gz/.tsv/')
            bash scripts/get_AC_AN_cases.sh "$in_file" "$out_file" {params.family}
        done

        if [ "{params.exons_and_genes}" = "True" ]; then
            echo ''
        else
            echo ''
            # rm {input.vcf}
        fi

        touch {output}
        """
#----------------------------------------------------------







# SECTION 4: Filter gnomAD and case data for CADD, get AC and mean(AN) for each region / gene and write to file
#----------------------------------------------------------
rule filter_CADD_cases:
    input:
        expand("workdir/tmp/cases_batch_{i}.done",i=range(cases_batch_map_length)),
        expand("workdir/tmp/ctl_batch_{i}.done",i=range(control_batch_map_length)),
    output:
        os.path.join(f"workdir/case_counts_CADD_{CADD_threshold}.tsv")  
    conda:
        'base'
    threads: lambda wildcards, attempt: ( 5 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 5000 * attempt ),
        runtime = lambda wildcards, attempt: ( 500 * attempt ),
        partition = 'epyc-mem',
    params:
        CADD_threshold = CADD_threshold,
        AF_treshold = AF_treshold,
    shell:
        """
        for case_file in workdir/tmp/cases_df_*.tsv; do
            # Extract the integer i from the filename
            i=$(echo "$case_file" | grep -oP '(?<=cases_df_)[^\.]+(?=\.tsv)')
            gnomad_file="workdir/tmp/gnomAD_df_${{i}}.tsv"
            case_counts_file="workdir/tmp/case_counts_CADD_${{i}}.tsv"
            
            python scripts/filter.py --input $case_file --output $case_counts_file --condition case --threshold {params.CADD_threshold} --af_threshold {params.AF_treshold} --region_ID $i --control_input $gnomad_file --annotation 'CADD' &

            sleep 0.15
        done
        wait

        for file in workdir/tmp/case_counts_CADD_*.tsv; do
            cat $file >> {output}
            rm $file
        done
        """


rule filter_CADD_gnomAD:
    input:
        expand("workdir/tmp/ctl_batch_{i}.done",i=range(control_batch_map_length)),
    output:
        os.path.join(f"workdir/gnomAD_counts_CADD_{CADD_threshold}.tsv")
    conda:
        'base'
    threads: lambda wildcards, attempt: ( 5 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 5000 * attempt ),
        runtime = lambda wildcards, attempt: ( 500 * attempt ),
        partition = 'epyc-mem',
    params:
        CADD_threshold = CADD_threshold,
        AF_treshold = AF_treshold,
    shell:
        """
        for control_file in workdir/tmp/gnomAD_df_*.tsv; do
            # Extract the integer i from the filename
            i=$(echo "$control_file" | grep -oP '(?<=gnomAD_df_)[^\.]+(?=\.tsv)')
            control_counts_file="workdir/tmp/gnomAD_counts_CADD_${{i}}.tsv"
            
            python scripts/filter.py --input $control_file --output $control_counts_file --condition control --threshold {params.CADD_threshold} --af_threshold {params.AF_treshold} --region_ID $i --annotation 'CADD' &

            sleep 0.15
        done
        wait

        for file in workdir/tmp/gnomAD_counts_CADD_*.tsv; do
            cat $file >> {output}
            rm $file
        done
        """
#----------------------------------------------------------




# SECTION 5: Statistical test
#----------------------------------------------------------
rule statistical_test_CADD:
    input:
        gnomAD = os.path.join(f"workdir/gnomAD_counts_CADD_{CADD_threshold}.tsv"),
        cases = os.path.join(f"workdir/case_counts_CADD_{CADD_threshold}.tsv"),
    output:
        os.path.join(f"workdir/results_CADD_{CADD_threshold}.tsv")  
    conda:
        'base'
    threads: lambda wildcards, attempt: ( 1 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 7000 * attempt ),
        runtime = lambda wildcards, attempt: ( 200 * attempt ),
        partition = 'epyc',
    params:
        check_region = check_region,
    shell:
        """
        python scripts/statistical_test.py --input1 {input.gnomAD} --input2 {input.cases} --output {output} --check_region "{params.check_region}"
        """
#----------------------------------------------------------





# Only executed, if rule all is adjusted.
rule remove_tmp_files:
    input:
        lambda wildcards: [f"workdir/results_CADD_{CADD_threshold}.tsv",f"workdir/results_exon_CADD_{CADD_threshold}.tsv",f"workdir/results_spliceAI_{spliceAI_threshold}.tsv",f"workdir/results_exon_spliceAI_{spliceAI_threshold}.tsv"] if exons_and_genes else [f"workdir/results_CADD_{CADD_threshold}.tsv",f"workdir/results_spliceAI_{spliceAI_threshold}.tsv"]
    output:
        "workdir/tmp/done.txt"
    threads: lambda wildcards, attempt: ( 1 * attempt ),
    resources: 
        mem_mb = lambda wildcards, attempt: ( 7000 * attempt ),
        runtime = lambda wildcards, attempt: ( 200 * attempt ),
        partition = 'epyc',
    shell:
        """
        for file in workdir/tmp/tmp*.bed; do
            rm $file
        done
        touch {output}
        """
