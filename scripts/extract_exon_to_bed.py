import pandas as pd
import argparse

def gtf_to_bed(gtf_file, gene_list, exons_or_genes, bp_extension, tmp_bed_pattern, chromosomes):
    """
    Converts a GTF file into a BED format for either genes or exons, with optional base-pair extension.

    Parameters:
    gtf_file (str): Path to the input GTF file.
    gene_list (list): List of gene names to filter.
    exons_or_genes (str): 'genes' or 'exons' to specify feature type.
    bp_extension (int): Number of base pairs to extend upstream and downstream.
    tmp_bed_pattern (str): Pattern for temporary BED file output names.
    chromosomes (str): Path to output file for unique chromosome names.
    """

    # Read the GTF file into a DataFrame, skipping comment lines
    df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, 
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    if exons_or_genes == 'genes':
        # Filter for gene features only
        genes = df[df['feature'] == 'gene']
        # Extract necessary columns for BED format
        bed_df = genes[['seqname', 'start', 'end', 'attribute']].copy()

    elif exons_or_genes == 'exons':
        # Filter for exon features only
        exons = df[df['feature'] == 'exon']
        # Extract necessary columns for BED format
        bed_df = exons[['seqname', 'start', 'end', 'attribute']].copy()

    else:
        # Raise an error for invalid input
        raise ValueError('exons_or_genes must be "exons" or "genes"')

    # Convert 1-based GTF start positions to 0-based BED start positions and apply base-pair extension
    bed_df['start'] = bed_df['start'].apply(lambda x: max(0, x - 1 - bp_extension))
    # Extend the end position by the specified number of base pairs
    bed_df['end'] = bed_df['end'] - 1 + bp_extension

    # Extract gene name or transcript ID from the 'attribute' column
    bed_df['name'] = bed_df['attribute'].apply(lambda x: extract_gene_name(x))

    # Retain only required BED format columns: chromosome, start, end, and name
    bed_df = bed_df[['seqname', 'start', 'end', 'name']]
    # Rename chromosome column to 'chr' for BED format compliance
    bed_df = bed_df.rename(columns={'seqname': 'chr'})

    # Filter the BED DataFrame for the specified gene list
    bed_df_filtered = bed_df[bed_df['name'].isin(gene_list)]

    # Extract unique chromosome names and save them to the specified file
    unique_chromosomes = bed_df_filtered['chr'].unique()
    # Remove 'chr' prefix from chromosome names, if present
    stripped_values = [value[3:] if value.startswith('chr') else value for value in unique_chromosomes]
    final_unique_values = pd.unique(stripped_values)
    with open(chromosomes, 'w') as file:
        for value in final_unique_values:
            file.write(f"{value}\n")

    # Write filtered BED data to separate files for each unique gene
    write_unique_files(bed_df_filtered, tmp_bed_pattern)

def extract_gene_name(attribute):
    """
    Extracts the gene name from the GTF attribute column.

    Parameters:
    attribute (str): GTF attribute string containing metadata.

    Returns:
    str: Extracted gene name or 'unknown' if not found.
    """
    for attr in attribute.split(';'):
        if 'gene_name' in attr:
            # Extract and clean the gene name value
            return attr.split(' ')[-1].replace('"', '')
    return 'unknown'

def write_unique_files(df, tmp_bed_pattern):
    """
    Writes unique BED files for each gene in the DataFrame.

    Parameters:
    df (DataFrame): Filtered DataFrame in BED format.
    tmp_bed_pattern (str): File name pattern for temporary BED files.
    """
    # Get the unique gene names from the DataFrame
    unique_names = df['name'].unique()

    # Loop over each unique name and save its corresponding BED data to a separate file
    for i, name in enumerate(unique_names, start=1):
        # Filter rows matching the current gene name
        df_filtered = df[df['name'] == name]
        # Define the output file path
        output_file = f"workdir/tmp/{tmp_bed_pattern}_{name}.bed"
        # Write the filtered rows to the output file
        df_filtered.to_csv(output_file, sep='\t', header=False, index=False)

def main():
    """
    Main function to parse arguments and run the GTF to BED conversion.
    """
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description='Convert GTF file exon coordinates to BED format.')
    parser.add_argument('--gtf_file', type=str, help='Input GTF file')
    parser.add_argument('--gene_list', type=str, help='Input gene list')
    parser.add_argument('--exons_or_genes', type=str, help='must be "exons" or "genes"')
    parser.add_argument('--bp_extension', type=int, help='Region extension in base pairs', default=0)
    parser.add_argument('--tmp_bed_pattern', required=False, type=str, help='Output file name pattern', default='tmp')
    parser.add_argument('--chromosomes', type=str, help='Output file for chromosome names')

    # Parse command-line arguments
    args = parser.parse_args()

    # Read the gene list file into a list of gene names
    with open(args.gene_list) as f:
        gene_list = f.read().splitlines()

    # Convert the GTF file to BED format
    gtf_to_bed(args.gtf_file, gene_list, args.exons_or_genes, args.bp_extension, args.tmp_bed_pattern, args.chromosomes)

if __name__ == '__main__':
    main()
