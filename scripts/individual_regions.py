import argparse

def create_temp_files(bed_file, necessary_chromosomes):
    """
    Splits a BED file into individual temporary files, one for each line, and extracts necessary chromosome information.

    Parameters:
    bed_file (str): Path to the input BED file.
    necessary_chromosomes (list): List to store extracted chromosome names.
    """
    with open(bed_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines, start=1):
        temp_file_name = f'workdir/tmp/tmp_{i}.bed'
        with open(temp_file_name, 'w') as temp_file:
            temp_file.write(line)
        # Extract chromosome name from the first column
        necessary_chromosomes.append(line.split('\t')[0].split('chr')[1])

def main():
    """
    Main function to convert a regions BED file into individual temporary files and output the necessary chromosome names.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert a region file into individual temporary files.')
    parser.add_argument('--regions', type=str, required=True, help='Input regions BED file.')
    parser.add_argument('--chromosomes', type=str, required=True, help='Output file for chromosome names.')

    # Parse arguments
    args = parser.parse_args()

    necessary_chromosomes = []

    # Create temporary files and collect chromosome names
    create_temp_files(args.regions, necessary_chromosomes)

    # Write unique chromosome names to the output file
    necessary_chromosomes = set(necessary_chromosomes)
    with open(args.chromosomes, 'w') as f:
        for chromosome in necessary_chromosomes:
            f.write(f"{chromosome}\n")

if __name__ == '__main__':
    main()
