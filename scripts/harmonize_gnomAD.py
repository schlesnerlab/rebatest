import pandas as pd
import vcf
import argparse
import os

def extract_info_field(info, key):
    """
    Extract a specific field from the INFO dictionary of a VCF record.

    Parameters:
    info (dict): INFO field of a VCF record.
    key (str): The key to extract.

    Returns:
    The value corresponding to the key, or None if the key is not present.
    """
    return info.get(key, None)

def check_and_get_only_element(variable):
    """
    Check if the variable is a list of length 1 and return the single element; otherwise, return the variable.

    Parameters:
    variable: The variable to check.

    Returns:
    The single element if the variable is a list of length 1, or the variable itself.
    """
    if isinstance(variable, list) and len(variable) == 1:
        return variable[0]
    return variable

def clean_vcf_in_place(vcf_file):
    """
    Replace occurrences of 'Infinity' in the VCF file with '100000', updating the file in place.

    Parameters:
    vcf_file (str): Path to the VCF file to clean.
    """
    temp_file = vcf_file + '.tmp'
    with open(vcf_file, 'r') as infile, open(temp_file, 'w') as outfile:
        for line in infile:
            line = line.replace('Infinity', '100000')
            outfile.write(line)
    os.replace(temp_file, vcf_file)

def process_vcf(vcf_file, output_file):
    """
    Process a VCF file to extract relevant fields and save them to a tab-delimited file.

    Parameters:
    vcf_file (str): Path to the input VCF file.
    output_file (str): Path to the output file.
    """
    try:
        vcf_reader = vcf.Reader(filename=vcf_file, compressed=True)
        records = []

        for record in vcf_reader:
            try:
                pos_ID = f"{record.CHROM}-{record.POS}-{record.REF}-{check_and_get_only_element(record.ALT)}"
                AN = extract_info_field(record.INFO, 'AN_genomes')
                AC = check_and_get_only_element(extract_info_field(record.INFO, 'AC_genomes'))
                AN_nfe = extract_info_field(record.INFO, 'AN_genomes_nfe')
                AC_nfe = check_and_get_only_element(extract_info_field(record.INFO, 'AC_genomes_nfe'))

                # Extract CADD and spliceAI scores from CSQ field
                CSQ = record.INFO.get('CSQ', '')
                CADD = None
                spliceAI = None

                if CSQ:
                    for annotation in CSQ:
                        annotation_fields = annotation.split('|')
                        annotation_fields = ['0' if field == '' else field for field in annotation_fields]
                        CADD = float(annotation_fields[0])  # Assuming CADD score is the first field
                        spliceAI = max(
                            float(annotation_fields[1]),
                            float(annotation_fields[2]),
                            float(annotation_fields[3]),
                            float(annotation_fields[4])
                        )
                        break

                records.append({
                    'pos_id': pos_ID,
                    'AN': AN,
                    'AC': AC,
                    'AN_nfe': AN_nfe,
                    'AC_nfe': AC_nfe,
                    'CADD': CADD,
                    'spliceAI': spliceAI
                })

            except Exception as e:
                print(f"Error processing record at {record.CHROM}-{record.POS}: {e}")
                continue

        # Save results to a DataFrame
        df = pd.DataFrame(records)
        df.to_csv(output_file, index=False, sep='\t')

    except Exception as e:
        print(f"Error reading VCF: {e}. Attempting to clean and reprocess.")
        clean_vcf_in_place(vcf_file)
        process_vcf(vcf_file, output_file)

def main():
    """
    Main function to process the VCF file and output the results.
    """
    parser = argparse.ArgumentParser(description='Harmonize gnomAD VCF data.')
    parser.add_argument('--vcf_input', type=str, required=True, help='Path to the input VCF file.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file.')

    args = parser.parse_args()
    process_vcf(args.vcf_input, args.output)

if __name__ == '__main__':
    main()
