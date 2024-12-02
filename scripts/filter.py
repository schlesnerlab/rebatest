import pandas as pd
import argparse
import numpy as np

def filter_and_process_data(input_file, output_file, condition, threshold, region_ID, control_input, annotation, af_threshold):
    """
    Filters data based on annotation thresholds and computes the sum of AC and mean of AN.

    Parameters:
    input_file (str): Path to the input data file.
    output_file (str): Path to the output file.
    condition (str): Whether filtering applies to 'case' or 'gnomAD' data.
    threshold (float): Annotation threshold for filtering.
    region_ID (str): Region identifier for output.
    control_input (str): Optional control input file for AF calculation.
    annotation (str): Annotation column name (e.g., CADD or spliceAI).
    af_threshold (float): AF threshold for filtering.
    """
    # Read input data file
    df = pd.read_csv(input_file, sep='\t', header=0)
    df.drop_duplicates(inplace=True)

    # Check if DataFrame is empty
    if len(df) == 0:
        with open(output_file, 'a') as f:
            f.write(f"{region_ID}\t0\tnp.nan\n")
        return

    if condition == 'case':
        # Process case data
        df = df[['pos_id', annotation, 'AN_cases', 'AC_cases']]
        df.rename(columns={'AN_cases': 'AN', 'AC_cases': 'AC'}, inplace=True)
        backup_mean = df['AN'].mean()

        # Process control data if provided
        control = pd.read_csv(control_input, sep='\t', header=0)
        control.drop_duplicates(inplace=True)
        control = control[['pos_id', annotation, 'AN_nfe', 'AC_nfe']]
        control.rename(columns={'AN_nfe': 'AN', 'AC_nfe': 'AC'}, inplace=True)
        control['AF'] = control['AC'] / control['AN']

        # Merge control AF data with input data
        merged_df = df.merge(control[['pos_id', 'AF']], on='pos_id', how='left')
        merged_df['AF'] = merged_df['AF'].fillna(0)
    else:
        # Process gnomAD data
        df = df[['pos_id', annotation, 'AN_nfe', 'AC_nfe']]
        df.rename(columns={'AN_nfe': 'AN', 'AC_nfe': 'AC'}, inplace=True)
        backup_mean = df['AN'].mean()
        df['AF'] = df['AC'] / df['AN']
        merged_df = df

    # Clean and prepare the annotation column
    merged_df[annotation] = merged_df[annotation].replace('unknown', np.nan)
    merged_df[annotation] = pd.to_numeric(merged_df[annotation], errors='coerce').fillna(0)

    # Filter data based on thresholds
    df_filtered = merged_df[(merged_df[annotation] > threshold) & (merged_df['AC'] > 0) & (merged_df['AF'] < af_threshold)]

    # Calculate results
    AC = df_filtered['AC'].sum() if len(df_filtered) > 0 else 0
    AN = df_filtered['AN'].mean() if len(df_filtered) > 0 else backup_mean

    # Write results to output file
    with open(output_file, 'w') as f:
        f.write(f"{region_ID}\t{AC}\t{AN}\n")

def main():
    """
    Main entry point for filtering and processing CADD or spliceAI data.
    """
    parser = argparse.ArgumentParser(description='Filter for CADD, get sum of AC and mean of AN, and write to file.')
    parser.add_argument('--input', type=str, help='Input file')
    parser.add_argument('--output', type=str, help='Output file')
    parser.add_argument('--condition', type=str, help='Specify whether to filter case or gnomAD data')
    parser.add_argument('--threshold', type=float, help='Threshold for filtering annotation')
    parser.add_argument('--region_ID', type=str, help='Region identifier for the output')
    parser.add_argument('--control_input', type=str, required=False, help='Optional control input file for AF calculation')
    parser.add_argument('--annotation', type=str, required=True, help='Annotation column (e.g., CADD or spliceAI)')
    parser.add_argument('--af_threshold', type=float, required=True, help='AF threshold for filtering')

    args = parser.parse_args()

    # Run the filter and process function
    filter_and_process_data(
        input_file=args.input,
        output_file=args.output,
        condition=args.condition,
        threshold=args.threshold,
        region_ID=args.region_ID,
        control_input=args.control_input,
        annotation=args.annotation,
        af_threshold=args.af_threshold
    )

if __name__ == '__main__':
    main()
