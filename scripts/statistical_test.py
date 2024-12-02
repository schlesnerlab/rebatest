import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
import argparse
import numpy as np

def replace_nan_strings(x):
    """
    Replace string representations of NaN with actual NaN values.

    Parameters:
    x: The value to check and replace if necessary.

    Returns:
    The original value or np.nan if it matches a NaN string.
    """
    if x in ['np.nan', 'nan', 'NaN']:
        return np.nan
    return x

def read_input(file):
    """
    Read and preprocess the input file into a DataFrame.

    Parameters:
    file (str): Path to the input file.

    Returns:
    DataFrame: Processed DataFrame with numeric columns.
    """
    df = pd.read_csv(file, sep='\t', header=None, names=['ID', 'Value1', 'Value2'])
    df['Value2'] = df['Value2'].apply(replace_nan_strings)
    df['Value1'] = pd.to_numeric(df['Value1']).astype(float)
    df['Value2'] = pd.to_numeric(df['Value2']).astype(float)
    return df

def perform_fishers_exact_test(df1, df2, empty_df):
    """
    Perform Fisher's exact test for matching rows in two DataFrames.

    Parameters:
    df1 (DataFrame): First input DataFrame.
    df2 (DataFrame): Second input DataFrame.
    empty_df (list): List to store IDs with missing data.

    Returns:
    tuple: IDs, p-values, and direction (less/more) for each test.
    """
    p_values = []
    less_or_more = []
    ids = df1['ID'].tolist()

    for id in ids:
        row1 = df1[df1['ID'] == id]
        row2 = df2[df2['ID'] == id]

        if pd.isna(row2['Value2'].iloc[0]):
            empty_df.append(id)
            p_values.append(1)
            less_or_more.append('-')
        elif pd.isna(row1['Value2'].iloc[0]):
            contingency_table = [[0, 76215], 
                                 [row2['Value1'].values[0], row2['Value2'].values[0]]]
            _, p_value = stats.fisher_exact(contingency_table)
            p_values.append(p_value)
            less_or_more.append('more')
        else:
            contingency_table = [[row1['Value1'].values[0], row1['Value2'].values[0]], 
                                 [row2['Value1'].values[0], row2['Value2'].values[0]]]
            _, p_value = stats.fisher_exact(contingency_table)
            p_values.append(p_value)
            if (row1['Value1'].values[0] / row1['Value2'].values[0]) > (row2['Value1'].values[0] / row2['Value2'].values[0]):
                less_or_more.append('less')
            else:
                less_or_more.append('more')

    return ids, p_values, less_or_more

def adjust_p_values(p_values):
    """
    Adjust p-values using the Benjamini-Hochberg procedure.

    Parameters:
    p_values (list): List of p-values to adjust.

    Returns:
    list: Adjusted p-values.
    """
    _, p_adjusted, _, _ = multitest.multipletests(p_values, method='fdr_bh')
    return p_adjusted

def main(input1, input2, output):
    """
    Main function to perform Fisher's exact tests and save results.

    Parameters:
    input1 (str): Path to the first input file.
    input2 (str): Path to the second input file.
    output (str): Path to the output file.
    """
    # Read input files
    df1 = read_input(input1)
    df2 = read_input(input2)

    empty_df = []

    # Perform Fisher's exact test
    ids_tmp, p_values_tmp, less_or_more_tmp = perform_fishers_exact_test(df1, df2, empty_df)

    # Filter out empty entries
    ids = []
    p_values = []
    less_or_more = []
    for id, p_val, direction in zip(ids_tmp, p_values_tmp, less_or_more_tmp):
        if id not in empty_df:
            ids.append(id)
            p_values.append(p_val)
            less_or_more.append(direction)

    # Adjust p-values if there are any
    p_adjusted = adjust_p_values(p_values) if p_values else []

    # Create the result DataFrame
    result_df = pd.DataFrame({
        'ID': ids,
        'p_value': p_values,
        'adjusted_p_value': p_adjusted,
        'direction': less_or_more
    }).sort_values(by='adjusted_p_value', ascending=True)

    # Save results to output file
    result_df.to_csv(output, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform Fisher's exact test on two files.")
    parser.add_argument('--input1', type=str, required=True, help='Path to the control input file.')
    parser.add_argument('--input2', type=str, required=True, help='Path to the case input file.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file.')

    args = parser.parse_args()
    main(args.input1, args.input2, args.output)
