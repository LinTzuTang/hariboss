import pandas as pd
import argparse

def correct_hariboss_data(input_file, error_file, revised_file, output_file):
    # Load the data from CSV files
    df = pd.read_csv(input_file)
    error_df = pd.read_csv(error_file)
    revised_df = pd.read_csv(revised_file, keep_default_na=False)  # Disable the default NaN recognition

    # Filter out rows in the main dataframe that are listed in the error dataframe
    df_filtered = df[~df.set_index(['pdbid', 'rna_chain', 'ligand_chain', 'ligand_id', 'ligand_resnum']).index.isin(error_df.set_index(['pdbid', 'rna_chain', 'ligand_chain', 'ligand_id', 'ligand_resnum']).index)]

    # Append the revised data to the filtered dataframe
    df_updated = df_filtered.append(revised_df, ignore_index=True)
    df_updated.reset_index(drop=True, inplace=True)

    # Save the updated dataframe to a new CSV file
    df_updated.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description='Correct Hariboss data by removing erroneous entries and appending revised entries.')
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('--error_file', type=str, required=True, help='Path to the error CSV file.')
    parser.add_argument('--revised_file', type=str, required=True, help='Path to the revised error CSV file.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output CSV file.')
    args = parser.parse_args()

    correct_hariboss_data(args.input_file, args.error_file, args.revised_file, args.output_file)
    print(f"Data correction completed successfully and saved to '{args.output_file}'.")

if __name__ == "__main__":
    main()
