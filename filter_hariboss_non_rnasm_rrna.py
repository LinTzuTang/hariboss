import pandas as pd
import argparse
import os
import shutil

def filter_complexes(hariboss_csv, non_rnasm_csv, rrna_csv, output_csv, input_dir, output_dir, filter_option, copy_files):
    # Read the data
    df = pd.read_csv(hariboss_csv, keep_default_na=False)
    df_non_rnasm = pd.read_csv(non_rnasm_csv)

    if filter_option == 'non_rnasm':
        # Remove non-RNA SM complexes
        df_filtered = df[~df["pdbid"].isin(df_non_rnasm['non_RNA_SM_complexes'])]
        removed_non_rnasm = df[df['pdbid'].isin(df_non_rnasm['non_RNA_SM_complexes'])]
        removed_non_rnasm.to_csv('removed_non_rna_sm_complexes.csv', index=False)
    elif filter_option == 'non_rnasm_rrna':
        # Remove non-RNA SM complexes
        df_filtered = df[~df["pdbid"].isin(df_non_rnasm['non_RNA_SM_complexes'])]
        removed_non_rnasm = df[df['pdbid'].isin(df_non_rnasm['non_RNA_SM_complexes'])]
        removed_non_rnasm.to_csv('removed_non_rna_sm_complexes.csv', index=False)

        # Read ribosomal complexes data
        df_ribosome = pd.read_csv(rrna_csv)
        
        # Remove ribosomal complexes
        removed_ribosome = df_filtered[df_filtered['pdbid'].isin(df_ribosome['id'])]
        removed_ribosome.to_csv('removed_ribosomal_complexes.csv', index=False)
        df_filtered = df_filtered[~df_filtered['pdbid'].isin(df_ribosome['id'])]
    else:
        print("Invalid filter option. Choose 'non_rnasm' or 'non_rnasm_rrna'.")
        return

    # Save the filtered data to a new CSV file
    df_filtered.to_csv(output_csv, index=False)

    if copy_files:
        # Create the output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Copy non-removed CIF files to the new directory
        non_removed_pdbids = set(df_filtered['pdbid'])
        for filename in os.listdir(input_dir):
            pdbid = filename.split('.')[0]
            if pdbid in non_removed_pdbids:
                src_path = os.path.join(input_dir, filename)
                dest_path = os.path.join(output_dir, filename)
                shutil.copy(src_path, dest_path)

def main():
    parser = argparse.ArgumentParser(description="Filter out non-RNA SM and/or ribosomal complexes from Hariboss dataset and optionally copy filtered CIF files to a new directory")
    parser.add_argument('hariboss_csv', type=str, help="Path to the Hariboss CSV file")
    parser.add_argument('non_rnasm_csv', type=str, help="Path to the non RNA SM complexes CSV file")
    parser.add_argument('--rrna_csv', type=str, help="Path to the ribosomal complexes CSV file (required if filter_option is 'non_rnasm_rrna')")
    parser.add_argument('output_csv', type=str, help="Path to the output CSV file")
    parser.add_argument('--input_dir', type=str, help="Directory containing the original CIF files (required if --copy_files is specified)")
    parser.add_argument('--output_dir', type=str, help="Directory to save the filtered CIF files (required if --copy_files is specified)")
    parser.add_argument('filter_option', type=str, help="Filter option: 'non_rnasm' or 'non_rnasm_rrna'")
    parser.add_argument('--copy_files', action='store_true', help="Option to copy filtered CIF files to the output directory")
    
    args = parser.parse_args()

    if args.filter_option == 'non_rnasm_rrna' and args.rrna_csv is None:
        parser.error("rrna_csv is required when filter_option is 'non_rnasm_rrna'")
    
    if args.copy_files and (args.input_dir is None or args.output_dir is None):
        parser.error("input_dir and output_dir are required when --copy_files is specified")
    
    filter_complexes(args.hariboss_csv, args.non_rnasm_csv, args.rrna_csv, args.output_csv, args.input_dir, args.output_dir, args.filter_option, args.copy_files)

if __name__ == "__main__":
    main()

# For filtering only non-RNA SM complexes without copying files:
# python filter_hariboss_non_rnasm_rrna.py hariboss_clean_20240621_error_revised.csv non_RNA_SM_complexes.csv hariboss_clean_20240621_error_revised_non_rnasm_removed.csv non_rnasm

# For filtering both non-RNA SM and ribosomal complexes without copying files:
# python filter_hariboss_non_rnasm_rrna.py hariboss_clean_20240621_error_revised.csv non_RNA_SM_complexes.csv --rrna_csv Ribosomal_Complexes.csv hariboss_clean_20240621_error_revised_non_rnasm_removed_rrna_removed.csv non_rnasm_rrna

# For filtering and copying the filtered CIF files to the new directory:
# python filter_hariboss_non_rnasm_rrna.py hariboss_clean_20240621_error_revised.csv non_RNA_SM_complexes.csv hariboss_clean_20240621_error_revised_non_rnasm_removed.csv --input_dir complex_mmcif_database --output_dir complex_mmcif_database_non_rnasm_removed non_rnasm --copy_files
# python filter_hariboss_non_rnasm_rrna.py hariboss_clean_20240621_error_revised.csv non_RNA_SM_complexes.csv --rrna_csv Ribosomal_Complexes.csv hariboss_clean_20240621_error_revised_non_rnasm_removed_rrna_removed.csv --input_dir complex_mmcif_database --output_dir complex_mmcif_database_non_rnasm_removed_rrna_removed --copy_files
