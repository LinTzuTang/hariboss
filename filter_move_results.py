import pandas as pd
import os
import shutil

def filter_non_rrna_files(input_csv, input_dir, sub_dirs, output_dir):
    # Read the non-rRNA PDB IDs from the CSV file
    df_rrna_removed = pd.read_csv(input_csv, keep_default_na=False)
    non_rrna_pdbids = set(df_rrna_removed['pdbid'])

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Filter and copy non-rRNA CIF files from the input directories to the new directory
    for sub_dir in sub_dirs:
        sub_input_dir = os.path.join(input_dir, sub_dir)
        sub_output_dir = os.path.join(output_dir, sub_dir)

        if not os.path.exists(sub_output_dir):
            os.makedirs(sub_output_dir)

        for filename in os.listdir(sub_input_dir):
            pdbid = filename.split('_')[0]  # Adjusting extraction method
            if pdbid in non_rrna_pdbids:
                src_path = os.path.join(sub_input_dir, filename)
                dest_path = os.path.join(sub_output_dir, filename)
                shutil.copy(src_path, dest_path)

def main():
    input_csv = 'hariboss_clean_20240621_error_revised_non_rnasm_removed_rrna_removed.csv'
    input_dir = 'hariboss_20240621_error_revised_result'
    sub_dirs = [
        'parsed_ligands_cif', 
        'parsed_ligands_pdb', 
        'parsed_ligands_pdb_add_H', 
        'parsed_ligands_pdb_remove_H', 
        'parsed_rna3db_rnas', 
        'parsed_rna3db_rnas_pdb', 
        'parsed_rna3db_rnas_pdb_add_H', 
        'parsed_rna3db_rnas_pdb_remove_H'
    ]
    output_dir = 'hariboss_20240621_error_revised_result_non_rnasm_removed_rrna_removed'

    filter_non_rrna_files(input_csv, input_dir, sub_dirs, output_dir)

if __name__ == "__main__":
    main()
