import pandas as pd
import argparse


def parse_ligand_info(row):
    pdbid = row["id"]
    ligand_info_list = eval(row["sm_ligand_ids"])

    parsed_data = []
    for ligand_info in ligand_info_list:
        # Split the ligand information
        ligand_parts = ligand_info.split(":")
        ligand_id = ligand_parts[0].split("_")[0]
        ligand_resnum_chain = ligand_parts[1]
        ligand_chain_and_resnum = ligand_resnum_chain.split("/")
        ligand_chain = ligand_chain_and_resnum[0]
        ligand_resnum = ligand_chain_and_resnum[1]
        # Ensure proper splitting of RNA chain
        rna_chain = ligand_parts[2]
        parsed_data.append([pdbid, ligand_id, ligand_resnum, rna_chain, ligand_chain])

    return parsed_data


def main(input_file, output_file):
    # Read the input CSV file
    df = pd.read_csv(input_file)
    # Select relevant columns
    df = df[["id", "sm_ligand_ids"]]
    # Apply the function to the dataframe and flatten the list of lists
    parsed_data = df.apply(parse_ligand_info, axis=1).explode().tolist()
    # Create a new dataframe from the parsed data
    new_df = pd.DataFrame(
        parsed_data,
        columns=["pdbid", "ligand_id", "ligand_resnum", "rna_chain", "ligand_chain"],
    )
    # Write the output CSV file
    new_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean HARIBOSS dataset.")
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Input CSV file path"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output CSV file path"
    )
    args = parser.parse_args()
    main(args.input, args.output)
