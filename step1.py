import sys

sys.path.append("..")
from rlclip.preprocessor import ComplexesDatabase
import pandas as pd


if __name__ == "__main__":
    hariboss_pd = pd.read_csv("./hariboss_clean.csv")
    # get pdbid list
    pdb_list = hariboss_pd["pdbid"].tolist()
    pdb_list = pdb_list[:5]

    # unique pdbid
    pdb_list = list(set(pdb_list))

    cd = ComplexesDatabase()
    # cd.clear_database()

    cd.pdb_batch_download(pdb_list, if_download=True)
    # this step is necessary for separating the complexes into rna and ligand
    # and parse RNA sequence and ligand SMILES
    # cd.parse_complex_database(filter_list=['HOH', 'MG', 'ALA', 'ARG', 'ASN', 'ASP'])
    cd.parse_complex_database(
        filter_list=["HOH"],
        rna_count_threshold=-1,
        distance_threshold=10.0,
        parse_id_list=pdb_list,
    )

    # build the dataset and save as npy files
    cd.build_dataset()

    print(pdb_list)