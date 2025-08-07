import sys

sys.path.append("..")

import pandas as pd
from rlclip.preprocessor import DatasetGenerator, ComplexesDatabase

if __name__ == "__main__":

    total_dir = "./complex_database_test"
    dataset_path = "./hariboss_clean_20240621.csv"
    result_path = "./example_result"
    # complex_database_path = "./complex_database"
    contact_threshold = 10
    filter_rna_length = 5
    fp_type_list = ["maccs", "fp2", "morgan"]
    use_rna3db = True
    write_rna3db_rna = True

    # read the hariboss dataset
    hariboss_pd = pd.read_csv(dataset_path)
    # get pdbid list
    pdb_list = hariboss_pd["pdbid"].tolist()
    pdb_list = pdb_list[:5]
    # unique pdbid
    pdb_list = list(set(pdb_list))

    # step1 download pdb files and parse the complex
    cd = ComplexesDatabase(total_dir=total_dir)
    cd.clear_database()

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

    # step2 dataset generation and RNA3DB praser and writer
    dg = DatasetGenerator(
        dataset_path=dataset_path,
        result_path=result_path,
        complex_database=total_dir,
        contact_threshold=contact_threshold,
        filter_rna_length=filter_rna_length,
        fp_type_list=fp_type_list,
        use_rna3db=use_rna3db,
        write_rna3db_rna=write_rna3db_rna,
    )

    # temporarily select 5 pdb files for example
    dg.clean_result_path()
    dg.select_pdb(pdb_list)
    dg.generate_dataset()
