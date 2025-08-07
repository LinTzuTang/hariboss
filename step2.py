import sys

sys.path.append("..")
from rlclip.preprocessor import DatasetGenerator

if __name__ == "__main__":
    dataset_path = "./hariboss_clean.csv"
    result_path = "./example_result"
    mmcif_path = "./complex_database/complexes_mmcif"
    contact_threshold = 10
    filter_rna_length = 5
    fp_type_list = ["fp2", "morgan"]
    molecule_filter_list = ['HOH']

    dg = DatasetGenerator(
        dataset_path=dataset_path,
        result_path=result_path,
        mmcif_path=mmcif_path,
        contact_threshold=contact_threshold,
        filter_rna_length=filter_rna_length,
        molecule_filter_list=molecule_filter_list,
        fp_type_list=fp_type_list,
    )

    # temporarily select 5 pdb files for example
    dg.select_pdb(["1aju", '6ymj'])
    dg.generate_dataset()
