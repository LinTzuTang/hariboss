import sys

sys.path.append("..")
from download_mmcif import MmcifDownloader
from rlclip.preprocessor import DatasetGenerator

import pandas as pd


def main(
    hariboss_file,
    mmcif_path,
    result_path,
    # num_pdb=None,
    # pdb_list=None,
    contact_threshold=10,
    # filter_rna_length=5,
    fp_type_list=None,
):

    if fp_type_list is None:
        fp_type_list = ["fp2", "morgan"]

    # 1. Download cif files by MmcifDownloader
    mmcifdownloader = MmcifDownloader(
        csv_file=hariboss_file,
        output_dir=mmcif_path,
        # num_pdb=num_pdb,
        # pdb_id_list=pdb_list,
    )
    #mmcifdownloader.clean_output_files()
    mmcifdownloader.download_mmcif()

    # 2. RNA Ligands split and RNA3DB parser by DatasetGenerator
    dg = DatasetGenerator(
        dataset_path=hariboss_file,
        result_path=result_path,
        mmcif_path=mmcif_path,
        contact_threshold=contact_threshold,
        # filter_rna_length=filter_rna_length,
        fp_type_list=fp_type_list
    )

    # # read the hariboss dataset
    # hariboss_pd = pd.read_csv("hariboss_clean.csv")
    # # get pdbid list
    # pdb_ids = hariboss_pd["pdbid"].tolist()
    # pdb_ids = pdb_ids[:5]
    # # unique pdbid
    # pdb_ids = list(set(pdb_ids))

    dg.clean_result_path()
    # dg.select_pdb(['4w2i'])
    dg.generate_dataset()


if __name__ == "__main__":
    #pdb_list = ["5hbw"]
    hariboss_file = 'hariboss_clean_20240621_error_revised.csv'
    mmcif_path = "./complex_mmcif_database"
    #num_pdb = 5  # number of pdb files to download
    result_path = "./hariboss_20240621_error_revised_result"

    main(
        hariboss_file,
        mmcif_path,
        result_path
    )


