# RL-CLIP Data Collection and Processing Pipeline

This is a un-refactor version of the pipeline used to collect and process data for the RL-CLIP.

## 1. Clean Hariboss data
#### 1. Download Hariboss then clean 
You can download the original and up-to-date HARIBOSS dataset in CSV format from the [HARIBOSS website](https://hariboss.pasteur.cloud/complexes/). To generate the new `hariboss_clean.csv` file, run the `clean_hariboss.py` script as follows:

```bash
python clean_hariboss.py -i Complexes.csv -o ./hariboss_clean_20240621.csv
```

#### 2. Revise error records 
We found that there are some error records listed in error_record.csv. These errors have been manually corrected and are listed in error_record_revised.csv.

The `correct_hariboss.py` script is to replace the error records with the revised ones. To run the script, use the following command:

```bash
python correct_hariboss.py --input_file hariboss_clean_20240621.csv --error_file error_record.csv --revised_file error_record_revised.csv --output_file hariboss_clean_20240621_error_revised.csv
```
#### 3. Filter Non-RNA small molecular complexes and ribosomal RNA
The script `filter_hariboss_non_rnasm_rrna.py` filters out non-RNA SM complexes (`non_RNA_SM_Complexes.csv`) and optionally ribosomal RNA (`Ribosomal_Complexes.csv`) from the Hariboss dataset. It can also copy the filtered CIF files to a new directory (if you already downloaded all complexes mmCIF files).

**Usuage**
```bash
python filter_hariboss_non_rnasm_rrna.py <hariboss_csv> <non_rnasm_csv> <output_csv> [--rrna_csv <rrna_csv>] [--input_dir <input_dir>] [--output_dir <output_dir>] <filter_option> [--copy_files]
```
**Arguments**
* hariboss_csv: Path to the Hariboss CSV file.
* non_rnasm_csv: Path to the non-RNA SM complexes CSV file.
* output_csv: Path to the output CSV file where filtered data will be saved.
* --rrna_csv: Path to the ribosomal complexes CSV file (required if filter_option is non_rnasm_rrna).
* --input_dir: Directory containing the original CIF files (required if --copy_files is specified).
* --output_dir: Directory to save the filtered CIF files (required if --copy_files is specified).
* filter_option: Filter option: non_rnasm or non_rnasm_rrna.
* --copy_files: Option to copy filtered CIF files to the output directory.

**Examples**
For filtering only non-RNA SM complexes:
```bash
python filter_hariboss_non_rnasm_rrna.py hariboss_clean_20240621_error_revised.csv non_RNA_SM_complexes.csv hariboss_clean_20240621_error_revised_non_rnasm_removed.csv non_rnasm
```
For filtering both non-RNA SM and ribosomal complexes:
```bash
python filter_hariboss_non_rnasm_rrna.py hariboss_clean_20240621_error_revised.csv non_RNA_SM_complexes.csv --rrna_csv Ribosomal_Complexes.csv hariboss_clean_20240621_error_revised_non_rnasm_removed_rrna_removed.csv non_rnasm_rrna
```
For filtering both non-RNA SM and ribosomal complexes and copying the filtered CIF files to a new directory:
```bash
python filter_hariboss_non_rnasm_rrna.py hariboss_clean_20240621_error_revised.csv non_RNA_SM_complexes.csv --rrna_csv Ribosomal_Complexes.csv hariboss_clean_20240621_error_revised_non_rnasm_removed_rrna_removed.csv --input_dir complex_mmcif_database --output_dir complex_mmcif_database_filtered non_rnasm_rrna --copy_files
```

## 2. Download mmCIF files

The `download_mmcif.py` script allows you to download mmCIF files for a list of PDB IDs from cleaned hariboss data(For example: `hariboss_clean.csv`). You can download all PDB IDs, limit the number of PDB IDs, or specify a list of PDB IDs.

#### To download all PDB IDs in the `hariboss_clean.csv`:

```bash
python download_mmcif.py ./hariboss_clean.csv --output_dir complex_mmcif_database
```

#### To limit the number of PDB IDs:

```bash
python download_mmcif.py ./hariboss_clean.csv --output_dir complex_mmcif_database_filtered --num_pdb 5
```

#### To specify a list of PDB IDs:

```bash
python download_mmcif.py ./hariboss_clean.csv --output_dir complex_mmcif_database_filtered --pdb_id_list 1aju 6ymj
```

#### To download all PDB IDs in the `hariboss_clean_20240621_error_revised_non_rnasm_removed_rrna_removed.csv`:

```bash
python download_mmcif.py ./hariboss_clean_20240621_error_revised_non_rnasm_removed_rrna_removed.csv --output_dir complex_mmcif_database_non_rnasm_removed_rrna_removed
```

## 3. RNA Liagnads split and Parsing

This main script is for processing of RNA ligand data by downloading necessary files and generating datasets with curated pocket information. It leverages the 'MmcifDownloader' and 'DatasetGenerator' modules to handle downloading and data preprocessing.

This step requires curated pocket information provided in a CSV file with the following columns:

- PDB ID
- Ligand ID
- Ligand Chain
- RNA Chain
- Ligand ResNum

The example file is `hariboss_clean.csv`. For RNA chain selection and the decoding algorithm in the code, we use the author-provided chain, while for ligand chain selection, we use the segment chain.

In this step, it will save a dictionary of the pocket information with every needed for the model training.
An example result directory is saved in `example_result`. </br>

### Parameter Configuration

#### (1) Download CIF Files by MmcifDownloader

To download the CIF files, the MmcifDownloader module requires the following parameters:

- csv_file: Path to the CSV file containing PDB IDs and associated information.
- output_dir: Directory where the downloaded CIF files will be stored.
- num_pdb (optional): Number of PDB files to download. If not specified, all files in the CSV will be downloaded.
- pdb_id_list (optional): List of specific PDB IDs to download. If not specified, all files in the CSV will be downloaded.

##### Example configuration:

```
csv_file = "hariboss_clean.csv"
output_dir = "complex_database/complexes_mmcif"
num_pdb = 5 (downloads 5 PDB files)
pdb_id_list = ["5hbw", "6m0j"] (downloads specified PDB files)
```

#### (2) RNA Ligands Split and RNA3DB Parser by DatasetGenerator

The DatasetGenerator module processes the downloaded CIF files and generates datasets. It requires the following parameters:

- dataset_path: Path to the CSV file containing dataset information.
- result_path: Directory where the processed results will be stored.
- mmcif_path: Directory containing the downloaded CIF files.
- contact_threshold (optional): Threshold for contact distance. Default is 10.
- filter_rna_length (optional): Minimum length of RNA to consider. Default is 5.
- fp_type_list (optional): List of fingerprint types to use. Default is ["fp2", "morgan"].
- molecule_filter_list (optional): List of molecule types to filter out. Default is ["HOH"].
- save_file (optional): File format for saving results. Default is "cif".
- use_rna3db (optional): Boolean indicating whether to use RNA3DB parser. Default is True.

##### Example configuration:

```
dataset_path = "hariboss_clean.csv"
result_path = "./example_result"
mmcif_path = "complex_database/complexes_mmcif"
contact_threshold = 10
filter_rna_length = 5
fp_type_list = ["fp2", "morgan"]
molecule_filter_list = ["HOH"]
save_file = "cif"
use_rna3db = True
```

#### Error Handling

While processing, the script will generate error files in the form of txt and csv files to indicate specific issues encountered.

#### Running the Script

```bash
python hariboss_process.py
```

