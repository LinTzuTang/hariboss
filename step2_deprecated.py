import pandas as pd
import numpy as np
from Bio import PDB
import matplotlib.pyplot as plt
import os
import warnings
import pickle
import tqdm
import sys

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
from prody import parsePDB, writePDB, parseMMCIF
from rdkit import Chem
from rdkit.Chem import AllChem

from rlclip.preprocessor.pdb_decoder.util import RLCLIPLogger


def extract_small_molecules_from_cif(cif_file, residue_name, complex_id, ligand_storage_path):
    parser = MMCIFParser()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure("id", cif_file)

    matching_residues = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0].strip() != "" and residue.get_resname() == residue_name:
                    matching_residues.append(residue)
                    # get the residue number
                    residue_number = residue.id[1]
                    # save the ligand file
                    ligand_file_name = complex_id + '_' + residue_name + '_' + str(residue_number)
                    io = PDB.PDBIO()
                    io.set_structure(residue)
                    io.save(ligand_storage_path + ligand_file_name + '.pdb')

    return matching_residues


def extract_small_molecules_with_prody(cif_file, residue_name, ligand_storage_path, ligand_chain_id, pdb_id, logger):
    prody_instance = parsePDB(cif_file)
    unique_resids = set(prody_instance.getResnums())

    qualified_small_molecules = []

    logger.write(f'Processing {pdb_id}\n')

    for resid in unique_resids:
        if resid < 0:
            logger.write(f'Negative resid {resid} found in {pdb_id}\n')
            continue

        small_molecule = prody_instance.select(f'resname {residue_name} and resnum {resid}')
        if small_molecule:
            qualified_small_molecules.append((small_molecule, resid))

    saved_ligands_path = []

    for small_molecule, resid in qualified_small_molecules:
        residue_name = small_molecule.getResnames()[0]
        file_path = f'{ligand_storage_path}{pdb_id}_{ligand_chain_id}_{residue_name}_{resid}.pdb'
        writePDB(file_path, small_molecule)
        saved_ligands_path.append((residue_name + '_' + str(resid), file_path))

    return saved_ligands_path


def extract_rna_chain_from_cif(cif_file, chain_id, logger):
    parser = MMCIFParser()
    structure = None

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", PDBConstructionWarning)
        try:
            structure = parser.get_structure("RNA", cif_file)
        except Exception as e:
            logger.write(f"Error occurred: {e}\n")
            return None

        if w:
            for warning in w:
                logger.write(f"Warning: {warning.message}\n")

    if structure is None:
        return None

    try:
        chain = structure[0][chain_id]
    except KeyError:
        logger.write(f"Chain ID {chain_id} not found in the structure.\n")
        return None

    for residue in list(chain):
        if residue.get_resname() not in ['A', 'U', 'C', 'G']:
            chain.detach_child(residue.id)
    return chain


MAX_BLOCKS = 50


def calculate_block_size(rna_length, max_blocks=MAX_BLOCKS):
    if rna_length <= max_blocks:
        return 1
    return (rna_length + max_blocks - 1) // max_blocks


def is_hydrogen(atom):
    return atom.name.startswith('H')


def blockwise_min_distance(rna_chain, ligand, block_size):
    num_rna_blocks = (len(rna_chain) + block_size - 1) // block_size
    num_ligand_atoms = len(ligand)

    distance_matrix = np.full((num_rna_blocks, num_ligand_atoms), float('inf'))
    original_distance_matrix = np.full((len(rna_chain), num_ligand_atoms), float('inf'))

    ligand_name = ligand.get_resname()
    rna_residues = list(rna_chain)

    for i, rna_residue in enumerate(rna_residues):
        rna_residue_name = rna_residue.get_resname()
        if rna_residue_name == ligand_name:
            continue
        for j, ligand_atom in enumerate(ligand):
            min_distance_in_block = float('inf')
            for rna_atom in rna_residue:
                if is_hydrogen(rna_atom):
                    continue
                dist = rna_atom - ligand_atom
                if dist < original_distance_matrix[i, j]:
                    original_distance_matrix[i, j] = dist
                if dist < min_distance_in_block:
                    min_distance_in_block = dist
            block_idx = i // block_size
            if min_distance_in_block < distance_matrix[block_idx, j]:
                distance_matrix[block_idx, j] = min_distance_in_block

    return distance_matrix, original_distance_matrix


def standardize_smiles(smiles_string, logger):
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        logger.write(f"Invalid SMILES string: {smiles_string}\n")
        return None
    return Chem.MolToSmiles(mol)


def generate_contact_map(distance_matrix, threshold=8):
    return (distance_matrix < threshold).astype(int)


def contact_map_generator(contact_threshold=8,
                          filter_rna_length=5,
                          file_path='./processed_data.csv',
                          logger=None):
    # dict that save all the information
    save_dict = {}

    og_hariboss_data = pd.read_csv(file_path,
                                   dtype={'pdbid': str, 'ligand_id': str, 'ligand_resnum': int, 'rna_chain': str,
                                          'ligand_chain': str})

    ############################################################
    # possible adjustment for hariboss data
    select_list = ['1eht', '1aju', '1am0', '1arj', '1akx']
    og_hariboss_data = og_hariboss_data[og_hariboss_data['pdbid'].isin(select_list)]
    og_hariboss_data = og_hariboss_data.reset_index(drop=True)
    ############################################################

    # add one more row to store contact map index
    og_hariboss_data[f'non_contact_map_index_{contact_threshold}A'] = np.nan
    # add contact percentage column
    og_hariboss_data[f'non_contact_percentage_{contact_threshold}A'] = np.nan
    # add canonical native smiles column
    og_hariboss_data['pipeline_canonical_smiles'] = np.nan

    new_processed_data = pd.DataFrame(columns=['pdbid', 'rna_chain', 'rna_chain_length',
                                               'ligand_chain', 'ligand_id', 'ligand_resnum', 'minimal_distance',
                                               f'non_contact_map_index_{contact_threshold}A', 'non_contact_percentage',
                                               'pipeline_canonical_smiles'])
    missing_data = pd.DataFrame(columns=['pdbid', 'rna_chain', 'rna_chain_length',
                                            'ligand_chain', 'ligand_id', 'ligand_resnum'])
    mismatch_resnum_data = pd.DataFrame(columns=['pdbid', 'rna_chain', 'rna_chain_length',
                                            'ligand_chain', 'ligand_id', 'ligand_resnum'])
    error_parse_data = pd.DataFrame(columns=['pdbid', 'rna_chain', 'rna_chain_length',
                                            'ligand_chain', 'ligand_id', 'ligand_resnum','error_type'])

    num_non_contact = 0
    parser = PDB.PDBParser()
    # using tqdm to show the progress bar
    for i in tqdm.tqdm(range(len(og_hariboss_data))):
        logger.write('\n---------------------------------------\n')
        logger.write(f'Length of og_hariboss_data: {len(og_hariboss_data)}\n')
        complex_id = str(og_hariboss_data.loc[i, 'pdbid'])
        rna_chain_id = str(og_hariboss_data.loc[i, 'rna_chain'])
        ligand_chain_id = str(og_hariboss_data.loc[i, 'ligand_chain'])
        ligand_id = str(og_hariboss_data.loc[i, 'ligand_id'])
        lignd_resnum = str(og_hariboss_data.loc[i, 'ligand_resnum'])

        # convert '6.00E+81' to '6e81'  # TODO: 6e81 downloaded in step 1 but seems failed, should be addressed
                                        # TODO: including more 6.00E+84, 6.00E+82,
        if '.00E+' in complex_id:
            # continue
            complex_id = complex_id.replace('.00E+', 'e')

        # in case the complex name is deprecated
        cif_path = "./complex_database/complexes_mmcif/" + complex_id + ".cif"

        logger.write(f"Complex ID: {complex_id}, RNA Chain ID: {rna_chain_id}, "
                     f"Ligand Chain ID: {ligand_chain_id}, Ligand ID: {ligand_id}, Ligand Resnum: {lignd_resnum}\n")
        pocket_name = f'{complex_id}_{rna_chain_id}_{ligand_chain_id}_{ligand_id}_{lignd_resnum}'

        if complex_id not in save_dict:
            save_dict[complex_id] = {}
        save_dict[complex_id][pocket_name] = {}
        save_dict[complex_id][pocket_name]['rna_chain_id'] = rna_chain_id
        save_dict[complex_id][pocket_name]['ligand_chain_id'] = ligand_chain_id
        save_dict[complex_id][pocket_name]['ligand_id'] = ligand_id
        save_dict[complex_id][pocket_name]['ligand_resnum'] = lignd_resnum
        save_dict[complex_id][pocket_name]['complex_id'] = complex_id

        if lignd_resnum == 'nan':
            lignd_resnum = '0'
            logger.write('No ligand residue number, set to 0\n')

        ligand_storage_path = './complex_database/ligands/'
        if ligand_chain_id == 'nan':
            ligand_chain_id = 'NA'

        ligand_file_name = complex_id + '_' + ligand_chain_id + '_' + ligand_id
        # find the ligand file in the ligand directory where the ligand_file_name is included in the file name
        ligand_file_list = [f for f in os.listdir(ligand_storage_path) if ligand_file_name in f]

        # there is also a proble that resnum and chain id are not match TODO: need to be addressed
        if len(ligand_file_list) == 0:
            ligand_file_list = [f for f in os.listdir(ligand_storage_path) if complex_id in f and ligand_id in f]
            if len(ligand_file_list) == 1:
                logger.write(
                    f"Ligand {ligand_id} in Complex {complex_id} is not match the resnum and chain id in the CIF File, "
                    f"but there is only one ligand file, so use it: {ligand_file_list[0]}\n")

        if len(ligand_file_list) == 0:
            logger.write(f"Ligand {ligand_file_name} in Complex {complex_id} is Missing, other information: "
                         f"RNA Chain ID: {rna_chain_id}, Ligand Chain ID: {ligand_chain_id}, Ligand ID: {ligand_id}, "
                         f"ligand_resnum: {lignd_resnum}\n")
            missing_data = missing_data.append({'pdbid': complex_id, 'rna_chain': rna_chain_id,
                                                'ligand_id': ligand_id,
                                                'ligand_resnum': lignd_resnum,
                                                'ligand_chain': ligand_chain_id,},
                                                ignore_index=True)
            missing_data.to_csv('./hariboss_complexes_missing_ligand_chain.csv', index=False)
            print('*******************************************')
            continue

        if len(ligand_file_list) > 1:
            logger.write(f"Ligand {ligand_file_name} in Complex {complex_id} is Ambiguous\n")
            logger.write(f"Ligand File List: {ligand_file_list}\n")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Discontinuous chains warning
            rna_chain = extract_rna_chain_from_cif(cif_path, rna_chain_id, logger)
            # save_dict[complex_id]['rna_chain_structure'] = rna_chain
            save_dict[complex_id][pocket_name]['rna_chain_file'] = cif_path.split('/')[-1]
            if rna_chain is None:
                logger.write(f"======= RNA Chain {rna_chain_id} in Complex {complex_id} encountered an error of parsing =======\n")
                error_parse_data = error_parse_data.append({'pdbid': complex_id,
                                                            'rna_chain': rna_chain_id,
                                                            'ligand_id': ligand_id,
                                                            'ligand_resnum': lignd_resnum,
                                                            'ligand_chain': ligand_chain_id,
                                                            'error_type': 'RNA Chain cannot be parsed'},
                                                            ignore_index=True)
                error_parse_data.to_csv('./hariboss_complexes_error_parse.csv', index=False)
                continue
            # num_rna_residues = len(rna_chain)
            logger.write(f"RNA Chain ID: {rna_chain_id}, RNA Chain Length: {len(rna_chain)}\n")

            if len(rna_chain) < filter_rna_length:
                logger.write(f"RNA Chain {rna_chain_id} in Complex {complex_id} is Too Short ({len(rna_chain)} < {filter_rna_length})\n")

            if len(rna_chain) == 0:
                logger.write(f"RNA Chain {rna_chain_id} in Complex {complex_id} has no RNA Residues, leading to an error of parsing\n")
                error_parse_data = error_parse_data.append({'pdbid': complex_id,
                                                            'rna_chain': rna_chain_id,
                                                            'ligand_id': ligand_id,
                                                            'ligand_resnum': lignd_resnum,
                                                            'ligand_chain': ligand_chain_id,
                                                            'error_type': 'RNA Chain has no RNA Residues'},
                                                            ignore_index=True)
                error_parse_data.to_csv('./hariboss_complexes_error_parse.csv', index=False)
                continue

        # get j from ligand_file_list where residue number is ligand-residue_num is in the file name
        j = 0
        if_match = False
        for k in range(len(ligand_file_list)):
            if lignd_resnum in ligand_file_list[k]:
                j = k
                if_match = True
                break
        if not if_match:
            """
            Mis match resnum part
            """
            if len(ligand_file_list) == 1:
                j = 0
                logger.write(f"Ligand {ligand_id} in Complex {complex_id} is not match the resnum in the CIF File, "
                             f"but there is only one ligand file, so use it: {ligand_file_list[j]}\n")
            else:
                logger.write(f"Ligand {ligand_id} in Complex {complex_id} is not match the resnum in the CIF File, other information:"
                             f"RNA Chain ID: {rna_chain_id}, Ligand Chain ID: {ligand_chain_id}, Ligand ID: {ligand_id}, "
                             f"ligand_resnum: {lignd_resnum}\n")
                mismatch_resnum_data = mismatch_resnum_data.append({'pdbid': complex_id, 'rna_chain': rna_chain_id,
                                                    'ligand_id': ligand_id,
                                                    'ligand_resnum': lignd_resnum,
                                                    'ligand_chain': ligand_chain_id,},
                                                    ignore_index=True)
                mismatch_resnum_data.to_csv('./hariboss_complexes_mismatch_resnum.csv', index=False)
                print('===========================================')
                continue

        ligand_file = ligand_storage_path + ligand_file_list[j]
        logger.write(f"Loaded Ligand File: {ligand_file}\n")

        # read structure of rna and ligand
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                ligand_structure = parser.get_structure(ligand_id, ligand_file)
                # save_dict[complex_id]['ligand_structure'] = ligand_structure
                save_dict[complex_id][pocket_name]['ligand_file'] = ligand_file.split('/')[-1]
            except Exception as e:
                logger.write(f"======= Ligand {ligand_id} in Complex {complex_id} encountered an error of parsing =======\n")
                error_parse_data = error_parse_data.append({'pdbid': complex_id,
                                                            'rna_chain': rna_chain_id,
                                                            'ligand_id': ligand_id,
                                                            'ligand_resnum': lignd_resnum,
                                                            'ligand_chain': ligand_chain_id,
                                                            'error_type': 'Ligand cannot be parsed'},
                                                            ignore_index=True)
                error_parse_data.to_csv('./hariboss_complexes_error_parse.csv', index=False)
                logger.write(f"Error occurred: {e}\n")
                continue
            ligand = list(ligand_structure.get_atoms())[0].parent
            # get SMILES from ligand file
            ligand_smiles = Chem.MolToSmiles(Chem.MolFromPDBFile(ligand_file))
            # standardize the SMILES
            standard_ligand_smiles = standardize_smiles(ligand_smiles, logger)

        block_size = calculate_block_size(len(rna_chain))
        distance_matrix, original_distance_matrix = blockwise_min_distance(rna_chain, ligand, block_size)
        # get the minimal distance
        min_distance = np.min(original_distance_matrix)

        # switch original distance matrix to contact map
        contact_map_org = generate_contact_map(original_distance_matrix, contact_threshold)
        # compress contact map towards y-axis (to 1D), if the row is all 0, then assign 0 to the place, else assign 1
        contact_map = np.where(np.sum(contact_map_org, axis=1) == 0, 0, 1)
        # get the index of zero elements
        zero_index = np.where(contact_map == 0)[0]
        # encode zero_index ndarray to string
        zero_index = ','.join(map(str, zero_index))
        # store the contact map index
        og_hariboss_data.loc[i, f'non_contact_map_index_{contact_threshold}A'] = zero_index
        # store the contact percentage
        og_hariboss_data.loc[i, f'non_contact_percentage_{contact_threshold}A'] = np.sum(zero_index == 0) / len(zero_index)
        if np.sum(contact_map == 0) / len(contact_map) == 1:
            logger.write(f'There is no contact between RNA and Ligand in this pair. {complex_id} {rna_chain_id} {ligand_id} {lignd_resnum}\n')
            num_non_contact += 1
            new_processed_data = new_processed_data.append({'pdbid': complex_id, 'rna_chain': rna_chain_id,
                                                            'rna_chain_length': len(rna_chain),
                                                            'ligand_id': ligand_id, 'minimal_distance': min_distance,
                                                            f'non_contact_map_index_{contact_threshold}A': zero_index,
                                                            'non_contact_percentage': np.sum(contact_map == 0) / len(contact_map),
                                                            'pipeline_canonical_smiles': standard_ligand_smiles,
                                                            'ligand_resnum': lignd_resnum,
                                                            'ligand_chain': ligand_chain_id,},
                                                            ignore_index=True)
            new_processed_data.to_csv('./processed_hariboss.csv', index=False)
            save_dict[complex_id][pocket_name]['contact_map_2d'] = contact_map_org
            save_dict[complex_id][pocket_name]['contact_map_1d'] = contact_map
            save_dict[complex_id][pocket_name]['contact_percentage'] = np.sum(contact_map == 1) / len(contact_map)
            save_dict[complex_id][pocket_name]['canonical_smiles'] = standard_ligand_smiles
            save_dict[complex_id][pocket_name]['minimal_distance'] = min_distance
            save_dict[complex_id][pocket_name]['rna_length'] = len(rna_chain)
            save_dict[complex_id][pocket_name]['rna_chain_sequence'] = ''.join([residue.get_resname() for residue in rna_chain])
            save_dict[complex_id][pocket_name]['ligand_SMILES'] = ligand_smiles
        else:
            new_processed_data = new_processed_data.append({'pdbid': complex_id, 'rna_chain': rna_chain_id,
                                                            'rna_chain_length': len(rna_chain),
                                                            'ligand_id': ligand_id, 'minimal_distance': min_distance,
                                                            f'non_contact_map_index_{contact_threshold}A': zero_index,
                                                            'non_contact_percentage': np.sum(contact_map == 0) / len(contact_map),
                                                            'pipeline_canonical_smiles': standard_ligand_smiles,
                                                            'ligand_resnum': lignd_resnum,
                                                            'ligand_chain': ligand_chain_id},
                                                            ignore_index=True)
            # save
            new_processed_data.to_csv('./processed_hariboss.csv', index=False)
            save_dict[complex_id][pocket_name]['contact_map_2d'] = contact_map_org
            save_dict[complex_id][pocket_name]['contact_map_1d'] = contact_map
            save_dict[complex_id][pocket_name]['contact_percentage'] = np.sum(contact_map == 1) / len(contact_map)
            save_dict[complex_id][pocket_name]['canonical_smiles'] = standard_ligand_smiles
            save_dict[complex_id][pocket_name]['minimal_distance'] = min_distance
            save_dict[complex_id][pocket_name]['rna_length'] = len(rna_chain)
            save_dict[complex_id][pocket_name]['rna_chain_sequence'] = ''.join([residue.get_resname() for residue in rna_chain])
            save_dict[complex_id][pocket_name]['ligand_SMILES'] = ligand_smiles

            # save the dictionary
            with open('./hariboss_processed_dict.pkl', 'wb') as f:
                pickle.dump(save_dict, f)

    # save the dictionary
    with open('./hariboss_processed_dict.pkl', 'wb') as f:
        pickle.dump(save_dict, f)

    og_hariboss_data.to_csv('./hariboss_complexes_exception.csv', index=False)
    logger.write(f'\nProcessed {len(og_hariboss_data)} pairs finished!\n')
    logger.write(f'Number of pairs with no contact: {num_non_contact} in total {len(og_hariboss_data)} pairs')


if __name__ == '__main__':

    logger = RLCLIPLogger('./data_build_log')

    contact_map_generator(contact_threshold=10, filter_rna_length=5, file_path='./hariboss_clean.csv',
                          logger=logger)
