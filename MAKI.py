'''
This is MAKI ฅ≽(•⩊ •マ≼  It's a ( Molecule Against Kinase Inspector) tool. It's an attempt to do massive inverse screening and present some visualisation of the results.
Written by: Sophie Bernard(-Doucet)
Produced within the Najmanovich Research Group (May 2025 -> August 2025).
'''
import argparse
import os
import re
import shutil
from pathlib import Path
from getcleftpy import run_getcleft
from nrgrank import process_target, process_ligands, nrgrank_main, generate_conformers
from tqdm import tqdm
import subprocess

# def __list_files_cleft(input_dir):
#     file_list = os.listdir(input_dir)
#     #id_list = [file.split('_')[1] for file in file_list]
#     # id_list = [".".join(file.split("-")[:2]) for file in file_list] #this pattern adaptation is for AF notation.
#     return id_list, file_list

def obabel_translate(input_dir, output_path):
    output_path.mkdir(parents=False, exist_ok=True)
    pdb_files=list(input_dir.glob("*.pdb"))
    for pdb_file in tqdm(pdb_files):
        mol2_nom =pdb_file.stem+".mol2"
        output_file = output_path / mol2_nom
        try:
            result = subprocess.run(
                ["obabel", "-ipdb", str(pdb_file), "-omol2", "-O", str(output_file)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if result.returncode != 0:
                print(f"obabel failed with {result.returncode}")
            else:
                continue
        except FileNotFoundError:
            print(f"File for the conversion not found at {pdb_file}")
            break
    return output_path

def getcleft(input_dir, output_clefts_proteome, verbose=False, delete_clfs=True ):
    os.makedirs(output_clefts_proteome, exist_ok=True)
    file_list = os.listdir(input_dir)
    for index, pdb_path in enumerate(file_list):
        pdb_id = pdb_path.split('_')[1]
        verbose and print(f"Structure n* : {index+1} of {len(file_list)}")
        pdb_path = os.path.join(input_dir, pdb_path)
        if os.path.isfile(pdb_path):
            try:
                output_dir_id = os.path.join(output_clefts_proteome, pdb_id)
                os.makedirs(output_dir_id, exist_ok=True)
                output_base_path = os.path.join(output_dir_id, f"{pdb_id}")
                run_getcleft(pdb_path, output_base=output_base_path)
                delete_clfs and remove_clfs(output_clefts_proteome, verbose=verbose)
            except Exception as e:
                verbose and print(f"≽(•⩊ •マ≼ GetCleft failed for {pdb_path}: {e} :(")
        else:
            verbose and print(f"≽(•⩊ •マ≼  ; file {pdb_path} for getcleft cannot be found, find it >:(")

def remove_clfs(clefts_id_folder, verbose=False):
    include_pattern = re.compile(r'_clf_([1-9][0-9]{0,2})\.pdb$', re.IGNORECASE)
    for root, _, files in os.walk(clefts_id_folder):
        for filename in files:
            if include_pattern.search(filename):
                source_path = os.path.join(root, filename)
                verbose and print(f"≽(•⩊ •マ≼: {source_path} deleted")
                try:
                    os.remove(source_path)
                except Exception as e:
                    verbose and print(f"Couldn't download ; {source_path}: {e}")
            else:
                continue

def nrgrank(ligands_dir,root_cleft_dir, output_dir,verbose=True):
    for ligand_file in os.listdir(ligands_dir):
        if ligand_file.endswith(".mol2"):
            conformers_per_molecule=1
            if not ligand_file:
                verbose and print("catastrophic failure,no mol2 files here, find the drugs ! ")
            ligand_name = Path(ligand_file).stem
            ligand_save_dir = process_ligands(
                ligand_path=os.path.join(ligands_dir,ligand_file),
                conformers_per_molecule=conformers_per_molecule,
                output_dir=os.path.join(ligands_dir,ligand_name),
                overwrite=False)
        else:
            continue

        #todo : arg/fix this, its not user friendly at all but here we go :
        mol2_kinases_path =(r"/lustre09/project/6049264/sbernarr/MAKI/AF2_pdb_to_mol2_kinases")
        mol2_kinases = sorted(next(os.walk(mol2_kinases_path), (None, None, []))[2])

        for kin_target in mol2_kinases[403:]:
            if kin_target.endswith(".mol2"):
                kinase_name = (kin_target.split("_")[1])
                kinase_dir = next((d for d in Path(root_cleft_dir).iterdir() if d.is_dir() and d.name == kinase_name), None)
                files_in_dir = [f for f in Path(kinase_dir).iterdir() if f.is_file()]
                bdsite_file_path = files_in_dir[0]

                target_save_dir = process_target(
                    target_mol2_path=str(os.path.join(mol2_kinases_path,kin_target)),
                    binding_site_file_path=str(bdsite_file_path),
                    create_new_dir=True,
                    overwrite=False)


                if target_save_dir:
                    result_file_path = nrgrank_main(target_name=kinase_name,
                                                    preprocessed_target_path=target_save_dir,
                                                    preprocessed_ligand_path=ligand_save_dir,
                                                    result_folder_path=os.path.join(output_dir,ligand_file,kin_target),
                                                    USE_CLASH=False,
                                                    POSES_SAVED_PER_MOLECULE=0)
                    verbose and print(f" ฅ≽(•⩊ •マฅ≼ NRGrank ran !  {result_file_path}")

if __name__ == '__main__':
    #id_list, file_list = __list_files_cleft(input_dir=r"C:\Users\sophi\Documents\NRG\KINASES\Kincore_AlphaFold2_ActiveHumanCatalyticKinases_v2\AF2_pdb_kinases")

    getcleft(input_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_pdb_kinases", output_clefts_proteome=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_kinase_clefts")
    nrgrank(ligands_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_ligands", root_cleft_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_kinase_clefts", output_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_nrgrank")



#------------------------------------------------------------------------- future argumentalisation
# # id_list, file_list, repoutput, flag poru verbose, flag for clfs
# parser = argparse.ArgumentParser(description="Maki pipeline to present the kinases against a chosen set of ligands using NRGRank and introduce further analysis.")
#
# parser.add_argument("--test_ligand_set", type=str, required=True,
#                     help="Path to repo containing all 10 test ligands (X 10  total)")
# parser.add_argument("--kin_targets_dir", type=str, required=True,
#                     help="Path to repo containing multiple .mol2 target files (154 total). Aims to do the treatment of one file / run")
# parser.add_argument("--pdb_kin_bdsite_folder", type=str, required=True,
#                     help="Path to repo containing PDB binding site files")
# args = parser.parse_args()
#
# # running it on narval
# test_ligand_set = Path(args.test_ligand_set)
# test_lig_files = list(test_ligand_set.glob("*.mol2"))
# if not test_lig_files:
#     raise FileNotFoundError("catastrophic failure, there are no ligands")
# test_lig_file = test_lig_files[0]
# current_ligand_name = test_lig_file.stem
#
# kin_targets_dir = Path(args.kin_targets_dir)
# kin_files = sorted(kin_targets_dir.glob("*.mol2"))
# if not kin_files:
#     raise FileNotFoundError("catastrophic failure, there are no mol2 files")
#
# pdb_kin_bdsite_folder = Path(args.pdb_kin_bdsite_folder)
# kin_bd_files = sorted(pdb_kin_bdsite_folder.glob("*.pdb"))
# if not kin_bd_files:
#     raise FileNotFoundError("catastrophic failure, there are no binding site files")