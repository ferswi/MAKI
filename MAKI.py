
import argparse
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import requests
import subprocess
import time

from collections import defaultdict
from getcleftpy import run_getcleft
from nrgrank import process_target, process_ligands, nrgrank_main, generate_conformers
from pathlib import Path
from sklearn.metrics import roc_curve, auc
from tqdm import tqdm

#personal functions :
def __list_files_cleft(input_dir):
    file_list = os.listdir(input_dir)
    id_list = [file.split('_')[1] for file in file_list]
    # id_list = [".".join(file.split("-")[:2]) for file in file_list] #AF2 notation
    return id_list, file_list

def obabel_translate(input_dir, output_path):
    output_path.mkdir(parents=False, exist_ok=True)
    pdb_files = list(input_dir.glob("*.pdb"))
    for pdb_file in tqdm(pdb_files):
        mol2_nom = pdb_file.stem + ".mol2"
        output_file = output_path / mol2_nom
        try:
            result = subprocess.run(
                ["obabel", "-ipdb", str(pdb_file), "-omol2", "-O", str(output_file)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if result.returncode != 0:
                print(f"Obabel failed with {result.returncode}")
        except FileNotFoundError:
            raise (FileNotFoundError(f"File for the conversion not found : {pdb_file}"))
    return output_path

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

def find_kins_overlapping_dfs(empiric_df, experimental_df, verbose):
    empiric_kinases = set(empiric_df["Kinase"])

    experimental_df = experimental_df.reset_index()
    experimental_kinases = set(experimental_df["Kinase"])

    overlapping_kinases = empiric_kinases.intersection(experimental_kinases)
    verbose and print("≽(•⩊ •マ≼ Overlapping kinases found : ", overlapping_kinases)
    return sorted(list(overlapping_kinases))

def getcleft(input_dir, output_clefts_proteome, verbose=True, delete_clfs=True, seperator="_"):
    os.makedirs(output_clefts_proteome, exist_ok=True)
    file_list = os.listdir(input_dir)
    for index, pdb_path in enumerate(file_list):
        pdb_id = pdb_path.split(seperator)[1]
        verbose and print(f"Structure n* : {index + 1} of {len(file_list)}")
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
            verbose and print(f"≽(•⩊ •マ≼  ; file {pdb_path} for GetCleft cannot be found, find it >:(")

def nrgrank(ligands_dir, root_cleft_dir, output_dir, kinase_mol_dir, seperator="_", verbose=True):
    for ligand_file in os.listdir(ligands_dir):
        if ligand_file.endswith(".mol2"):
            conformers_per_molecule = 1
            if not ligand_file:
                ligand_name = Path(ligand_file).stem
                ligand_save_dir = process_ligands(
                    ligand_path=os.path.join(ligands_dir, ligand_file),
                    conformers_per_molecule=conformers_per_molecule,
                    output_dir=os.path.join(ligands_dir, ligand_name),
                    overwrite=False)
        else:
            continue

        mol2_kinases = sorted(next(os.walk(), (None, None, []))[2])
        for kin_target in mol2_kinases:
            if kin_target.endswith(".mol2"):
                kinase_name = (kin_target.split(seperator)[1])
                kinase_dir = next((d for d in Path(root_cleft_dir).iterdir() if d.is_dir() and d.name == kinase_name),
                                  None)  ## PRAISE: I really like this, this is clever
                files_in_dir = [f for f in Path(kinase_dir).iterdir() if f.is_file()]
                bdsite_file_path = files_in_dir[0]

                target_save_dir = process_target(
                    target_mol2_path=str(os.path.join(kinase_mol_dir, kin_target)),
                    binding_site_file_path=str(bdsite_file_path),
                    create_new_dir=True,
                    overwrite=True)

                if target_save_dir:
                    result_file_path = nrgrank_main(target_name=kinase_name,
                                                    preprocessed_target_path=target_save_dir,
                                                    preprocessed_ligand_path=ligand_save_dir,
                                                    result_folder_path=os.path.join(output_dir,
                                                                                    ligand_name.split(".mol2")[0],
                                                                                    kin_target),
                                                    USE_CLASH=False,
                                                    POSES_SAVED_PER_MOLECULE=0)
                    verbose and print(f" ฅ≽(•⩊ •マฅ≼ NRGrank ran !  {result_file_path}")

def creating_nrgrank_df(input_dir, verbose, output_path):
    data = defaultdict(dict)
    ligands = []
    for ligand in sorted(os.listdir(input_dir)):
        ligand_path = os.path.join(input_dir, ligand)
        ligands.append(ligand)
        verbose and print("≽(•⩊ •マ≼ Ligands found for MAKI: ", ligand)

        for kinase_repo in sorted(os.listdir(ligand_path)):
            kinase_repo_path = (os.path.join(ligand_path, kinase_repo))
            for kinase_file in sorted(os.listdir(kinase_repo_path)):
                if kinase_file.endswith(".csv"):
                    csv_file = os.path.join(kinase_repo_path, kinase_file)
                    df = pd.read_csv(csv_file, header=None)
                    nrg_value = df.iloc[1, 1]
                    data[kinase_file.strip(".csv")][ligand] = nrg_value

    nrgrank_df = pd.DataFrame.from_dict(data, orient='index')
    nrgrank_df.index.name = 'légume'
    nrgrank_df = nrgrank_df[ligands]
    nrgrank_df.to_csv(output_path)
    return nrgrank_df, ligands

def merging_dfs(nrgrank_df, davis_df, output_path, overlapping_kinases):
    ligands = nrgrank_df.columns.tolist()[1:]
    merged_df = pd.DataFrame({"Kinase": overlapping_kinases})
    davis_df = davis_df.fillna(20000) # todo : do i change this logic and keep the nan in both, or fill the nrgrank ones as wlel w 2000 ?

    nrg_index = 1
    dav_index = 2

    for ligand in ligands:
        nrgrank_sub = nrgrank_sub.rename(columns={ligand: f"{ligand}_nrg"}, inplace=False)
        davis_sub = davis_sub.rename(columns={ligand: f"{ligand}_dav"}, inplace=False)

        merged_df = pd.merge(merged_df, nrgrank_sub, on="Kinase", suffixes=("_x", "_dav"), how="inner")
        merged_df = pd.merge(merged_df, davis_sub, on="Kinase", suffixes=("_y", "_nrg"), how="inner")
        merged_df.to_csv(output_path, index=False)

        # rank organising per ligand
        nrg_sorted = merged_df.sort_values(by=merged_df.columns[nrg_index])
        nrgnames = nrg_sorted.iloc[:, 0].tolist()
        dav_sorted = merged_df.sort_values(by=merged_df.columns[dav_index])
        davnames = dav_sorted.iloc[:, 0].tolist()

        ranked_df = merged_df.copy()
        ranked_df['nrgnames_sorted'] = pd.Series(nrgnames)
        ranked_df['davnames_sorted'] = pd.Series(davnames)

        # skipping columns to maintain good ligand indexing.
        nrg_index += 2
        dav_index += 2
        ranked_df.to_csv(output_path, index=False)
    return ligands

def plotting(ligands,ranked_df,values_plotting, ranks_plotting, roc_curving, contingency_table):
    # ini figure of ROC curves
    roc_data = []
    colors = cm.tab20(np.linspace(0, 1, len(ligands)))
    plt.figure(figsize=(8, 8), dpi=200)
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')
    for ligand in ligands:
        #todo; should these lines be in the loop ?
        nrg = ranked_df.iloc[:, ranked_df.nrg_index]
        dav = ranked_df.iloc[:, ranked_df.davindex]

        # ROC curve implementation & ini of the roc figure, one figure w all ligands within
        davis_true = (dav < 20000).astype(int)

        nrg_inverted_scores = -nrg  # roc curve logic assumes bigger is more likely to be +.

        fpr, tpr, _ = roc_curve(davis_true, nrg_inverted_scores)
        roc_auc = auc(fpr, tpr)
        roc_data.append((ligand, fpr, tpr, roc_auc, colors[len(roc_data)]))
        plt.plot(fpr, tpr, color=colors[len(roc_data) - 1],
                 label=f"{ligand} (AUC = {roc_auc:.2f})")

        if values_plotting:
            plt.figure(figsize=(8, 8), dpi=200)
            plt.scatter(nrg, dav)
            for i, label in enumerate(ligands):  # todo : labels ?, ligands  ## CONCERN: Aren't the point kinases? Why would you be annotating with ligands then? I'm a bit confused
                plt.annotate(label, (nrg[i], dav[i]), fontsize=5)
            plt.yscale("log")  # cant do log for xscale bc negative values
            plt.title(f"Scores value comparison for {ligand}")
            plt.xlabel(ranked_df.columns[ranked_df.nrg_index])
            plt.ylabel(ranked_df.columns[ranked_df.davindex])
            plt.show()

        if ranks_plotting:
            nrg_rank_map = {name: i + 1 for i, name in enumerate(ranked_df.nrgnames)}  # mapping them as indexes
            dav_rank_map = {name: i + 1 for i, name in enumerate(ranked_df.davnames)}
            # new rank columns ฅ≽(•⩊ •マฅ≼)
            ranked_df["nrg_rank"] = ranked_df["Kinase"].map(nrg_rank_map)
            ranked_df["dav_rank"] = ranked_df["Kinase"].map(dav_rank_map)
            plt.figure(figsize=(8, 8), dpi=200)
            plt.scatter(ranked_df['nrg_rank'], ranked_df['dav_rank'])
            for i in range(len(ranked_df)):
                plt.annotate(ranked_df["Kinase"].iloc[i], (ranked_df["nrg_rank"].iloc[i],
                                                           ranked_df["dav_rank"].iloc[i]), fontsize=8, ha='center')
            plt.title(f"Ranks comparison for each kinase against {ligand}")
            plt.xlabel("nrg ranks")
            plt.ylabel("dav ranks")
            plt.tight_layout()
            plt.show()

        if roc_curving:
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(f"ROC curving for all ligands combined")
            plt.legend(loc="lower right", fontsize=8)
            plt.ylim([0.0, 1])
            plt.xlim([0.0, 1])

        # not tested this so far
        if contingency_table:
            davis_col = ranked_df.iloc[:, ranked_df.davindex]
            nrg_col = ranked_df.iloc[:, ranked_df.nrg_index]

            mask_davis = (davis_col < 20000).values
            mask_nrgrank = (nrg_col < 100000000).values

            contingency = pd.crosstab(mask_davis, mask_nrgrank,
                                      rownames=["Davis"], colnames=["nrgrank"])
            contingency.index = contingency.index.map({True: "binding", False: "no-binding"})
            contingency.columns = contingency.columns.map({True: "binding", False: "no-binding"})

            print("\nContingency Table:")
            print(contingency)
    return ranked_df


# -------------------------------------------------------------------------
if __name__ == '__main__':
    id_list, file_list = __list_files_cleft(
        input_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_pdb_kinases")

    #obabel_translate(input_dir=r"/lustre09/project/6049264/sbernarr/MAKI/DUAL_kinases_b4trans",output_path=r"/lustre09/project/6049264/sbernarr/MAKI/DUAL_kinases_translated")

    davis_df = pd.read_excel(r"/lustre09/project/6049264/sbernarr/MAKI/data/filtered2_davis_data.xlsx")

    nrgrank_df, ligands = creating_nrgrank_df(input_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_nrgrank",
                                              verbose=True,
                                              output_path=r"/lustre09/project/6049264/sbernarr/MAKI/data/nrgrank_dataframe.csv")
    #nrgrank_df = pd.read_csv(r"/lustre09/project/6049264/sbernarr/MAKI/data/nrgrank_dataframe.csv")

    overlapping_kinases = find_kins_overlapping_dfs(empiric_df=davis_df, experimental_df=nrgrank_df, verbose=False)



    #getcleft(input_dir=r"/lustre09/project/6049264/sbernarr/MAKI/DUAL_kinases_b4trans",
    #         output_clefts_proteome=r"/lustre09/project/6049264/sbernarr/MAKI/DUAL_kinase_clefts")

    # nrgrank(ligands_dir=r"/lustre09/project/6049264/sbernarr/MAKI/DUAL_ligands",
    #         root_cleft_dir=r"/lustre09/project/6049264/sbernarr/MAKI/DUAL_kinase_clefts",
    #         output_dir=r"/lustre09/project/6049264/sbernarr/MAKI/DUAL_nrgrank",
    #         kinase_mol_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_pdb_to_mol2_kinases")

    merging_dfs(nrgrank_df, davis_df, output_path=r"/lustre09/project/6049264/sbernarr/MAKI/",overlapping_kinases=overlapping_kinases)

    #plotting(ligands,ranked_df,values_plotting=True, ranks_plotting=True,roc_curving=True, contingency_table=True)



# ------------------------------------------------------------------------- future argumentalisation
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


# def download_alphafold_structures(output_dir=r"/lustre09/project/6049264/sbernarr/MAKI/data/DUAL_250_notKinases",
#                                    taxon_id=9606, # homosapiens
#                                    min_len=200,
#                                    max_len=300,
#                                    exclude_kinase=True,
#                                    limit=250):
#     os.makedirs(output_dir, exist_ok=True)
#
#     query = f"(organism_id:{taxon_id}) AND (length:[{min_len} TO {max_len}])"
#     if exclude_kinase:
#         query += " NOT annotation:(type:keyword kinase)"
#
#     url = "https://rest.uniprot.org/uniprotkb/search"
#     params = {
#         "query": query,
#         "format": "json",
#         "fields": "accession,proteinName,length",
#         "size": limit * 2 # ??
#     }
#
#     r = requests.get(url, params=params)
#     r.raise_for_status()
#     results = r.json()
#
#     count = 0
#     for entry in results.get("results", []):
#         uniprot_id = entry["primaryAccession"]
#         pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
#         pdb_path = os.path.join(output_dir, f"{uniprot_id}.pdb")
#
#         try:
#             pdb_response = requests.get(pdb_url)
#             if pdb_response.status_code == 200 and "ATOM" in pdb_response.text:
#                 with open(pdb_path, "w") as f:
#                     f.write(pdb_response.text)
#                 count += 1
#                 print(f"≽(•⩊ •マ≼  for {uniprot_id})")
#             else:
#                 print(f"not found ; skipping {uniprot_id}")
#         except Exception as e:
#             print(f"problem w {uniprot_id}: {e}")
#
#         if count >= limit:
#             break
#         time.sleep(0.2) # otherwise it fucks
#     print(f"files : {output_dir}")