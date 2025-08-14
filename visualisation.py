import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.metrics import roc_curve, auc
import numpy as np
from collections import defaultdict
import openpyxl
import csv

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
                else:
                    continue
    nrgrank_df = pd.DataFrame.from_dict(data, orient='index')
    nrgrank_df.index.name = 'Kinase'
    nrgrank_df = nrgrank_df[ligands]
    nrgrank_df.to_csv(output_path)
    return nrgrank_df, ligands

def overlap_theo_exp(empiric_df, experimental_df, verbose):
    empiric_kinases = set(empiric_df["Kinase"])

    experimental_df = experimental_df.reset_index()
    experimental_kinases = set(experimental_df["Kinase"])

    overlapping_kinases = empiric_kinases.intersection(experimental_kinases)
    verbose and print("≽(•⩊ •マ≼ Overlapping kinases found : ", overlapping_kinases)
    return sorted(list(overlapping_kinases))

def merging_df_and_plotting(nrgrank_df, davis_df, output_path, overlapping_kinases):
    ligands = nrgrank_df.columns.tolist()[1:]
    merged_df = pd.DataFrame({"Kinase": overlapping_kinases})
    #davis_df = davis_df.fillna(20000)

    nrg_index = 1
    dav_index = 2

    # ini figure of ROC curves
    roc_data = []
    colors = cm.tab20(np.linspace(0, 1, len(ligands)))
    plt.figure(figsize=(8, 8), dpi=200)
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')

    for ligand in ligands:
        #avoiding error
        nrgrank_sub = nrgrank_df[["Kinase", ligand]].copy()
        davis_sub = davis_df[["Kinase", ligand]].copy()

        nrgrank_sub.rename(columns={ligand: f"{ligand}_nrg"}, inplace=False)
        davis_sub.rename(columns={ligand: f"{ligand}_dav"}, inplace=False)

        #todo: outer here
        merged_df = pd.merge(merged_df, nrgrank_sub, on="Kinase", suffixes=("_x","_dav"),  how="inner")
        merged_df = pd.merge(merged_df, davis_sub, on="Kinase", suffixes=("_y","_nrg"), how="inner")
        merged_df.to_csv(output_path, index=False)

        # rank organising per ligand
        nrg_sorted = merged_df.sort_values(by=merged_df.columns[nrg_index])
        nrgnames = nrg_sorted.iloc[:, 0].tolist()
        dav_sorted = merged_df.sort_values(by=merged_df.columns[dav_index])
        davnames = dav_sorted.iloc[:, 0].tolist()

        ranked_df = merged_df.copy()
        ranked_df['nrgnames_sorted'] = pd.Series(nrgnames)
        ranked_df['davnames_sorted'] = pd.Series(davnames)

        nrg = ranked_df.iloc[:, nrg_index]
        dav = ranked_df.iloc[:, dav_index]
        #skipping columns to maintain good indexing
        nrg_index += 2
        dav_index += 2

        # ROC curve implementation
        davis_true = (dav < 20000).astype(int)
        nrg_inverted_scores = -nrg
        fpr, tpr, _ = roc_curve(davis_true, nrg_inverted_scores)
        roc_auc = auc(fpr, tpr)
        roc_data.append((ligand, fpr, tpr, roc_auc, colors[len(roc_data)]))
        plt.plot(fpr, tpr, color=colors[len(roc_data) - 1],
                 label=f"{ligand} (AUC = {roc_auc:.2f})")
    #------------------------------------------------------------------------- PLOTTING
        # plt.figure(figsize=(8, 8), dpi=200)
        # plt.scatter(nrg, dav)
        # for i, label in enumerate(labels):
        #     plt.annotate(label, (nrg[i], dav[i]), fontsize=5)
        # plt.yscale("log")  # cant do log for xscale bc negative values
        # plt.title(f"Scores value comparison for {ligand}")
        # plt.xlabel(ranked_df.columns[nrg_index])
        # plt.ylabel(ranked_df.columns[dav_index])
        # plt.show()

    #     # plot the ranks ≽(•⩊ •マ≼
    #     nrg_rank_map = {name: i + 1 for i, name in enumerate(nrgnames)} # mapping them as indexes
    #     dav_rank_map = {name: i + 1 for i, name in enumerate(davnames)}
    #     # new rank columns ฅ≽(•⩊ •マฅ≼)
    #     ranked_df["nrg_rank"] = ranked_df["Kinase"].map(nrg_rank_map)
    #     ranked_df["dav_rank"] = ranked_df["Kinase"].map(dav_rank_map)
    #     plt.figure(figsize=(8, 8), dpi=200)
    #     plt.scatter(ranked_df['nrg_rank'], ranked_df['dav_rank'])
    #     for i in range(len(ranked_df)):
    #         plt.annotate(ranked_df["Kinase"].iloc[i], (ranked_df["nrg_rank"].iloc[i],
    #                                                    ranked_df["dav_rank"].iloc[i]), fontsize=8, ha='center')
    #     plt.title(f"Ranks comparison for each kinase against {ligand}")
    #     plt.xlabel("nrg ranks")
    #     plt.ylabel("dav ranks")
    #     plt.tight_layout()
    #     plt.show()

    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC curving for all ligands combined")
    plt.legend(loc="lower right", fontsize=8)
    plt.ylim([0.0, 1])
    plt.xlim([0.0, 1])
    plt.show()

    return ranked_df

if __name__ == '__main__':
    # davis_df = pd.read_excel(r"/lustre09/project/6049264/sbernarr/MAKI/data/filtered2_davis_data.xlsx")
    # nrgrank_df,ligands = creating_nrgrank_df(input_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_nrgrank",verbose=True,output_path=r"/lustre09/project/6049264/sbernarr/MAKI/data/nrgrank_dataframe.csv")
    # nrgrank_df = pd.read_csv(r"/lustre09/project/6049264/sbernarr/MAKI/data/nrgrank_dataframe.csv")
    # merged_df = merging_dataframes(nrgrank_df=nrgrank_df,davis_df=davis_df,output_path=r"/lustre09/project/6049264/sbernarr/MAKI/data/merged_df.csv",overlapping_kinases=overlapping_kinases)

    davis_df = pd.read_excel(r"C:\Users\sophi\Documents\NRG\KINASES\filtered2_davis_data.xlsx")
    nrgrank_df\
        = pd.read_csv(r"C:\Users\sophi\Documents\NRG\KINASES\nrgrank_dataframe.csv")
    overlapping_kinases = overlap_theo_exp(empiric_df=davis_df, experimental_df=nrgrank_df, verbose=False)
    merged_df = merging_df_and_plotting(nrgrank_df=nrgrank_df, davis_df=davis_df,
                                   output_path=r"C:\Users\sophi\Documents\NRG\KINASES\merged_df.csv",
                                   overlapping_kinases=overlapping_kinases)

