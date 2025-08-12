import os
import pandas as pd
import matplotlib.pyplot as plt
import csv
import openpyxl
from collections import defaultdict


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


def merging_dataframes(nrgrank_df, davis_df, output_path,
                       overlapping_kinases):  # resetting index so Kinase goes back to column
    ligands = nrgrank_df.columns.tolist()[1:]
    merged_df = pd.DataFrame({"Kinase": overlapping_kinases})
    davis_df = davis_df.fillna(15000)

    for ligand in ligands:  # ≽(•⩊ •マ≼ ≽(•⩊ •マ≼
        nrgrank_sub = nrgrank_df[["Kinase", ligand]].copy()
        davis_sub = davis_df[["Kinase", ligand]].copy()
        nrgrank_sub.rename(columns={ligand: f"{ligand}_nrg"}, inplace=True)
        davis_sub.rename(columns={ligand: f"{ligand}_dav"}, inplace=True)
        merged_df = pd.merge(merged_df, nrgrank_sub, on="Kinase", how="outer")
        merged_df = pd.merge(merged_df, davis_sub, on="Kinase", how="outer")
    merged_df.to_csv(output_path, index=False)

    # rank organising per ligand
    nrg_sorted = merged_df.sort_values(by=merged_df.columns[1])
    nrgnames = nrg_sorted.iloc[:, 0].tolist()
    dav_sorted = merged_df.sort_values(by=merged_df.columns[2])
    davnames = dav_sorted.iloc[:, 0].tolist()
    ranked_df = merged_df.copy()
    ranked_df['nrgnames_sorted'] = pd.Series(nrgnames)
    ranked_df['davnames_sorted'] = pd.Series(davnames)

    # initial plot with all of the ligands in the same plot, color per ligand
    nrg = merged_df.iloc[:, 1]
    dav = merged_df.iloc[:, 2]
    labels = merged_df.iloc[:, 0]
    plt.figure(figsize=(8, 5), dpi=200)
    plt.scatter(nrg, dav)
    for i, label in enumerate(labels):
        plt.annotate(label, (nrg[i], dav[i]), fontsize=5)
    plt.yscale("log")  # cant do log for xscale bc negative values
    plt.title("maaaaa")
    plt.xlabel(merged_df.columns[1])
    plt.ylabel(merged_df.columns[2])
    plt.show()

    # plot the ranks ≽(•⩊ •マ≼
    # map them as ranks indexes
    nrg_rank_map = {name: i + 1 for i, name in enumerate(nrgnames)}
    dav_rank_map = {name: i + 1 for i, name in enumerate(davnames)}
    # new rank columsn ฅ≽(•⩊ •マฅ≼)
    ranked_df["nrg_rank"] = ranked_df["Kinase"].map(nrg_rank_map)
    ranked_df["dav_rank"] = ranked_df["Kinase"].map(dav_rank_map)
    plt.figure(figsize=(8, 5), dpi=200)
    plt.scatter(ranked_df['nrg_rank'], ranked_df['dav_rank'])
    for i in range(len(ranked_df)):
        plt.annotate(ranked_df["Kinase"].iloc[i], (ranked_df["nrg_rank"].iloc[i],
                                                   ranked_df["dav_rank"].iloc[i]), fontsize=8, ha='center')
    plt.xlabel("nrg ranks")
    plt.ylabel("dav ranks")
    plt.tight_layout()
    plt.show()
    return ranked_df


if __name__ == '__main__':
    # davis_df = pd.read_excel(r"/lustre09/project/6049264/sbernarr/MAKI/data/filtered2_davis_data.xlsx")
    # nrgrank_df,ligands = creating_nrgrank_df(input_dir=r"/lustre09/project/6049264/sbernarr/MAKI/AF2_nrgrank",verbose=True,output_path=r"/lustre09/project/6049264/sbernarr/MAKI/data/nrgrank_dataframe.csv")
    # nrgrank_df = pd.read_csv(r"/lustre09/project/6049264/sbernarr/MAKI/data/nrgrank_dataframe.csv")
    # merged_df = merging_dataframes(nrgrank_df=nrgrank_df,davis_df=davis_df,output_path=r"/lustre09/project/6049264/sbernarr/MAKI/data/merged_df.csv",overlapping_kinases=overlapping_kinases)

    davis_df = pd.read_excel(r"/lustre09/project/6049264/sbernarr/MAKI/data/filtered2_davis_data.xlsx")
    nrgrank_df = pd.read_csv(r"/lustre09/project/6049264/sbernarr/MAKI/data/nrgrank_dataframe.csv")
    overlapping_kinases = overlap_theo_exp(empiric_df=davis_df, experimental_df=nrgrank_df, verbose=False)
    merged_df = merging_dataframes(nrgrank_df=nrgrank_df, davis_df=davis_df,
                                   output_path=r"/lustre09/project/6049264/sbernarr/MAKI/data/merged_df.csv",
                                   overlapping_kinases=overlapping_kinases)

