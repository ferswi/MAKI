import os
import pandas as pd
import matplotlib.pyplot as plt
import csv


def creating_nrgrank_df(input_dir):
    print("ฅ≽(•⩊ •マฅ≼")

def overlap_theo_exp(empiric_df, experimental_df):
    empiric_kinases = set(empiric_df["Kinase"])
    experimental_kinases = set(experimental_df["Kinase"])
    overlapping_kinases = empiric_kinases.intersection(experimental_kinases)
    return sorted(list(overlapping_kinases))

def merging_dataframes(nrgrank_df, davis_df, output_path,overlapping_kinases):
    #todo implement multi ligand treatment and prduction  w {ligand instead of mastinib } + loop
    nrgrank_df = pd.read_csv(nrgrank_df)[["Kinase","results_mastinib"]]
    nrgrank_df["Kinase"]=nrgrank_df["Kinase"].str.split("_").str[1]
    nrgrank_df.rename(columns={'results_mastinib': 'res_mastinib_nrg'}, inplace=True)
    davis_df = pd.read_csv(davis_df)[["Kinase","results_mastinib"]]
    davis_df.rename(columns={'results_mastinib': 'res_mastinib_davis'}, inplace=True)
    merged_df = pd.merge(nrgrank_df, davis_df,on="Kinase")

    # rank organising
    nrg_sorted = merged_df.sort_values(by=merged_df.columns[1])
    nrgnames = nrg_sorted.iloc[:, 0].tolist()
    dav_sorted = merged_df.sort_values(by=merged_df.columns[2])
    davnames = dav_sorted.iloc[:, 0].tolist()
    ranked_df = merged_df.copy()
    ranked_df['nrgnames_sorted'] = pd.Series(nrgnames)
    ranked_df['davnames_sorted'] = pd.Series(davnames)

    # initial plot
    nrg = merged_df.iloc[:, 1]
    dav = merged_df.iloc[:, 2]
    labels = merged_df.iloc[:, 0]
    plt.figure(figsize=[5,5])
    plt.scatter(nrg, dav)
    for i, label in enumerate(labels):
        plt.annotate(label, (nrg[i], dav[i]),fontsize=5)
    plt.yscale("log")# cant do log for xscale bc negative values
    plt.title("mastinib")
    plt.xlabel(merged_df.columns[1])
    plt.ylabel(merged_df.columns[2])
    plt.show()

    #plot the ranks ≽(•⩊ •マ≼
    #map them as ranks indexes
    nrg_rank_map = {name:i+1 for i, name in enumerate(nrgnames)}
    dav_rank_map = {name:i+1 for i, name in enumerate(davnames)}
    # new rank columsn ≽(•⩊ •マ≼ ≽(•⩊ •マ≼
    ranked_df["nrg_rank"] = ranked_df["Kinase"].map(nrg_rank_map)
    ranked_df["dav_rank"] = ranked_df["Kinase"].map(dav_rank_map)

    plt.figure(figsize=(8, 5), dpi=200)
    plt.scatter(ranked_df['nrg_rank'],ranked_df['dav_rank'])
    for i in range(len(ranked_df)):
        plt.annotate(ranked_df["Kinase"].iloc[i], (ranked_df["nrg_rank"].iloc[i],
                                                    ranked_df["dav_rank"].iloc[i]), fontsize = 8, ha = 'center')
    plt.xlabel("nrg ranks")
    plt.ylabel("dav ranks")
    plt.tight_layout()
    plt.show()

    #plot punnet squares of blanks ฅ≽(•⩊ •マฅ≼

    return ranked_df

if __name__ == '__main__':
    davis_df = pd.read_csv(r"C:\Users\sophi\Documents\NRG\KINASES\csv_davis_to_compare_sample_to.csv")
    nrgrank_df = pd.read_csv(r"C:\Users\sophi\Documents\NRG\KINASES\dataframe_cf_kinases.csv")

    overlapping_kinases = overlap_theo_exp(empiric_df=davis_df, experimental_df=nrgrank_df)

    merged_df = merging_dataframes(
        nrgrank_df=r"C:\Users\sophi\Documents\NRG\KINASES\dataframe_cf_kinases.csv",
        davis_df=r"C:\Users\sophi\Documents\NRG\KINASES\csv_davis_to_compare_sample_to.csv",
        output_path=r"C:\Users\sophi\Documents\NRG\KINASES\dataframe_cf_kinases_merged.csv",
        overlapping_kinases=overlapping_kinases
    )
