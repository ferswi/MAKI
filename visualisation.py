import os
import pandas as pd
import matplotlib.pyplot as plt

def extract_kinase_names(directory_path):
    kinase_names = set()
    for filename in os.listdir(directory_path):
        if "_" in filename:
            parts = filename.split("_")
            if len(parts) > 2:
                kinase_name = parts[1]  # between first and second underscore
                kinase_names.add(kinase_name)
    print(list(kinase_names))
    print(len(kinase_names))
    return list(kinase_names)

def overlap_theo_exp(empiric_df, experimental_df):



def merging_dataframes(nrgrank_df, davis_df, output_path):
    nrgrank_df = pd.read_csv(nrgrank_df)[["Kinase","results_mastinib"]]
    nrgrank_df["Kinase"]=nrgrank_df["Kinase"].str.split("_").str[1]
    nrgrank_df.rename(columns={'results_mastinib': 'res_mastinib_nrg'}, inplace=True)
    davis_df = pd.read_csv(davis_df)[["Kinase","results_mastinib"]]
    davis_df.rename(columns={'results_mastinib': 'res_mastinib_davis'}, inplace=True)
    merged_df = pd.merge(nrgrank_df, davis_df,on="Kinase")

    nrg_sorted = merged_df.sort_values(by=merged_df.columns[1])
    nrgnames = nrg_sorted.iloc[:, 0].tolist()
    dav_sorted = merged_df.sort_values(by=merged_df.columns[2])
    davnames = dav_sorted.iloc[:, 0].tolist()

    ranked_df = merged_df.copy()
    ranked_df['nrgnames_sorted'] = pd.Series(nrgnames)
    ranked_df['davnames_sorted'] = pd.Series(davnames)
    # plot the ranks ≽(•⩊ •マ≼
    # map them as ranks indexes
    nrg_rank_map = {name: i + 1 for i, name in enumerate(nrgnames)}
    dav_rank_map = {name: i + 1 for i, name in enumerate(davnames)}

    # plot without ranking
    nrg = merged_df.iloc[:, 1]
    dav = merged_df.iloc[:, 2]
    labels = merged_df.iloc[:, 0]
    plt.figure(figsize=[5, 5])
    plt.scatter(nrg, dav)
    for i, label in enumerate(labels):
        plt.annotate(label, (nrg[i], dav[i]), fontsize=5)
    plt.yscale("log")  # cant do log for xscale bc negative values
    plt.title("mastinib")
    plt.xlabel(merged_df.columns[1])
    plt.ylabel(merged_df.columns[2])
    plt.show()

    #plot with ranking
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

if __name__ == '__main__':
    kinase_names_AF2 = extract_kinase_names(directory_path=r"C:\Users\sophi\Documents\NRG\KINASES\Kincore_AlphaFold2_ActiveHumanCatalyticKinases_v2\Kincore_AF2_HumanCatalyticKinases")
    overlap_theo_exp(empiric_df=r"C:\Users\sophi\Documents\NRG\KINASES\csv_davis_to_compare_sample_to.csv", experimental_df=)