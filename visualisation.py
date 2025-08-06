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
    return list(kinase_names)

def creating_nrgrank_df(input_dir):
    print("ฅ≽(•⩊ •マฅ≼")

def overlap_theo_exp(empiric_df, experimental_df):
    empiric_kinases = set(empiric_df['Kinase'])
    experimental_kinases = set(experimental_df['Kinase'])
    overlapping_kinases = empiric_kinases.intersection(experimental_kinases)
    return sorted(list(overlapping_kinases))

#nrgrank serves as experimental vs davis is empirical
def merging_dataframes(nrgrank_df, davis_df, output_path,overlapping_kinases):
    #todo; implement this logic for multiple ligands. creates plots / vis tool per ligand
    nrgrank_df = pd.read_csv(nrgrank_df)[["Kinase", "results_mastinib"]]
    nrgrank_df["Kinase"] = nrgrank_df["Kinase"].str.split("_").str[1] # name treatment, remove group
    nrgrank_df.rename(columns={'results_mastinib': 'res_mastinib_nrg'}, inplace=True)

    davis_df = pd.read_csv(davis_df)[["Kinase", "results_mastinib"]]
    davis_df.rename(columns={'results_mastinib': 'res_mastinib_davis'}, inplace=True)

    overlapping_kinases=overlapping_kinases
    nrgrank_df = nrgrank_df[nrgrank_df["Kinase"].isin(overlapping_kinases)]
    davis_df = davis_df[davis_df["Kinase"].isin(overlapping_kinases)]
    merged_df = pd.merge(nrgrank_df, davis_df, on="Kinase", how="inner")


    nrg_sorted = merged_df.sort_values(by='res_mastinib_nrg')
    dav_sorted = merged_df.sort_values(by='res_mastinib_davis')
    nrgnames = nrg_sorted['Kinase'].tolist()
    davnames = dav_sorted['Kinase'].tolist()

    ranked_df = merged_df.copy()
    ranked_df['nrgnames_sorted'] = pd.Series(nrgnames)
    ranked_df['davnames_sorted'] = pd.Series(davnames)
    ranked_df.to_csv(output_path, index=False)

    # # map them as ranks indexes
    # nrg_rank_map = {name: i + 1 for i, name in enumerate(nrgnames)}
    # dav_rank_map = {name: i + 1 for i, name in enumerate(davnames)}

    # # plot without ranking ≽(•⩊ •マ≼
    # nrg = merged_df.iloc[:, 1]
    # dav = merged_df.iloc[:, 2]
    # labels = merged_df.iloc[:, 0]
    # plt.figure(figsize=[5, 5])
    # plt.scatter(nrg, dav)
    # for i, label in enumerate(labels):
    #     plt.annotate(label, (nrg[i], dav[i]), fontsize=5)
    # plt.yscale("log")  # cant do log for xscale bc negative values
    # plt.title("mastinib")
    # plt.xlabel(merged_df.columns[1])
    # plt.ylabel(merged_df.columns[2])
    # plt.show()
    #
    # #plot with ranking ฅ≽(•⩊ •マฅ≼
    # ranked_df["nrg_rank"] = ranked_df["Kinase"].map(nrg_rank_map)
    # ranked_df["dav_rank"] = ranked_df["Kinase"].map(dav_rank_map)
    # plt.figure(figsize=(8, 5), dpi=200)
    # plt.scatter(ranked_df['nrg_rank'], ranked_df['dav_rank'])
    # for i in range(len(ranked_df)):
    #     plt.annotate(ranked_df["Kinase"].iloc[i], (ranked_df["nrg_rank"].iloc[i],
    #                                                ranked_df["dav_rank"].iloc[i]), fontsize=8, ha='center')
    # plt.xlabel("nrg ranks")
    # plt.ylabel("dav ranks")
    # plt.tight_layout()
    # plt.show()

    #punnet square creating (binairy value)

    return ranked_df

if __name__ == '__main__':
    #kinase_names_AF2 = extract_kinase_names(directory_path=r"C:\Users\sophi\Documents\NRG\KINASES\Kincore_AlphaFold2_ActiveHumanCatalyticKinases_v2\Kincore_AF2_HumanCatalyticKinases")

    overlapping_kinases= overlap_theo_exp(empiric_df=r"C:\Users\sophi\Documents\NRG\KINASES\filtered_sample_davis_data.xlsx", experimental_df=r"C:\Users\sophi\Documents\NRG\dataframe_cf_kinases.csv")

    merging_dataframes(davis_df=r"C:\Users\sophi\Documents\NRG\KINASES\filtered_sample_davis_data.xlsx", nrgrank_df=r"C:\Users\sophi\Documents\NRG\dataframe_cf_kinases.csv",output_path=r"C:\Users\sophi\Documents\NRG\KINASES\merged_df.csv")