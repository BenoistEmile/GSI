#%%
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

root_dir = Path.cwd().parent
data_dir = root_dir / "data"
sol_dir = root_dir / "solution"
digestion_dir = data_dir / "digestion"

def prediction_category(row):
    if row["protein_prediction"] == True:
        if row["protein_truth"] == True:
            return "TP"
        else:
            return "FP"
    else:
        if row["protein_truth"] == True:
            return "FN"
        else:
            return "TN"

class Results_Analysis():

    def __init__(self, upper_edges: pd.DataFrame, lower_edges: pd.DataFrame, ref: pd.DataFrame, sol: pd.DataFrame, min_proba: float = 0) -> None:
        self.upper_edges = upper_edges.loc[upper_edges["Prob"] >= min_proba]
        self.lower_edges = lower_edges
        self.ref = ref
        self.sol = sol
        self.N_prot = len(self.upper_edges.groupby("accession"))
        self.protein_to_spectra = pd.merge(self.upper_edges, self.lower_edges, left_on = "peptide_id", right_on = "Peptide", how = 'left').drop("Peptide", axis = 1)
        self.protein_to_spectra["protein_prediction"] = self.protein_to_spectra["protein_id"].isin(self.sol['id'])
        self.protein_to_spectra["protein_truth"] = self.protein_to_spectra["protein_id"].isin(self.ref['protein_id'])
        self.protein_to_spectra["has_spectrum"] = self.protein_to_spectra["protein_id"].isin(self.protein_to_spectra.dropna()["protein_id"])
        self.protein_to_spectra["prediction_category"] = self.protein_to_spectra.apply(prediction_category, axis = 1)

    def print_stats_proteins(self) -> None:
        N_edges = len(self.upper_edges)
        median_prob = self.upper_edges["Prob"].median()
        max_prob = self.upper_edges["Prob"].max()
        min_prob = self.upper_edges["Prob"].min()
        mean_prob = self.upper_edges["Prob"].mean()
        std_prob = self.upper_edges["Prob"].std()
        self.upper_edges["Prob"].hist(bins = 50)
        N_prot_spectra = len(self.protein_to_spectra.dropna()["protein_id"].drop_duplicates())
        N_path_prot_spectra = len(self.protein_to_spectra[["protein_id","Spectrum"]].dropna())
        print(f"""Number of proteins : {N_prot}
Number_of edges : {N_edges}
Median probability : {round(median_prob, 2)}
Min probability : {round(min_prob, 2)}
Max probability : {round(max_prob, 2)}
Mean probability (std) : {round(mean_prob, 2)} ({round(std_prob, 2)})
Number of proteins with spectra : {N_prot_spectra}
Number of proteins-spectra paths : {N_path_prot_spectra}""")
        fig, axs = plt.subplots(1,2, sharey = True, figsize = (20,6))
        self.upper_edges.groupby("protein_id").count()["accession"].hist(bins = self.upper_edges.groupby("protein_id").count()["accession"].max(), ax = axs[0])
        self.protein_to_spectra[["protein_id","Spectrum"]].dropna().groupby("protein_id").count().hist(bins = self.protein_to_spectra[["protein_id","Spectrum"]].dropna().groupby("protein_id").count()["Spectrum"].max(), ax = axs[1])
        axs[0].set_title("Number of peptides per protein")
        axs[0].set_xlabel("Number of peptides")
        axs[0].set_ylabel("Number of proteins")
        axs[1].set_title("Number of spectra par protein")
        axs[1].set_xlabel("Number of spectra")
        plt.show()
    
    def print_stats_true_proteins(self) -> None:
        true_proteins = self.protein_to_spectra.loc[self.protein_to_spectra["protein_truth"]]
        N_proteins = len(true_proteins["protein_id"].drop_duplicates())
        N_spectra = len(true_proteins.dropna())
        N_prot_spectra = len(true_proteins.dropna()["protein_id"].drop_duplicates())
        print(f"""Number of true proteins : {N_proteins}
Number of associated spectra : {N_spectra}
Number of proteins with spectra : {N_prot_spectra}""")
        fig, axs = plt.subplots(1, 2, sharey = False, figsize = (20,6))
        true_proteins.dropna()["Score"].hist(bins = 100, ax = axs[0])
        true_proteins.dropna()[["protein_id","Score"]].groupby("protein_id").count().hist(bins = true_proteins.dropna()[["protein_id","Score"]].groupby("protein_id").count()["Score"].max(), ax = axs[1])
        axs[0].set_title("Scores distribution")
        axs[0].set_xlabel("Score")
        axs[0].set_ylabel("Number of spectra")
        axs[1].set_title("Number of spectra par protein")
        axs[1].set_xlabel("Number of spectra")
        axs[1].set_ylabel("Number of proteins")
        plt.show()

    def print_stats_false_proteins(self) -> None:
        false_proteins = self.protein_to_spectra.loc[self.protein_to_spectra["protein_truth"] == False]
        N_proteins = len(false_proteins["protein_id"].drop_duplicates())
        N_spectra = len(false_proteins.dropna())
        N_prot_spectra = len(false_proteins.dropna()["protein_id"].drop_duplicates())
        print(f"""Number of false proteins : {N_proteins}
Number of associated spectra : {N_spectra}
Number of proteins with spectra : {N_prot_spectra}""")
        fig, axs = plt.subplots(1, 2, sharey = False, figsize = (20,6))
        false_proteins.dropna()["Score"].hist(bins = 100, ax = axs[0])
        false_proteins.dropna()[["protein_id","Score"]].groupby("protein_id").count().hist(bins = false_proteins.dropna()[["protein_id","Score"]].groupby("protein_id").count()["Score"].max(), ax = axs[1])
        axs[0].set_title("Scores distribution")
        axs[0].set_xlabel("Score")
        axs[0].set_ylabel("Number of spectra")
        axs[1].set_title("Number of spectra par protein")
        axs[1].set_xlabel("Number of spectra")
        axs[1].set_ylabel("Number of proteins")
        plt.show()

    def print_stats_scores(self) -> None:
        N_edges = len(self.lower_edges)
        N_peptides = len(self.lower_edges["Peptide"].drop_duplicates())
        median_score = self.lower_edges["Score"].median()
        min_score = self.lower_edges["Score"].min()
        max_score = self.lower_edges["Score"].max()
        mean_score = self.lower_edges["Score"].mean()
        std_score = self.lower_edges["Score"].std()
        print(f"""Number of edges : {N_edges}
Number of peptides with spectra : {N_peptides}
Median score : {round(median_score,2)}
Min score : {round(min_score,2)}
Max score : {round(max_score,2)}
Mean score (std) : {round(mean_score,2)} ({round(std_score,2)})""")
        fig, axs = plt.subplots(1,2, sharey = False, figsize = (20,6))
        self.lower_edges["Score"].hist(bins = 50, ax = axs[0])
        self.lower_edges.groupby("Peptide").count()["Spectrum"].hist(bins = self.lower_edges.groupby("Peptide").count()["Spectrum"].max(), ax = axs[1])
        axs[0].set_title("Scores distribution")
        axs[0].set_xlabel("Score")
        axs[0].set_ylabel("Number of spectra")
        axs[1].set_title("Number of spectra per peptide")
        axs[1].set_xlabel("Number of spectra")
        axs[1].set_ylabel("Number of peptide")
        plt.show()

    def print_stats_predictions(self):
        df = self.protein_to_spectra[["accession","prediction_category"]].drop_duplicates()
        N_prot = len(df)
        TP = len(df[df["prediction_category"] == "TP"])
        FP = len(df[df["prediction_category"] == "FP"])
        TN = len(df[df["prediction_category"] == "TN"])
        FN = len(df[df["prediction_category"] == "FN"])
        accuracy = (TP + TN) / N_prot
        specificity = TN/ (TN + FP)
        sensitivity = TP / (TP + FN)
        FNR = 1 - sensitivity
        PPV = TP / (TP + FP)
        NPV = TN / (TN + FN)
        print(f"""True positives : {TP}
True negatives : {TN}
False positives : {FP}
False negatives : {FN}
Accuracy : {round(accuracy, 3)}
Sensitivity : {round(sensitivity, 3)}
Specificity : {round(specificity, 3)}
FNR : {round(FNR, 3)}
PPV : {round(PPV, 3)}
NPV : {round(NPV, 3)}""")
#%%
for threshold in [1000,3000,5000]:
    for n_edges in [4]:
        upper_edges = pd.read_csv(digestion_dir / "digestion_yeast+ups1_result.csv")[["accession", "protein_id", "peptide_id", "Prob"]]
        lower_edges = pd.read_csv(sol_dir / f"lower_edges{threshold}_{n_edges}.csv")
        N_prot = len(upper_edges.groupby("accession"))

        ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
        sol = pd.read_csv(sol_dir / f"results_yeast_10fmol{threshold}_{n_edges}.csv")
        ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)
        results = Results_Analysis(upper_edges, lower_edges, ref, sol)

        print(f"==============================================================\nResults for threshold {threshold} and max edges {n_edges}")
        # results.print_stats_proteins()
        results.print_stats_scores()
        # results.print_stats_true_proteins()
        # results.print_stats_false_proteins()
        # results.print_stats_predictions()
        print("==============================================================")
# %%
upper_edges = pd.read_csv(digestion_dir / "digestion_yeast+ups1_result.csv")[["accession", "protein_id", "peptide_id", "Prob"]]
lower_edges = pd.read_csv(sol_dir / f"lower_edges100_5.csv")
N_prot = len(upper_edges.groupby("accession"))

ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
sol = pd.read_csv(sol_dir / f"results_yeast_10fmol100_5.csv")
ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)

results = Results_Analysis(upper_edges, lower_edges, ref, sol, min_proba = 0)

print(f"==============================================================\nResults for threshold {threshold} and max edges {n_edges}")
results.print_stats_proteins()
results.print_stats_scores()
results.print_stats_true_proteins()
results.print_stats_false_proteins()
results.print_stats_predictions()
print("==============================================================")