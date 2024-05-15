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

class Results_Analysis:

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
    
    def analyse_df(self, dataset_name = "", notes = "", thresh = 0, max_edges = 0) -> pd.DataFrame:
        true_proteins = self.protein_to_spectra.loc[self.protein_to_spectra["protein_truth"]]
        false_proteins = self.protein_to_spectra.loc[self.protein_to_spectra["protein_truth"] == False]
        df = self.protein_to_spectra[["accession","prediction_category"]].drop_duplicates()
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
        return pd.DataFrame({
            "Dataset": dataset_name,
            "Notes": notes,
            "Threshold": thresh,
            "Max Score Edges": max_edges,
            "Proteins number": self.N_prot,
            "Peptides number": len(self.lower_edges["Peptide"].drop_duplicates()),
            "Protein Edges number": len(self.upper_edges),
            "Proteins with Spectra number": len(self.protein_to_spectra.dropna()["protein_id"].drop_duplicates()),
            "Proteins Spectra Paths": len(self.protein_to_spectra[["protein_id","Spectrum"]].dropna()),
            "Target Proteins": len(true_proteins["protein_id"].drop_duplicates()),
            "Target Proteins associated Spectra": len(true_proteins.dropna()),
            "Target Proteins with Spectra": len(true_proteins.dropna()["protein_id"].drop_duplicates()),
            "Decoy Proteins": len(false_proteins["protein_id"].drop_duplicates()),
            "Decoy Proteins associated Spectra": len(false_proteins.dropna()),
            "Decoy Proteins with Spectra": len(false_proteins.dropna()["protein_id"].drop_duplicates()),
            "Score Edges Number": len(self.lower_edges),
            "Peptides with Spectra number": len(self.lower_edges["Peptide"].drop_duplicates()),
            "True Positives": TP,
            "False Positives": FP,
            "True Negatives": TN,
            "False Negatives": FN,
            "Accuracy": accuracy,
            "Specificity": specificity,
            "Sensitivity": sensitivity,
            "False Negatives Rate": FNR,
            "Positive Predictive Value": PPV,
            "Negative Predictive Value": NPV},
            index = [0])

class Model_Analyses:

    def __init__(self) -> None:
        self.analyses_df = pd.DataFrame(columns = [
            "Dataset",
            "Notes",
            "Threshold",
            "Max Score Edges",
            "Proteins",
            "Peptides"
            "Protein Edges",
            "Proteins with Spectra",
            "Proteins Spectra Paths",
            "Target Proteins",
            "Target Proteins associated Spectra",
            "Target Proteins with Spectra",
            "Decoy Proteins",
            "Decoy Proteins associated Spectra",
            "Decoy Proteins with Spectra",
            "Score Edges Number",
            "Peptides with Spectra number",
            "True Positives",
            "False Positives",
            "True Negatives",
            "False Negatives",
            "Accuracy",
            "Specificity",
            "Sensitivity",
            "False Negatives Rate",
            "Positive Predictive Value",
            "Negative Predictive Value"])
    
    def Add_Analysis(self, analyse: pd.DataFrame) -> None:
        self.analyses_df = self.analyses_df.append(analyse)
#%%
for threshold in [7500]:
    for n_edges in [4,10,0]:
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
lower_edges = pd.read_csv(sol_dir / f"lower_edges2000_10.csv")
N_prot = len(upper_edges.groupby("accession"))

ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
sol = pd.read_csv(sol_dir / f"results_yeast_10fmol2000_10.csv")
ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)

results = Results_Analysis(upper_edges, lower_edges, ref, sol, min_proba = 0)

print(f"==============================================================\nResults for threshold {threshold} and max edges {n_edges}")
# results.print_stats_proteins()
results.print_stats_scores()
results.print_stats_true_proteins()
results.print_stats_false_proteins()
results.print_stats_predictions()
print("==============================================================")
# %%
fig, axs = plt.subplots(2,1, sharex = True, figsize = (7,10))
thresh, n_edge, n_prot, n_scores = [], [], [], []
for (threshold, n_edges) in [(7,2),(7,3),(7,4),(7,10),(12,2),(12,3),(12,4),(17,2),(17,3),(17,4),(17,10),(27,2),(27,3),(27,4),(27,10),(37,2),(37,3),(37,4),
                             (100,4),(100,5),
                             (200,4),
                             (500,4),(500,10),
                             (1000,4),
                             (2000,4),
                             (3000,4),(3000,10),
                             (4000,10),
                             (5000,4),(5000,10),(5000,20),
                             (7500,4),(7500,10)]:
    upper_edges = pd.read_csv(digestion_dir / "digestion_yeast+ups1_result.csv")[["accession", "protein_id", "peptide_id", "Prob"]]
    lower_edges = pd.read_csv(sol_dir / f"lower_edges{threshold}_{n_edges}.csv")
    N_prot = len(upper_edges.groupby("accession"))

    ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
    sol = pd.read_csv(sol_dir / f"results_yeast_10fmol{threshold}_{n_edges}.csv")
    ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)
    results = Results_Analysis(upper_edges, lower_edges, ref, sol)
    thresh.append(threshold)
    n_edge.append(n_edges)
    true_proteins = results.protein_to_spectra.loc[results.protein_to_spectra["protein_truth"]]
    n_prot.append(len(true_proteins.dropna()["protein_id"].drop_duplicates()))
    n_scores.append(len(results.lower_edges))
colormap = {2:'k',3:'k',4:'b',5:'b',10:'r',20:'g'}
color = [colormap[i] for i in n_edge]
axs[0].scatter(thresh, n_prot, c = color, marker = '+')
axs[0].set_xscale('log')
axs[1].scatter(thresh, n_scores, c = color, marker = '+')
axs[1].set_xscale('log')
plt.show()
fig,axs = plt.subplots(1,2, sharey = True, figsize = (14,5))
axs[0].scatter(n_prot, n_scores, c = thresh, marker = '+')
axs[1].scatter(n_prot, n_scores, c = color, marker = '+')
plt.show()
# %%
threshold = 2000
n_edges = 10
fig, axs = plt.subplots(2,1, sharex = True, figsize = (7,10))
upper_edges = pd.read_csv(digestion_dir / "digestion_yeast+ups1_result.csv")[["accession", "protein_id", "peptide_id", "Prob"]]
lower_edges = pd.read_csv(sol_dir / f"lower_edges{threshold}_{n_edges}.csv")
N_prot = len(upper_edges.groupby("accession"))

ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
sol = pd.read_csv(sol_dir / f"results_yeast_10fmol{threshold}_{n_edges}.csv")
ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)
results = Results_Analysis(upper_edges, lower_edges, ref, sol)
thresh.append(threshold)
n_edge.append(n_edges)
true_proteins = results.protein_to_spectra.loc[results.protein_to_spectra["protein_truth"]]
n_prot.append(len(true_proteins.dropna()["protein_id"].drop_duplicates()))
n_scores.append(len(results.lower_edges))
colormap = {2:'k',3:'k',4:'b',5:'b',10:'r',20:'g'}
color = [colormap[i] for i in n_edge]
axs[0].scatter(thresh, n_prot, c = color, marker = '+')
axs[0].set_xscale('log')
axs[1].scatter(thresh, n_scores, c = color, marker = '+')
axs[1].set_xscale('log')
plt.show()
fig,axs = plt.subplots(1,2, sharey = True, figsize = (14,5))
axs[0].scatter(n_prot, n_scores, c = thresh, marker = '+')
axs[1].scatter(n_prot, n_scores, c = color, marker = '+')
plt.show()
# %%
Analyses = Model_Analyses()
for threshold in [7500]:
    for n_edges in [4,10,0]:
        upper_edges = pd.read_csv(digestion_dir / "digestion_yeast+ups1_result.csv")[["accession", "protein_id", "peptide_id", "Prob"]]
        lower_edges = pd.read_csv(sol_dir / f"lower_edges{threshold}_{n_edges}.csv")
        N_prot = len(upper_edges.groupby("accession"))

        ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
        sol = pd.read_csv(sol_dir / f"results_yeast_10fmol{threshold}_{n_edges}.csv")
        ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)
        results = Results_Analysis(upper_edges, lower_edges, ref, sol)
        Analyses.Add_Analysis(results.analyse_df(dataset_name = "UPS+Yeast", thresh = threshold, max_edges = n_edges))
# %%
