from pathlib import Path
import pandas as pd

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

    def __init__(self, upper_edges: pd.DataFrame, lower_edges: pd.DataFrame, ref: pd.DataFrame, sol: pd.DataFrame) -> None:
        self.upper_edges = upper_edges
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
        N_prot_spectra = len(self.protein_to_spectra.dropna()["protein_id"].drop_duplicates())
        print(f"""Number of proteins : {N_prot}
Number_of edges : {N_edges}
Number of proteins with spectra : {N_prot_spectra}""")
        self.protein_to_spectra.dropna().groupby("protein_id").count().groupby("accession").size().hist(bins = self.protein_to_spectra.dropna().groupby("protein_id").count().groupby("accession").size().max()//100)


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
        self.lower_edges.groupby("Peptide").count()["Spectrum"].hist(bins = lower_edges.groupby("Peptide").count()["Spectrum"].max())

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
        
upper_edges = pd.read_csv(digestion_dir / "digestion_yeast+ups1_result.csv")[["accession", "protein_id", "peptide_id", "Prob"]]
lower_edges = pd.read_csv(sol_dir / "lower_edges12_4.csv")
N_prot = len(upper_edges.groupby("accession"))

ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
sol = pd.read_csv(sol_dir / "results_yeast_10fmol12_4.csv")
ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)

results = Results_Analysis(upper_edges, lower_edges, ref, sol)

results.print_stats_scores()
results.print_stats_predictions()