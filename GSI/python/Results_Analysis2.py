from pathlib import Path
import pandas as pd

root_dir = Path.cwd().parent
data_dir = root_dir / "data"
sol_dir = root_dir / "solution"
digestion_dir = data_dir / "digestion"

upper_edges = pd.read_csv(digestion_dir / "digestion_yeast+ups1_result.csv")[["accession", "protein_id", "peptide_id", "Prob"]]
lower_edges = pd.read_csv(sol_dir / "lower_edges.csv")
N_prot = len(upper_edges.groupby("accession"))

ref = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";")
sol = pd.read_csv(sol_dir / "results_yeast_10fmol100.csv")
ref = pd.merge(ref, upper_edges[["accession","protein_id"]].drop_duplicates(), left_on = "Accession", right_on="accession", how = 'left').drop('accession', axis = 1)

protein_to_spectra = pd.merge(upper_edges, lower_edges, left_on = "peptide_id", right_on = "Peptide", how = 'left').drop("Peptide", axis = 1)
protein_to_spectra["protein_prediction"] = protein_to_spectra["protein_id"].isin(sol['id'])
protein_to_spectra["protein_truth"] = protein_to_spectra["protein_id"].isin(ref['protein_id'])
protein_to_spectra["has_spectrum"] = protein_to_spectra["protein_id"].isin(protein_to_spectra.dropna()["protein_id"])

TP_prot = protein_to_spectra.loc[(protein_to_spectra["protein_prediction"] == True) & (protein_to_spectra["protein_truth"] == True)]
FP_prot = protein_to_spectra.loc[(protein_to_spectra["protein_prediction"] == True) & (protein_to_spectra["protein_truth"] == False)]
TN_prot = protein_to_spectra.loc[(protein_to_spectra["protein_prediction"] == False) & (protein_to_spectra["protein_truth"] == False)]
FN_prot = protein_to_spectra.loc[(protein_to_spectra["protein_prediction"] == False) & (protein_to_spectra["protein_truth"] == True)]

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

protein_to_spectra["prediction_category"] = protein_to_spectra.apply(prediction_category, axis = 1)

def print_stats(df):
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