from pathlib import Path
import pandas as pd

gsi_path = Path.cwd()

proteins = pd.read_csv(gsi_path / 'data' / 'digestion' / 'digestion_file.csv')[["accession", "protein_id"]].drop_duplicates()
ap3_results = pd.read_csv(gsi_path / 'local' / 'ap3_results' / 'DetectabilitiesOfPeptides.txt', delimiter = "\t")

output = pd.merge(proteins, ap3_results, left_on = ["accession"], right_on = ["Protein id"])[["Peptide sequence", "protein_id", "Peptide detectability"]]
output.rename(columns = {"Peptide detectability": "Prob", "Peptide sequence": "peptide"}, inplace = True)
output["U"] = output["peptide"].str.contains("U")
output["B"] = output["peptide"].str.contains("B")
output["Z"] = output["peptide"].str.contains("Z")
output["exclude seq"] = output[["U","B","Z"]].any(axis = 1)
output.drop(output[output["exclude seq"]].index, inplace = True)
output[["peptide", "protein_id", "Prob"]].to_csv(gsi_path / 'data' / 'digestion' / 'output_file.csv', index = False)