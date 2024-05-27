import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("ap3_result_file")
parser.add_argument("digestion_file")
args = parser.parse_args()

digestion = pd.read_csv(Path.cwd().parent / 'data' / 'digestion' / f"{args.digestion_file}.csv")[["peptide","accession","protein_id","peptide_id"]]
ap3_res = pd.read_csv(Path.cwd().parent / 'data' / 'digestion' / f"{args.ap3_result_file}.txt", delimiter = '\t')

merged = pd.merge(digestion, ap3_res, left_on = ["peptide", "accession"], right_on = ["Peptide sequence", "Protein id"]).drop_duplicates()
merged.rename({"Peptide detectability": "Prob"})
merged.to_csv(Path.cwd().parent / 'data' / 'digestion'/ f"{args.ap3_result_file}_result.csv", index = False)