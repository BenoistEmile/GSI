import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main(args):
    score_dir = Path(args.score_dir)
    input_file_name = args.input_file_name
    FDR_rate = args.FDR_rate
    input_df = pd.read_csv(score_dir / f"{input_file_name}.csv", delimiter=";")
    input_df.sort_values(["score after specfit"], inplace=True, ascending=False)
    input_df.reset_index(inplace=True, drop=True)
    last_FDR = 0
    FDR_cutoff_index = 0
    for index, row in input_df.iterrows(): # Computation of the FDR for every PSM
        if row["presence in target database"] == "target":
            last_FDR = last_FDR*index / (index + 1)
            input_df.loc[index, "FDR"] = last_FDR
        else:
            last_FDR = (last_FDR*index + 1) / (index + 1)
            input_df.loc[index, "FDR"] = last_FDR
        if last_FDR > FDR_rate and FDR_cutoff_index == 0:
            FDR_cutoff_index = index - 1
    input_df.loc[input_df["presence in decoy database"] == "target"][:FDR_cutoff_index].to_csv(score_dir / f"{input_file_name}_FDR_{FDR_rate}.csv", sep=";")

    fig, ax = plt.subplots(1, 1)
    (n, bins, patches) = ax.hist(input_df.loc[input_df["presence in decoy database"] == "target", "score after specfit"], bins=75, alpha=0.5, label="True")
    ax.hist(input_df.loc[input_df["presence in decoy database"] == "decoy", "score after specfit"], bins=bins, alpha=0.5, label="False")
    ax.vlines(input_df.loc[FDR_cutoff_index, "score after specfit"], 0, 1500, label="FDR cutoff")
    ax.set_title("Scores distributions of True and False identifications")
    ax.set_xlabel("Scores")
    ax.set_ylabel("Number of identifications")
    ax.legend()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--score_dir', type=str, default="C:\\Users\\berthier-a-1\\code\\GSI\\GSI\\data\\scores", help='Path of the score directory')
    parser.add_argument('--input_file_name', type=str, default="specoms_output_HeLa_human", help="Name of the input file (without extension)")
    parser.add_argument("--FDR_rate", type=float, default=0.01, help="FDR cutoff")
    args = parser.parse_args()
    print(args)
    main(args)
