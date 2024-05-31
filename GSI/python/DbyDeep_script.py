import argparse
from pathlib import Path
import pandas as pd
import dbydeep_model

if __name__ == "__main__":
    gsi_path = Path.cwd()
    opt = argparse.Namespace(retrain_flag = False,
                              data_path = gsi_path / 'data' / 'digestion' / 'digestion_file.csv',
                              model_path = gsi_path / 'python' / 'DbyDeep.h5',
                              save_path = (gsi_path / 'data' / 'digestion').absolute().as_posix() + '/',
                              job_name = "output_DbyDeep")
    print(opt)
    # dbydeep_model.main(opt)

    detectabilities = pd.read_csv(gsi_path / 'data' / 'digestion' / 'output_DbyDeep.csv')[["peptide", "protein_id", "Prob"]]
    detectabilities.to_csv(gsi_path / 'data' / 'digestion' / 'output_file.csv', index = False)