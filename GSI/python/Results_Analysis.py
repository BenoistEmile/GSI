#%%
from pathlib import Path
import pandas as pd

#%%
root_dir = Path.cwd().parent
data_dir = root_dir / "data"
sol_dir = root_dir / "solution"

#%%
ref_df = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";", index_col = 0)
sol_df = [pd.read_csv(sol_dir / ("results_yeast0." + str(i) + "00000.csv"), index_col = 0) for i in range(1,10)]
# %%
