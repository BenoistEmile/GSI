# %%
from pathlib import Path
import pandas as pd

# %%
root_dir = Path.cwd().parent
data_dir = root_dir / "data"
sol_dir = root_dir / "solution"

# %%
ref_df = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";", index_col = 0)
sol_df = [pd.read_csv(sol_dir / ("results_yeast0." + str(i) + "00000.csv"), index_col = 0) for i in range(1,10)]
# %%
def str_to_float(string : str) -> float:
    return float(string.replace(',', '.'))

ref_df['110714_yeast_ups1_2fmol_r3'] = ref_df['110714_yeast_ups1_2fmol_r3'].apply(str_to_float)

# %%
sorted_ref = ref_df.sort_values('110714_yeast_ups1_2fmol_r3', ascending = False)['110714_yeast_ups1_2fmol_r3']
sorted_sol = [df.sort_values('abundance', ascending = False) for df in sol_df]
# %%
missed_lines = [sorted_ref.index.get_loc(df.iloc[-1].name) - df['abundance'].size for df in sorted_sol]
# %%
pears_corr = [ref_df['110714_yeast_ups1_2fmol_r3'].corr(df['abundance'], 'pearson') for df in sol_df]
spear_corr = [ref_df['110714_yeast_ups1_2fmol_r3'].corr(df['abundance'], 'spearman') for df in sol_df]
# %%
sol_df[0].join(sol_df[1], how = 'inner', rsuffix='2')
# %%
for i in range(9):
    print(len(ref_df.join(sol_df[i], how = 'inner')))
# %%
