# %%
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
root_dir = Path.cwd().parent
data_dir = root_dir / "data"
sol_dir = root_dir / "solution"
# %%
N_prot = 25008
ref_df = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";", index_col = 0)
sol_df = [pd.read_csv(sol_dir / (f"results_yeast_10fmol0.{i}00000.csv"), index_col = 0) for i in range(1,10)]
#%%
N_prot = 25008
ref_df = pd.read_csv(data_dir / 'YEAST-Data-NonNormalized.csv', sep = ";", index_col = 0)
sol_df = pd.read_csv(sol_dir / "results_yeast_2fmol100.csv", index_col = 0)
# %%
ref_series = ref_df['110616_yeast_ups_10fmol']
sol_series = sol_df[-1]['abundance']
# %%
ref_series = ref_df['110714_yeast_ups1_2fmol_r3']
sol_series = sol_df['abundance']
# %%
def str_to_float(string : str) -> float:
    return float(string.replace(',', '.'))

ref_series = ref_series.apply(str_to_float)
# %%
plt.hist(sol_series.to_numpy(), bins = np.logspace(np.log10(1), np.log10(9), 50))
# %%
sorted_ref = ref_series.sort_values(ascending = False)
sorted_sol = sol_series.sort_values(ascending = False)
# %%
pears_corr = [ref_series.corr(df['abundance'], 'pearson') for df in sol_df]
spear_corr = [ref_series.corr(df['abundance'], 'spearman') for df in sol_df]
# %%
pears_corr = ref_series.corr(sol_series, 'pearson')
spear_corr = ref_series.corr(sol_series, 'spearman')
# %%
sol_df[0].join(sol_df[-1], how = 'inner', rsuffix='2')
# %%
for i in range(9):
    print(len(ref_df.join(sol_df[i], how = 'inner')))
# %% Stats
TP = len(pd.merge(ref_series, sol_series, left_index = True, right_index = True, how = 'inner'))
FP = len(sol_series) - TP
FN = len(ref_series) - TP
TN = N_prot - TP - FP - FN

accuracy = (TP + TN) / N_prot
specificity = TN/ (TN + FN)
sensitivity = TP / (TP + FP)
print(f"""True positives : {TP}
True negatives : {TN}
False positives : {FP}
False negatives : {FN}
Accuracy : {round(accuracy, 3)}
Sensitivity : {round(sensitivity, 3)}
Specificity : {round(specificity, 3)}
Pearson correlation : {round(pears_corr, 3)}
Spearman correlation : {round(spear_corr, 3)}""")
# %%
missed_lines = [sorted_ref.index.get_loc(df.iloc[-1].name) - df['abundance'].size for df in sorted_sol]
# %%
