from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import sys

import pandas as pd
counts = pd.read_csv("counts_df.csv", index_col=0)
clinical = pd.read_csv("clinical_df.csv", index_col=0)

dds = DeseqDataSet(
    counts=counts,
    clinical=clinical,
    design_factors="condition",
    refit_cooks=False,
    n_cpus=1
)

dds.deseq2()
pd.DataFrame({"gene": dds.varm["genewise_dispersions"], "fitted": dds.varm["fitted_dispersions"], "map": dds.varm["MAP_dispersions"]}).to_csv("fitted_disp_table.csv")

dds.uns["trend_coeffs"].to_csv("fitted_trend_coefs.csv")

# printing an ArrayView to a file
original_stdout = sys.stdout
f = open("fitted_disp_prior_var.txt", "w")
sys.stdout = f
print(dds.uns["prior_disp_var"])
sys.stdout = original_stdout
f.close()

stat_res = DeseqStats(dds, n_cpus=1)
stat_res.summary()

stat_res.results_df.to_csv("fitted_results_table.csv")

stat_res.lfc_shrink(coeff="condition_B_vs_A")
stat_res.results_df.to_csv("fitted_shrunken_lfc.csv")
