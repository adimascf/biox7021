import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

input_file=snakemake.input.csv
output_file=snakemake.output.csv

df = pd.read_csv(input_file)

summary = df.groupby(["trimmer", "VAR_TYPE"])[["PREC", "RECALL", "F1_SCORE"]].mean().reset_index()

summary = summary.round(4)

summary.to_csv(output_file, index=False)
