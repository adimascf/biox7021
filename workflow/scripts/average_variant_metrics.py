import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

input_file=snakemake.input.csv
output_file=snakemake.output.csv

df = pd.read_csv(input_file)

summary = df.groupby(["trimmer", "depth", "VAR_TYPE"])[["PREC", "RECALL", "F1_SCORE"]].mean().reset_index()
summary = summary.round(5)
# sort 100x, 50x, and 20x
summary["depth"] = pd.Categorical(summary["depth"], categories=["100x", "50x", "20x"], ordered= True)
summary = summary.sort_values(by=["trimmer", "depth", "VAR_TYPE"])

summary.to_csv(output_file, index=False)
