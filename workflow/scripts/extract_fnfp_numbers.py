from pathlib import Path
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


input_files=snakemake.input.pr
output_file=snakemake.output.csv

frames = []
for p in map(Path, input_files):
    df = pd.read_csv(p, sep="\t")
    df["sample"] = p.parts[-2]
    df["depth"] = p.parts[-3]
    df["trimmer"] = p.parts[-4]
    frames.append(df)

pr_df = pd.concat(frames)
pr_df.reset_index(inplace=True, drop=True)

metrics = ["TRUTH_FN", "QUERY_FP"]
x = "trimmer"
hue = "depth"
rows = ("SNP", "INDEL")
nrows = len(rows)

# Get the row with the maximum F1 score for each combination
dataix = pr_df.groupby([x, hue, "VAR_TYPE", "sample"])["F1_SCORE"].idxmax()
data = pr_df.iloc[dataix]

# order it to make it easy to compare trimmer performance
data_fnfp = data[['trimmer', 'depth', 'sample', 'VAR_TYPE', 'TRUTH_FN', 'QUERY_FP', 'TRUTH_TOTAL']].copy()
data_fnfp = data_fnfp[(data_fnfp["VAR_TYPE"] == "SNP") | (data_fnfp["VAR_TYPE"] == "INDEL")]
depth_order = {"100x": 1, "50x": 2, "20x": 3}
data_fnfp['depth_sort'] = data_fnfp['depth'].map(depth_order)

trimmer_order = {t: i for i, t in enumerate(sorted(data_fnfp['trimmer'].unique())) if t != 'untrimmed'}
trimmer_order['untrimmed'] = 999 # force to the end
data_fnfp['trimmer_sort'] = data_fnfp['trimmer'].map(trimmer_order)

# sort the data_fnfp: sample -> depth -> variant type -> trimmer
data_fnfp = data_fnfp.sort_values(
    by=['sample', 'depth_sort', 'VAR_TYPE', 'trimmer_sort']
).drop(columns=['depth_sort', 'trimmer_sort'])

# reset the index so it looks clean
data_fnfp.reset_index(drop=True, inplace=True)
data_fnfp.to_csv(output_file, index=False)
