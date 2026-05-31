from pathlib import Path
import pandas as pd

# redirect both stderr and stdout to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")


input_files=snakemake.input.pr
output_file=snakemake.output.csv

frames = []
for p in map(Path, input_files):
    df = pd.read_csv(p, sep="\t")
    # .../<combo>/<depth>/<model>/<sample>.file.tsv
    df["model"] = p.parts[-2]
    df["depth"] = p.parts[-3]
    df["combo"] = p.parts[-4]
    df["sample"] = p.name.split('.')[0]
    frames.append(df)

pr_df = pd.concat(frames)
pr_df.reset_index(inplace=True, drop=True)

# Get the row with the maximum F1 score for each combination
dataix = pr_df.groupby(["combo", "depth", "VAR_TYPE", "sample", "model"])["F1_SCORE"].idxmax()
data = pr_df.iloc[dataix]

# order it to make it easy to compare combo performance
data_fnfp = data[['combo', 'depth', 'sample', 'model', 'VAR_TYPE', 'TRUTH_FN', 'QUERY_FP', 'TRUTH_TOTAL']].copy()
data_fnfp = data_fnfp[(data_fnfp["VAR_TYPE"] == "SNP") | (data_fnfp["VAR_TYPE"] == "INDEL")]

# Sort mappings
depth_order = {"100x": 1, "50x": 2, "20x": 3}
data_fnfp['depth_sort'] = data_fnfp['depth'].map(depth_order)

model_order = {"sup": 1, "hac": 2}
data_fnfp['model_sort'] = data_fnfp['model'].map(model_order).fillna(3)

trimmer_order = {t: i for i, t in enumerate(sorted(data_fnfp['combo'].unique())) if t != 'unprocessed-untrimmed'}
trimmer_order['unprocessed-untrimmed'] = 999 # force to the end
data_fnfp['trimmer_sort'] = data_fnfp['combo'].map(trimmer_order)

# sort the data_fnfp: sample -> model -> depth -> variant type -> combo
data_fnfp = data_fnfp.sort_values(
    by=['sample', 'model_sort', 'depth_sort', 'VAR_TYPE', 'trimmer_sort']
).drop(columns=['depth_sort', 'model_sort', 'trimmer_sort'])

# reset the index so it looks clean
data_fnfp.reset_index(drop=True, inplace=True)
data_fnfp.to_csv(output_file, index=False)

print(f"Extracted FN/FP numbers and saved to {output_file}")
