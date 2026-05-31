import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Redirect both stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Processing {len(snakemake.input.pr)} precision-recall files...")

# Load all TSV files
frames = []
for p in map(Path, snakemake.input.pr):
    df = pd.read_csv(p, sep="\t")
    # Path: .../<combo>/<depth>x/<model>/<sample>.precision-recall.tsv
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

# Copy relevant columns
data_fnfp = data[['combo', 'depth', 'sample', 'model', 'VAR_TYPE', 'TRUTH_FN', 'QUERY_FP', 'TRUTH_TOTAL', "F1_SCORE", "PREC", "RECALL"]].copy()
data_fnfp = data_fnfp[(data_fnfp["VAR_TYPE"] == "SNP") | (data_fnfp["VAR_TYPE"] == "INDEL")]

# Sorting maps
depth_order = {"100x": 1, "50x": 2, "20x": 3}
data_fnfp['depth_sort'] = data_fnfp['depth'].map(depth_order)

model_order = {"sup": 1, "hac": 2}
data_fnfp['model_sort'] = data_fnfp['model'].map(model_order).fillna(3)

trimmer_order = {t: i for i, t in enumerate(sorted(data_fnfp['combo'].unique())) if t != 'unprocessed-untrimmed'}
trimmer_order['unprocessed-untrimmed'] = 999 # force to the end
data_fnfp['trimmer_sort'] = data_fnfp['combo'].map(trimmer_order)

# Sort the data_fnfp
data_fnfp = data_fnfp.sort_values(
    by=['sample', 'model_sort', 'depth_sort', 'VAR_TYPE', 'trimmer_sort']
).drop(columns=['depth_sort', 'model_sort', 'trimmer_sort'])

data_fnfp.reset_index(drop=True, inplace=True)

data_fnfp.to_csv(snakemake.output.csv, index=False)
print(f"Saved compiled FN/FP metrics to {snakemake.output.csv}")

# Data Processing: Calculate Deltas
df_plot = data_fnfp.copy()

# Isolate the unprocessed-untrimmed baseline
baseline = df_plot[df_plot['combo'] == 'unprocessed-untrimmed'].set_index(['sample', 'depth', 'model', 'VAR_TYPE'])[['TRUTH_FN', 'QUERY_FP']]

# Join the baseline back to the main dataframe
df_merged = df_plot.set_index(['sample', 'depth', 'model', 'VAR_TYPE']).join(baseline, rsuffix='_baseline').reset_index()

# Calculate the Deltas (baseline - combo)
# Positive value = tool PREVENTED errors compared to baseline
df_merged['FN_Prevented'] = df_merged['TRUTH_FN_baseline'] - df_merged['TRUTH_FN']
df_merged['FP_Prevented'] = df_merged['QUERY_FP_baseline'] - df_merged['QUERY_FP']

# Remove the unprocessed-untrimmed rows
df_final = df_merged[df_merged['combo'] != 'unprocessed-untrimmed'].copy()

# Plotting Setup
sns.set_theme(style="whitegrid")

# Map metrics to their respective output files
metrics_to_plot = {
    "FN_Prevented": snakemake.output.fn_plot,
    "FP_Prevented": snakemake.output.fp_plot
}

rows = ["SNP", "INDEL"]
models = ["sup", "hac"]
hue_order = ["100x", "50x", "20x"]
order = sorted(df_final['combo'].unique())

# Loop through both metrics and generate a separate 2x2 grid for each
for metric_col, output_file in metrics_to_plot.items():
    
    fig, axes = plt.subplots(
        nrows=len(rows), ncols=len(models), figsize=(14, 10), dpi=300, sharex=True, sharey='row'
    )

    for r, vartype in enumerate(rows):
        for c, model in enumerate(models):
            ax = axes[r, c]
            df_sub = df_final.query("VAR_TYPE == @vartype and model == @model")

            # Draw the zero-baseline (represents 'unprocessed-untrimmed')
            ax.axhline(0, color='red', linestyle='--', linewidth=1.5, zorder=0, label="Untrimmed Baseline")

            if not df_sub.empty:
                sns.boxplot(
                    data=df_sub, x="combo", y=metric_col, hue="depth", 
                    order=order, hue_order=hue_order, ax=ax, fliersize=0, gap=0.2
                )
                sns.stripplot(
                    data=df_sub, x="combo", y=metric_col, hue="depth", 
                    order=order, hue_order=hue_order, ax=ax, 
                    alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=False
                )

            # Clean formatting
            ylabel_clean = metric_col.replace('_', ' ')
            ax.set_ylabel(f"{vartype} {ylabel_clean}" if c == 0 else "")
            ax.set_title(f"{model} model" if r == 0 else "")
            ax.set_xlabel("")
            
            # X-axis ticks only on bottom row
            if r == len(rows) - 1:
                ax.tick_params(axis="x", labelrotation=45, labelsize=10)
            else:
                ax.tick_params(axis="x", bottom=False)

            # Legend handling (only top right)
            if r == 0 and c == 1 and ax.get_legend() is not None:
                ax.legend(title="Depth", loc="upper left", bbox_to_anchor=(1.02, 1))
            elif ax.get_legend() is not None:
                ax.get_legend().remove()

    fig.tight_layout()
    fig.savefig(output_file, bbox_inches='tight')
    print(f"Saved {metric_col} plot to {output_file}")
