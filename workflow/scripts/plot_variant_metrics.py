# originally from Michael Hall https://github.com/mbhall88/NanoVarBench/blob/main/workflow/scripts/best_f1.py

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# redirect both stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

sns.set_theme(style="whitegrid")

output = {}
for path in snakemake.output.figures:
    p = Path(path)
    metric = p.stem.split("_")[-1] # extracts 'f1', 'recall', or 'precision'
    if metric == "f1":
        metric = "F1_SCORE"
    elif metric == "precision":
        metric = "PREC"
    elif metric == "recall":
        metric = "RECALL"
    output[metric] = p

# Load all TSV files 
frames = []
for p in map(Path, snakemake.input.pr):
    df = pd.read_csv(p, sep="\t")
    # Path: ../<trimmer>/<depth>/<sample>/<model>/<file.tsv>
    df["model"] = p.parts[-2]
    df["sample"] = p.parts[-3]
    df["depth"] = p.parts[-4]
    df["trimmer"] = p.parts[-5]
    frames.append(df)

pr_df = pd.concat(frames)
pr_df.reset_index(inplace=True, drop=True)

metrics = ["F1_SCORE", "PREC", "RECALL"]
vartypes = ["SNP", "INDEL"]
models = ["sup", "hac"]

# Group by the new 'model' column as well to get the max F1 per specific condition
dataix = pr_df.groupby(["trimmer", "depth", "VAR_TYPE", "sample", "model"])["F1_SCORE"].idxmax()
data = pr_df.iloc[dataix].copy()

order = sorted(set(data["trimmer"]))
hue_order = ["100x", "50x", "20x"]

# Generate the plots
for y in metrics:
    # Set up a 2x2 grid: rows are VAR_TYPE, cols are MODEL
    fig, axes = plt.subplots(
        nrows=len(vartypes),
        ncols=len(models),
        figsize=(14, 10), # Slightly wider to accommodate two columns comfortably
        dpi=300,
        sharex=True,
        sharey=True,
    )

    # Loop through rows (SNP, INDEL) and columns (sup, hac)
    for r, vartype in enumerate(vartypes):
        for c, model in enumerate(models):
            ax = axes[r, c]
            
            # Filter data for this specific row/col combination
            df = data.query("VAR_TYPE == @vartype and model == @model").copy()

            cap = 0.99999
            df.loc[:, y] = df[y].apply(lambda v: cap if v > cap else v)
            yticks = [0.01, 0.1, 0.5, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
            yticklabels = [f"{yval:.2%}" for yval in yticks]

            # Only plot if data exists for this combination
            if not df.empty:
                sns.boxplot(
                    data=df, x="trimmer", y=y, hue="depth", 
                    order=order, hue_order=hue_order, 
                    ax=ax, fliersize=0, gap=0.2
                )
                sns.stripplot(
                    data=df, x="trimmer", y=y, hue="depth", 
                    order=order, hue_order=hue_order, ax=ax,
                    alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=False
                )

            # Format Y-axis (Shared across columns)
            ax.set_yscale("logit", nonpositive="clip")
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)

            ylabel = {"F1_SCORE": "F1 score", "PREC": "Precision", "RECALL": "Recall"}[y]
            
            # Only put Y-axis labels on the leftmost column
            if c == 0:
                ax.set_ylabel(f"{vartype} {ylabel}")
            else:
                ax.set_ylabel("")

            # Only put Titles on the top row
            if r == 0:
                ax.set_title(f"{model} model", fontsize=14)
            else:
                ax.set_title("")

            # X-axis formatting (Labels only on bottom row)
            ax.set_xlabel("")
            if r == len(vartypes) - 1:
                ax.tick_params(axis="x", labelrotation=45, labelsize=10)
            else:
                ax.tick_params(axis="x", bottom=False)

            # Legend Handling: Show it only once, on the top-right plot
            if r == 0 and c == 1:
                if ax.get_legend() is not None:
                    # Move legend outside the plot to avoid overlapping data points
                    ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
            else:
                if ax.get_legend() is not None:
                    ax.get_legend().remove()

    fig.tight_layout()
    fig.savefig(output[y], bbox_inches='tight')

# Save the subsetted data to CSV
data = data.query("VAR_TYPE not in ('ALL', 'SV')")
data.sort_values(by=["trimmer", "sample", "model", "depth"], inplace=True)
data.to_csv(snakemake.output.csv, index=False)
