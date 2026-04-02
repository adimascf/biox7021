# originally from Michael Hall https://github.com/mbhall88/NanoVarBench/blob/main/workflow/scripts/best_f1.py

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sys.stderr = open(snakemake.log[0], "w")

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
    # In results/assess/mutref/{trimmer}/{sample}/file.tsv
    # p.parts[-2] is the sample, p.parts[-3] is the trimmer
    df["sample"] = p.parts[-2]
    df["depth"] = p.parts[-3]
    df["trimmer"] = p.parts[-4]
    frames.append(df)

pr_df = pd.concat(frames)
pr_df.reset_index(inplace=True, drop=True)

metrics = ["F1_SCORE", "PREC", "RECALL"]
x = "trimmer"
hue = "depth"
rows = ("SNP", "INDEL")
nrows = len(rows)

# Get the row with the maximum F1 score for each combination
dataix = pr_df.groupby([x, hue, "VAR_TYPE", "sample"])["F1_SCORE"].idxmax()
data = pr_df.iloc[dataix]
order = sorted(set(data[x]))
hue_order = ["100x", "50x", "20x"]

# Generate the plots
for y in metrics:
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=1,
        figsize=(12, 10),
        dpi=300,
        sharex=True,
        sharey=True,
    )

    i = 0
    for vartype in rows:
        ax = axes.flatten()[i]
        df = data.query("VAR_TYPE == @vartype").copy()

        cap = 0.99999
        df.loc[:, y] = df[y].apply(lambda v: cap if v > cap else v)
        yticks = [0.01, 0.1, 0.5, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
        yticklabels = [f"{yval:.2%}" for yval in yticks]

        sns.boxplot(
            data=df, x=x, y=y, hue=hue, order=order, hue_order=hue_order, 
            ax=ax, fliersize=0, gap=0.2
        )
        sns.stripplot(
            data=df, x=x, y=y, hue=hue, order=order, hue_order=hue_order, ax=ax,
            alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=False
        )

        ax.set_yscale("logit", nonpositive="clip")
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)

        ylabel = {"F1_SCORE": "F1 score", "PREC": "Precision", "RECALL": "Recall"}[y]
        ax.set_ylabel(f"{vartype} {ylabel}")
        ax.tick_params(axis="x", labelrotation=45, labelsize=10)
        ax.set_xlabel("")
        ax.set_title(f"{vartype} {ylabel}")

        # only show the legend on the top plot to save space
        if i == 0:
            ax.legend(title="Depth", loc="lower right")
        else:
            if ax.get_legend() is not None:
                ax.get_legend().remove()

        i += 1

    fig.tight_layout()
    fig.savefig(output[y])

data = data.query("VAR_TYPE not in ('ALL', 'SV')")
data.sort_values(by=["trimmer", "sample"], inplace=True)
data.to_csv(snakemake.output.csv, index=False)
