import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sys.stderr = open(snakemake.log[0], "w")

frames = []
for p in map(Path, snakemake.input.qends):
    df = pd.read_csv(p, sep="\t")
    # extract wildcard from Path
    df["sample"] = p.stem.split(".")[0]
    df["trimmer"] = p.parts[-2]
    frames.append(df)

df = pd.concat(frames)
df.reset_index(inplace=True, drop=True)

df_melt = df.melt(
        id_vars=["Position", "sample", "trimmer"],
        value_vars=["Start_MedianQ", "End_MedianQ"],
        var_name="End_Type", value_name="Median_Quality")

# for labeling
df_melt["End_Type"] = df_melt["End_Type"].replace({
    "Start_MedianQ": "Read Start (First 150bp)",
    "End_MedianQ": "Read End (Last 150bp)"})

# Generating the plot
sns.set_theme(style="whitegrid")
# I am using relationship plot phred score to base position
# https://seaborn.pydata.org/generated/seaborn.relplot.html
plot = sns.relplot(
        data=df_melt,
        x="Position",
        y="Median_Quality",
        hue="trimmer",
        col="End_Type",
        kind="line",
        height=6,
        aspect=1.2,
        linewidth=2,
        facet_kws={'sharex': False, 'sharey': True})

plot.set_axis_labels("Position from Read End (bp)", "Phred Quality Score (Median)")
plot.set_titles("{col_name}")
plot.set(ylim=(0, 50)) # cap at Q50

# move legends to the bottom
sns.move_legend(plot, 
                "lower center", 
                bbox_to_anchor=(0.45, -0.10), 
                ncol=3, title=None, frameon=False)
plot.figure.subplots_adjust(bottom=0.2)

plt.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")

