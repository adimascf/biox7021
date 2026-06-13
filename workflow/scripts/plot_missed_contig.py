import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Compiling full and partial missed contig data...")

# colourblind-friendly palette from colour universal design (CUD)
named_colors = {
    "black": "#000000",
    "orange": "#e69f00",
    "skyblue": "#56b4e9",
    "vermilion": "#d55e00",
    "bluish green": "#009e73",
    "yellow": "#f0e442",
    "blue": "#0072b2",
    "reddish purple": "#cc79a7",
}
cud_palette = list(named_colors.values())

def cud(n: int = len(cud_palette), start: int = 0) -> List[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]

# Compile the Data on the fly
records = []
input_files = snakemake.input.tsv if isinstance(snakemake.input.tsv, list) else [snakemake.input.tsv]

for f in input_files:
    path = Path(f)
    model = path.parent.name                 
    depth = path.parent.parent.name          
    combo = path.parent.parent.parent.name   
    sample = path.name.split('.')[0]         

    try:
        df = pd.read_csv(f, sep='\t')
        
        if not df.empty and 'coverage' in df.columns:
            # get the sequence name column '#rname'
            rname_col = '#rname' if '#rname' in df.columns else df.columns[0]

            full_mask = df['coverage'] == 0
            partial_mask = (df['coverage'] > 0) & (df['coverage'] < 50)

            full_missed = full_mask.sum()
            partial_missed = partial_mask.sum()

            # Extract names of the contigs that meet the missed criteria
            missed_names = df[full_mask | partial_mask][rname_col].tolist()
            missed_names_str = "; ".join(missed_names) if missed_names else ""
        else:
            full_missed, partial_missed = 0, 0
            missed_names_str = ""
            
    except pd.errors.EmptyDataError:
        full_missed, partial_missed = 0, 0
        missed_names_str = ""

    records.append({
        "sample": sample,
        "combo": combo,
        "depth": depth,
        "model": model,
        "full_missed": full_missed,
        "partial_missed": partial_missed,
        "total_missed": full_missed + partial_missed,
        "missed_contig_names": missed_names_str  
    })

# Convert to pandas DataFrame and save the compiled table
missed_df = pd.DataFrame(records)
missed_df.to_csv(snakemake.output.table, index=False)
print(f"Saved compiled data table to {snakemake.output.table}")

# Plot 1: 2x2 Grid (Individual Sample Points)
melted_df = missed_df.melt(
    id_vars=["sample", "combo", "depth", "model"],
    value_vars=["full_missed", "partial_missed"],
    var_name="miss_type",
    value_name="count"
)

sns.set_theme(style="whitegrid")
models = ["sup", "hac"]
hue_order = ["100x", "50x", "20x"]
miss_types = ["full_missed", "partial_missed"]
miss_labels = {"full_missed": "Full Miss (0% Cov)", "partial_missed": "Partial Miss (<50% Cov)"}

# Sort combos by lowest mean total missed contigs (Best to Worst overall)
combo_perf = missed_df.groupby("combo")["total_missed"].mean()
order_abs = combo_perf.sort_values(ascending=True).index.tolist()

fig, axes = plt.subplots(nrows=2, ncols=len(models), figsize=(14, 10), dpi=300, sharex=True, sharey=True)

for r, m_type in enumerate(miss_types):
    for c, model in enumerate(models):
        ax = axes[r, c]
        df_sub = melted_df.query("model == @model and miss_type == @m_type").copy()

        if not df_sub.empty:
            draw_legend = (r == 0 and c == len(models) - 1) 
            sns.stripplot(
                data=df_sub, x="combo", y="count", hue="depth", 
                order=order_abs, hue_order=hue_order,
                palette=cud(len(hue_order), start=2), 
                ax=ax, alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=draw_legend
            )

        ax.yaxis.set_major_formatter(plt.ScalarFormatter())
        
        if c == 0:
            ax.set_ylabel(miss_labels[m_type])
        else:
            ax.set_ylabel("")

        if r == 0:
            ax.set_title(f"{model.upper()} Model")
        
        ax.set_xlabel("")
        ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

        if r == 1:
            ax.set_xticks(range(len(order_abs)))
            ax.set_xticklabels(order_abs, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
        else:
            ax.tick_params(labelbottom=False)

        if r == 0 and c == len(models) - 1 and ax.get_legend() is not None:
            ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')

fig.tight_layout()
fig.savefig(snakemake.output.figures, bbox_inches='tight')
plt.close(fig)
print(f"Saved original 2x2 plot to {snakemake.output.figures}")


# Plot 2: 1x2 Grid (Pooled across the samples and missed contig types)
# Calculate the sum of total missed contigs across all samples for each tool combo/depth
pooled_df = missed_df.groupby(["combo", "model", "depth"])["total_missed"].sum().reset_index()

fig2, axes2 = plt.subplots(nrows=1, ncols=len(models), figsize=(14, 6), dpi=300, sharex=True, sharey=True)

# Safety check for single axes
if not isinstance(axes2, (np.ndarray, list)):
    axes2 = [axes2]

for c, model in enumerate(models):
    ax = axes2[c]
    df_sub = pooled_df.query("model == @model").copy()

    if not df_sub.empty:
        # sns.barplot(
        #     data=df_sub, x="combo", y="total_missed", hue="depth", 
        #     order=order_abs, hue_order=hue_order,
        #     palette=cud(len(hue_order), start=2), 
        #     ax=ax, edgecolor="black", linewidth=0.8
        # )

        sns.stripplot(
            data=df_sub, x="combo", y="total_missed", hue="depth", 
            order=order_abs, hue_order=hue_order,
            palette=cud(len(hue_order), start=2), 
            ax=ax, alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black"
        )

    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    
    if c == 0:
        ax.set_ylabel("Total Missed Contigs")
    else:
        ax.set_ylabel("")

    ax.set_title(f"Total Missed Contigs - {model.upper()} Model")
    ax.set_xlabel("")
    ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

    ax.set_xticks(range(len(order_abs)))
    ax.set_xticklabels(order_abs, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)

    # Manage legends
    if c == len(models) - 1 and ax.get_legend() is not None:
        ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
    elif ax.get_legend() is not None:
        ax.get_legend().remove()

fig2.tight_layout()
fig2.savefig(snakemake.output.figures_total, bbox_inches='tight')
plt.close(fig2)
print(f"Saved pooled totals plot to {snakemake.output.figures_total}")
