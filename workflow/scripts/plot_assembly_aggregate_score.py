import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from typing import List

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Generating sample-level aggregate assembly performance figures and tables...")

# File paths from snakemake
MASTER_CSV = snakemake.input.master_csv
SURVEY_CSV = snakemake.input.survey_csv

FIG_SUP = snakemake.output.fig_sup
FIG_HAC = snakemake.output.fig_hac
SCORES_CSV = snakemake.output.scores_csv
SUMMARY_CSV = snakemake.output.summary_csv

# Colourblind-friendly palette from colour universal design (CUD)
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

def main():
    sns.set_theme(style="whitegrid")

    # 1. Load data
    df = pd.read_csv(MASTER_CSV)
    survey_df = pd.read_csv(SURVEY_CSV)

    # 2. Calculate community average survey weights
    mean_wa = float(survey_df["w_accuracy"].mean())
    mean_wc = float(survey_df["w_contiguity"].mean())
    mean_wd = float(survey_df["w_decontam"].mean())
    mean_wr = float(survey_df["w_replicon"].mean())
    total_w = mean_wa + mean_wc + mean_wd + mean_wr
    if total_w == 0:
        total_w = 1e-5

    print(f"Community weights: Accuracy={mean_wa:.2f}, Contiguity={mean_wc:.2f}, Decontam={mean_wd:.2f}, Replicon={mean_wr:.2f}")

    # 3. Determine active depths dynamically
    if hasattr(snakemake, 'config') and 'depth' in snakemake.config:
        config_depths = snakemake.config['depth']
        hue_order = [f"{str(d).rstrip('x')}x" for d in sorted([int(str(d).rstrip('x')) for d in config_depths], reverse=True)]
    else:
        unique_depths = df["depth"].unique()
        depth_ints = [int(str(d).replace('x', '')) for d in unique_depths]
        hue_order = [f"{d}x" for d in sorted(depth_ints, reverse=True)]

    # Standardize depth format in df (ensure 'x' suffix)
    df["depth"] = df["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")

    # Filter for active depths
    df = df[df["depth"].isin(hue_order)].copy()
    print(f"Active depths: {hue_order}")

    palette = cud(len(hue_order), start=2)

    # 4. Compute Sample-Level Normalized Scores and Composite Score
    aunga_col = "auNGA_ratio" if "auNGA_ratio" in df.columns else ("auNGA_norm" if "auNGA_norm" in df.columns else "auNGA")
    missed_col = "total_missed" if "total_missed" in df.columns else "full_missed"

    sample_score_records = []

    # Normalize per (model, depth, sample) across all tested tool combinations
    for (model, depth, sample), group in df.groupby(['model', 'depth', 'sample']):
        g = group.copy()

        # Score 1: Contiguity / Structural Accuracy (target ratio = 1.0)
        g['s_contiguity'] = np.maximum(0, 1.0 - np.abs(g[aunga_col] - 1.0)) * 100.0

        # Score 2: Sequence Accuracy (Mismatches and Indels, lower is better)
        for err_col in ['Mismatches per 100kbp', 'Indels per 100kbp']:
            mx, mn = g[err_col].max(), g[err_col].min()
            suffix = 'mismatches' if 'Mismatches' in err_col else 'indels'
            g[f's_{suffix}'] = (((mx - g[err_col]) / (mx - mn) * 100.0) if mx != mn else 100.0)
        g['s_accuracy'] = (g['s_mismatches'] + g['s_indels']) / 2.0

        # Score 3: Decontamination (lower contamination count is better)
        mx_con, mn_con = g['contamination_count'].max(), g['contamination_count'].min()
        g['s_decontam'] = (((mx_con - g['contamination_count']) / (mx_con - mn_con) * 100.0) if mx_con != mn_con else 100.0)

        # Score 4: Replicon Completeness (fewer missed contigs is better)
        mx_mis, mn_mis = g[missed_col].max(), g[missed_col].min()
        g['s_replicon'] = (((mx_mis - g[missed_col]) / (mx_mis - mn_mis) * 100.0) if mx_mis != mn_mis else 100.0)

        # Weighted Geometric Mean
        eps = 1e-5
        s_contig_safe = np.where(g['s_contiguity'] <= 0, eps, g['s_contiguity'])
        s_acc_safe = np.where(g['s_accuracy'] <= 0, eps, g['s_accuracy'])
        s_decont_safe = np.where(g['s_decontam'] <= 0, eps, g['s_decontam'])
        s_rep_safe = np.where(g['s_replicon'] <= 0, eps, g['s_replicon'])

        final_scores = (
            (s_contig_safe ** mean_wc) *
            (s_acc_safe ** mean_wa) *
            (s_decont_safe ** mean_wd) *
            (s_rep_safe ** mean_wr)
        ) ** (1.0 / total_w)

        g['composite_score'] = final_scores
        sample_score_records.append(g)

    scores_df = pd.concat(sample_score_records, ignore_index=True)

    # Save detailed sample-level scores table
    Path(SCORES_CSV).parent.mkdir(parents=True, exist_ok=True)
    scores_df.to_csv(SCORES_CSV, index=False)
    print(f"Saved sample-level composite scores ({len(scores_df)} rows) to: {SCORES_CSV}")

    # 5. Generate Summary Ranking Table
    summary_df = scores_df.groupby(["model", "depth", "combo"], as_index=False).agg(
        mean_score=("composite_score", "mean"),
        median_score=("composite_score", "median"),
        std_score=("composite_score", "std"),
        min_score=("composite_score", "min"),
        max_score=("composite_score", "max"),
        sample_count=("sample", "count"),
        s_contiguity_mean=("s_contiguity", "mean"),
        s_accuracy_mean=("s_accuracy", "mean"),
        s_decontam_mean=("s_decontam", "mean"),
        s_replicon_mean=("s_replicon", "mean"),
    )

    # Calculate ranks per model and depth
    summary_df["rank"] = summary_df.groupby(["model", "depth"])["mean_score"].rank(ascending=False, method="min").astype(int)
    summary_df.sort_values(by=["model", "depth", "rank"], inplace=True)
    summary_df.to_csv(SUMMARY_CSV, index=False)
    print(f"Saved summary ranking table to: {SUMMARY_CSV}")

    # 6. Generate Standalone Figures for SUP and HAC
    fig_outputs = {"sup": FIG_SUP, "hac": FIG_HAC}
    num_samples = len(scores_df["sample"].unique())

    for model, out_path in fig_outputs.items():
        df_model = scores_df.query("model == @model").copy()
        if df_model.empty:
            print(f"Warning: No data for model {model}, skipping figure.")
            continue

        # Order combos descending by overall mean score across depths for this model
        order = (
            df_model.groupby("combo")["composite_score"]
            .mean()
            .sort_values(ascending=False)
            .index.tolist()
        )

        fig, ax = plt.subplots(figsize=(14, 7), dpi=300)

        # Background stripplot for the 13 individual biological samples
        sns.stripplot(
            data=df_model,
            x="combo",
            y="composite_score",
            hue="depth",
            order=order,
            hue_order=hue_order,
            palette=palette,
            alpha=0.5,
            dodge=True,
            linewidth=0.5,
            edgecolor="black",
            zorder=1,
            size=5,
            ax=ax
        )

        # Foreground pointplot with mean and 95% CI error bars
        sns.pointplot(
            data=df_model,
            x="combo",
            y="composite_score",
            hue="depth",
            order=order,
            hue_order=hue_order,
            palette=palette,
            dodge=0.3,
            errorbar=('ci', 95),
            capsize=0.1,
            err_kws={'linewidth': 1.2},
            linewidth=1.2,
            markersize=6,
            estimator=np.mean,
            legend=False,
            zorder=2,
            ax=ax
        )

        ax.set_title(
            f"Survey-Weighted Assembly Score - {model.upper()} Model",
            fontsize=14,
            pad=12,
            fontweight="bold"
        )
        ax.set_ylabel("Composite Performance Score", fontsize=12)
        ax.set_xlabel("")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
        ax.set_ylim(0, 105)

        ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)
        ax.yaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

        # Unified depth legend outside top right
        custom_handles = [
            Patch(facecolor=palette[i], edgecolor='black', label=hue_order[i])
            for i in range(len(hue_order))
        ]
        ax.legend(
            handles=custom_handles,
            title="Depth",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=True,
            framealpha=0.9,
            facecolor="white",
            edgecolor="lightgrey",
            fontsize=11,
            title_fontsize=12
        )

        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved figure for {model} to: {out_path}")

    print("All aggregate assembly assessment processing complete.")

if __name__ == "__main__":
    main()
