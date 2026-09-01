import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from typing import List, Tuple

# Redirect stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Generating sample-level aggregate assembly performance figures and tables (Global and Per-depth normalisations)...")

# File paths from snakemake
MASTER_CSV = snakemake.input.master_csv
SURVEY_CSV = snakemake.input.survey_csv

FIG_GLOBAL_SUP = getattr(snakemake.output, 'fig_global_sup', None) or snakemake.output.get('fig_sup')
FIG_GLOBAL_HAC = getattr(snakemake.output, 'fig_global_hac', None) or snakemake.output.get('fig_hac')
FIG_PERDEPTH_SUP = getattr(snakemake.output, 'fig_perdepth_sup', None)
FIG_PERDEPTH_HAC = getattr(snakemake.output, 'fig_perdepth_hac', None)

SCORES_CSV_GLOBAL = getattr(snakemake.output, 'scores_csv_global', None) or snakemake.output.get('scores_csv')
SUMMARY_CSV_GLOBAL = getattr(snakemake.output, 'summary_csv_global', None) or snakemake.output.get('summary_csv')
SCORES_CSV_PERDEPTH = getattr(snakemake.output, 'scores_csv_perdepth', None)
SUMMARY_CSV_PERDEPTH = getattr(snakemake.output, 'summary_csv_perdepth', None)

# Colourblind-friendly palette from colour universal design (CUD)
named_colours = {
    "black": "#000000",
    "orange": "#e69f00",
    "skyblue": "#56b4e9",
    "vermilion": "#d55e00",
    "bluish green": "#009e73",
    "yellow": "#f0e442",
    "blue": "#0072b2",
    "reddish purple": "#cc79a7",
}
cud_palette = list(named_colours.values())

def cud(n: int = len(cud_palette), start: int = 0) -> List[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]

def calculate_scores(df_input: pd.DataFrame, group_cols: List[str], weights: Tuple[float, float, float, float, float],
                     aunga_col: str, missed_col: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mean_wa, mean_wc, mean_wd, mean_wr, total_w = weights
    sample_score_records = []

    for _, group in df_input.groupby(group_cols):
        g = group.copy()

        # Score 1: Contiguity / Structural Accuracy (target ratio = 1.0)
        g['s_contiguity'] = np.maximum(0.0, 1.0 - np.abs(g[aunga_col] - 1.0)) * 100.0

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

    summary_df = scores_df.groupby(["model", "depth", "combo"], as_index=False, observed=True).agg(
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

    summary_df["rank"] = summary_df.groupby(["model", "depth"])["mean_score"].rank(ascending=False, method="min").astype(int)
    summary_df.sort_values(by=["model", "depth", "rank"], inplace=True)

    return scores_df, summary_df

def plot_global_figure(scores_df: pd.DataFrame, model: str, hue_order: List[str], palette: List[str], out_path: str):
    df_model = scores_df.query("model == @model").copy()
    if df_model.empty:
        print(f"Warning: No data for model {model}, skipping global figure.")
        return

    order = (
        df_model.groupby("combo")["composite_score"]
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    fig, ax = plt.subplots(figsize=(14, 7), dpi=300)

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
        f"Survey-Weighted Assembly Score (Global Normalisation) - {model.upper()} Model",
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
    print(f"Saved global normalisation figure for {model} to: {out_path}")

def plot_perdepth_faceted_figure(scores_df: pd.DataFrame, model: str, hue_order: List[str], palette: List[str], out_path: str):
    df_model = scores_df.query("model == @model").copy()
    if df_model.empty:
        print(f"Warning: No data for model {model}, skipping per-depth figure.")
        return

    n_depths = len(hue_order)
    fig, axes = plt.subplots(nrows=1, ncols=n_depths, figsize=(10 * n_depths, 7), dpi=300, sharey=True)
    if n_depths == 1:
        axes = [axes]

    for d_idx, depth_val in enumerate(hue_order):
        ax = axes[d_idx]
        df_d = df_model.query("depth == @depth_val").copy()
        if df_d.empty:
            continue

        # Sort combos independently for this depth panel (best to worst)
        order_d = (
            df_d.groupby("combo")["composite_score"]
            .mean()
            .sort_values(ascending=False)
            .index.tolist()
        )

        depth_colour = palette[d_idx]

        sns.stripplot(
            data=df_d,
            x="combo",
            y="composite_score",
            order=order_d,
            color=depth_colour,
            alpha=0.5,
            dodge=False,
            linewidth=0.5,
            edgecolor="black",
            zorder=1,
            size=5,
            ax=ax
        )

        sns.pointplot(
            data=df_d,
            x="combo",
            y="composite_score",
            order=order_d,
            color=depth_colour,
            errorbar=('ci', 95),
            capsize=0.1,
            err_kws={'linewidth': 1.2},
            linewidth=1.2,
            markersize=6,
            estimator=np.mean,
            zorder=2,
            ax=ax
        )

        ax.set_title(f"{depth_val} Sequencing Depth", fontsize=13, pad=10, fontweight="bold")
        if d_idx == 0:
            ax.set_ylabel("Composite Performance Score", fontsize=12)
        else:
            ax.set_ylabel("")

        ax.set_xlabel("")
        ax.set_xticks(range(len(order_d)))
        ax.set_xticklabels(order_d, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
        ax.set_ylim(0, 105)
        ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)
        ax.yaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

    fig.suptitle(
        f"Survey-Weighted Assembly Score (Within-Depth Relative Ranking) - {model.upper()} Model",
        fontsize=15,
        fontweight="bold",
        y=1.02
    )

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved per-depth faceted figure for {model} to: {out_path}")

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

    weights = (mean_wa, mean_wc, mean_wd, mean_wr, total_w)
    print(f"Community weights: Accuracy={mean_wa:.2f}, Contiguity={mean_wc:.2f}, Decontam={mean_wd:.2f}, Replicon={mean_wr:.2f}")

    # 3. Determine active depths dynamically
    if hasattr(snakemake, 'config') and 'depth' in snakemake.config:
        config_depths = snakemake.config['depth']
        hue_order = [f"{str(d).rstrip('x')}x" for d in sorted([int(str(d).rstrip('x')) for d in config_depths], reverse=True)]
    else:
        unique_depths = df["depth"].unique()
        depth_ints = [int(str(d).replace('x', '')) for d in unique_depths]
        hue_order = [f"{d}x" for d in sorted(depth_ints, reverse=True)]

    df["depth"] = df["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df = df[df["depth"].isin(hue_order)].copy()
    print(f"Active depths: {hue_order}")

    palette = cud(len(hue_order), start=2)

    aunga_col = "auNGA_ratio" if "auNGA_ratio" in df.columns else ("auNGA_norm" if "auNGA_norm" in df.columns else "auNGA")
    missed_col = "total_missed" if "total_missed" in df.columns else "full_missed"

    # 4. Version 1: Global / Cross-Depth Normalisation (grouped by model, sample)
    scores_df_global, summary_df_global = calculate_scores(
        df, group_cols=['model', 'sample'], weights=weights, aunga_col=aunga_col, missed_col=missed_col
    )
    if SCORES_CSV_GLOBAL:
        Path(SCORES_CSV_GLOBAL).parent.mkdir(parents=True, exist_ok=True)
        scores_df_global.to_csv(SCORES_CSV_GLOBAL, index=False)
        print(f"Saved global sample-level scores table to: {SCORES_CSV_GLOBAL}")
    if SUMMARY_CSV_GLOBAL:
        Path(SUMMARY_CSV_GLOBAL).parent.mkdir(parents=True, exist_ok=True)
        summary_df_global.to_csv(SUMMARY_CSV_GLOBAL, index=False)
        print(f"Saved global summary ranking table to: {SUMMARY_CSV_GLOBAL}")

    # 5. Version 2: Within-Depth Relative Normalisation (grouped by model, depth, sample)
    scores_df_perdepth, summary_df_perdepth = calculate_scores(
        df, group_cols=['model', 'depth', 'sample'], weights=weights, aunga_col=aunga_col, missed_col=missed_col
    )
    if SCORES_CSV_PERDEPTH:
        Path(SCORES_CSV_PERDEPTH).parent.mkdir(parents=True, exist_ok=True)
        scores_df_perdepth.to_csv(SCORES_CSV_PERDEPTH, index=False)
        print(f"Saved per-depth sample-level scores table to: {SCORES_CSV_PERDEPTH}")
    if SUMMARY_CSV_PERDEPTH:
        Path(SUMMARY_CSV_PERDEPTH).parent.mkdir(parents=True, exist_ok=True)
        summary_df_perdepth.to_csv(SUMMARY_CSV_PERDEPTH, index=False)
        print(f"Saved per-depth summary ranking table to: {SUMMARY_CSV_PERDEPTH}")

    # 6. Generate Global Figures
    if FIG_GLOBAL_SUP:
        plot_global_figure(scores_df_global, "sup", hue_order, palette, FIG_GLOBAL_SUP)
    if FIG_GLOBAL_HAC:
        plot_global_figure(scores_df_global, "hac", hue_order, palette, FIG_GLOBAL_HAC)

    # 7. Generate Per-Depth Faceted Figures (1x2 grid with independent sorting)
    if FIG_PERDEPTH_SUP:
        plot_perdepth_faceted_figure(scores_df_perdepth, "sup", hue_order, palette, FIG_PERDEPTH_SUP)
    if FIG_PERDEPTH_HAC:
        plot_perdepth_faceted_figure(scores_df_perdepth, "hac", hue_order, palette, FIG_PERDEPTH_HAC)

    print("All aggregate assembly assessment processing complete.")

if __name__ == "__main__":
    main()
