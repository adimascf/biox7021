import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_rel, shapiro
from scipy.special import logit
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# Redirect all prints and errors to the Snakemake log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

baseline_combo = "unprocessed-untrimmed" 
out_dir = Path(snakemake.output[0]).parent
out_dir.mkdir(parents=True, exist_ok=True)

# Define metrics, their directionality, and the required data transformation
metric_properties = {
    "F1_SCORE": {"higher_is_better": True, "transform": "logit"},
    "PREC": {"higher_is_better": True, "transform": "logit"},
    "RECALL": {"higher_is_better": True, "transform": "logit"},
    "NGA50_norm": {"higher_is_better": True, "transform": "logit"},
    "Mismatches per 100kbp": {"higher_is_better": False, "transform": "log1p"},
    "Indels per 100kbp": {"higher_is_better": False, "transform": "log1p"}
}

print(f"Loading data from {snakemake.input.csv}...")
df = pd.read_csv(snakemake.input.csv)

if "VAR_TYPE" in df.columns:
    group_cols = ["VAR_TYPE", "model", "depth"]
    df = df.query("VAR_TYPE not in ('ALL', 'SV')").copy()
else:
    group_cols = ["model", "depth"]
    if "NGA50_norm" in df.columns:
        dataix = df.groupby(["combo", "depth", "sample", "model"])["NGA50_norm"].idxmax()
        df = df.loc[dataix].copy()

for metric, props in metric_properties.items():
    if metric not in df.columns:
        continue 
        
    higher_is_better = props["higher_is_better"]
    transform_type = props["transform"]
        
    print(f"\n--- Processing {metric} (Transform: {transform_type}) ---")
    results = []
    groups = df.groupby(group_cols)
    
    for group_keys, group_data in groups:
        group_dict = dict(zip(group_cols, group_keys if isinstance(group_keys, tuple) else (group_keys,)))
        
        baseline_data = group_data[group_data["combo"] == baseline_combo]
        if baseline_data.empty:
            continue
            
        for combo in group_data["combo"].unique():
            if combo == baseline_combo:
                continue 
                
            target_data = group_data[group_data["combo"] == combo]
            merged = pd.merge(
                baseline_data[["sample", metric]], 
                target_data[["sample", metric]], 
                on="sample", 
                suffixes=('_base', '_target')
            )
            
            if len(merged) < 2:
                p_value = np.nan
                median_diff = np.nan
                direction = "Unknown"
            else:
                base_raw = merged[f"{metric}_base"]
                target_raw = merged[f"{metric}_target"]
                
                diff = target_raw - base_raw
                median_diff = diff.median()
                
                if median_diff == 0:
                    direction = "No Change"
                elif (median_diff > 0 and higher_is_better) or (median_diff < 0 and not higher_is_better):
                    direction = "Improved"
                else:
                    direction = "Worsened"
                
                if np.all(diff == 0):
                    p_value = 1.0 
                else:
                    try:
                        if transform_type == "logit":
                            epsilon = 1e-5
                            base_trans = logit(np.clip(base_raw, epsilon, 1 - epsilon))
                            target_trans = logit(np.clip(target_raw, epsilon, 1 - epsilon))
                        elif transform_type == "log1p":
                            base_trans = np.log1p(base_raw)
                            target_trans = np.log1p(target_raw)
                        
                        # Calculate differences of the transformed data for normality checking
                        trans_diff = target_trans - base_trans
                        
                        # Shapiro-Wilk requires at least 3 data points
                        if len(trans_diff) >= 3:
                            stat_shapiro, p_shapiro = shapiro(trans_diff)
                            normality_status = "Normal (p > 0.05)" if p_shapiro > 0.05 else "Not Normal (p <= 0.05)"
                            print(f"Normality Check [{combo} vs {baseline_combo} | {group_keys}]: p={p_shapiro:.4f} -> {normality_status}")
                        else:
                            print(f"Normality Check [{combo} vs {baseline_combo} | {group_keys}]: Skipped (N < 3)")

                        # Paired t-test
                        stat, p_value = ttest_rel(base_trans, target_trans)
                    except ValueError:
                        p_value = np.nan
                        
            result_row = {
                "combo": combo,
                **group_dict,
                "median_diff": median_diff,
                "direction": direction,
                "p_value": p_value
            }
            results.append(result_row)
            
    if not results:
        continue
        
    results_df = pd.DataFrame(results)
    
    valid_mask = results_df['p_value'].notna()
    results_df['p_adj_BH'] = np.nan 
    
    if valid_mask.any():
        _, p_adj, _, _ = multipletests(results_df.loc[valid_mask, 'p_value'], method='fdr_bh')
        results_df.loc[valid_mask, 'p_adj_BH'] = p_adj
        
    results_df['significant'] = results_df['p_adj_BH'] < 0.05
    
    sort_cols = group_cols + ["p_adj_BH"]
    results_df = results_df.sort_values(by=sort_cols)
    
    metric_safe_name = metric.replace(" ", "_").replace("/", "_").lower()
    output_file = out_dir / f"combo_{metric_safe_name}_significance.csv"
    results_df.to_csv(output_file, index=False)
    print(f"Saved statistical results to {output_file.name}")

print("Statistical testing complete.")
