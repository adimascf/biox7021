import sys
import os
import pandas as pd
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")

data_frames = []
summary_records = []
col_names = ['assembly_file', 'barcode_id', 'contig', 'start', 'end', 'cost', 'source_file']

print(f"Aggregating {len(snakemake.input.contam_report)} contamination files...")

for sample_report in snakemake.input.contam_report:
    p = Path(sample_report)
    
    # Path: .../contaminant/<trimmer>/<depth>/<sample>/<file.tsv>
    trimmer = p.parts[-4]
    depth = p.parts[-3]
    sample = p.parts[-2]

    # If the file is 0 bytes, record a count of 0
    if os.path.getsize(p) == 0:
        summary_records.append({
            'trimmer': trimmer, 'depth': depth, 'sample': sample, 'contamination_count': 0
        })
    else:
        df = pd.read_csv(p, sep='\t', names=col_names)
        df['trimmer'] = trimmer
        df['depth'] = depth
        df['sample'] = sample
        
        data_frames.append(df)
        
        # Record the actual length of the dataframe as the count
        summary_records.append({
            'trimmer': trimmer, 'depth': depth, 'sample': sample, 'contamination_count': len(df)
        })

if data_frames:
    combined_df = pd.concat(data_frames, ignore_index=True)
    final_cols = ['trimmer', 'depth', 'sample'] + col_names
    combined_df = combined_df[final_cols]
    combined_df.to_csv(snakemake.output.details, index=False)
    print(f"Saved detailed results to: {snakemake.output.details}")

else:
    # no contamination is found anywhere
    pd.DataFrame(columns=['trimmer', 'depth', 'sample'] + col_names).to_csv(snakemake.output.details, index=False)
    print("No detailed contaminations found. Empty template created.")


if summary_records:
    summary_df = pd.DataFrame(summary_records)
    
    # apply custom sorting logic
    depth_order = {"100x": 1, "50x": 2, "20x": 3}
    summary_df['depth_sort'] = summary_df['depth'].map(depth_order)
    
    trimmer_order = {t: i for i, t in enumerate(sorted(summary_df['trimmer'].unique())) if t != 'untrimmed'}
    trimmer_order['untrimmed'] = 999 # force to the end
    summary_df['trimmer_sort'] = summary_df['trimmer'].map(trimmer_order)
    
    # Sort the summary_df: sample -> depth -> trimmer
    summary_df = summary_df.sort_values(
        by=['sample', 'depth_sort', 'trimmer_sort']
    ).drop(columns=['depth_sort', 'trimmer_sort'])
    
    summary_df.reset_index(drop=True, inplace=True)
    
    # Save as CSV
    summary_df.to_csv(snakemake.output.summary, index=False)
    print(f"Saved summary counts to: {snakemake.output.summary}")

print("Aggregation complete.")
