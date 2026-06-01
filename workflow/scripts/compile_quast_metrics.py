import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Pre-calculating genome sizes from FASTA indices...")
genome_sizes = {}

# Ensure it's treated as a list even if there's only one sample
fai_files = snakemake.input.fai if isinstance(snakemake.input.fai, list) else [snakemake.input.fai]

for fai_file in fai_files:
    fai_path = Path(fai_file)
    # Extracts the sample name (everything before the first '.')
    sample_name = fai_path.name.split('.')[0]
    
    try:
        with open(fai_path, "r") as f:
            # The second column (index 1) in a .fai file is the sequence length
            genome_sizes[sample_name] = sum(int(line.split("\t")[1]) for line in f)
    except FileNotFoundError:
        print(f"Warning: .fai file not found: {fai_file}")

assembly_metrics = []
for file_path in snakemake.input.reports:
    try:
        p = Path(file_path)
        
        # Path: ../quast/<combo>/<depth>x/<model>/<sample>.<combo>.report.tsv
        model = p.parts[-2]
        depth = p.parts[-3]
        combo = p.parts[-4]
        
        # Extract the sample name from the file name
        sample = p.name.split('.')[0]

        df_quast = pd.read_csv(p, sep='\t', index_col=0)
        col_name = df_quast.columns[0]

        mismatches = float(df_quast.loc['# mismatches per 100 kbp', col_name])
        indels = float(df_quast.loc['# indels per 100 kbp', col_name])
        nga50 = int(df_quast.loc['NGA50', col_name])
        misassemblies = int(df_quast.loc['# misassemblies', col_name])

        genome_size = genome_sizes.get(sample, np.nan)
        nga50_norm = nga50 / genome_size if pd.notna(genome_size) else np.nan
        
        assembly_metrics.append({
            'combo': combo, 'depth': depth, 'sample': sample, 'model': model,
            'Mismatches per 100kbp': mismatches,
            'Indels per 100kbp': indels,
            'NGA50': nga50,
            'NGA50_norm': nga50_norm,
            'misassemblies': misassemblies
        })

    except Exception as e:
        print(f"Skipping {file_path} due to error: {e}")

df = pd.DataFrame(assembly_metrics)
df.sort_values(by=["combo", "sample", "model", "depth"], inplace=True)
df.to_csv(snakemake.output.csv, index=False)
print(f"Saved compiled metrics to {snakemake.output.csv}")
