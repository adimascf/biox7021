import pandas as pd

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Compiling all assembly metrics")

quast_metrics = pd.read_csv(snakemake.input.quast) 
contam_metrics = pd.read_csv(snakemake.input.contaminants) 
miscontig_metrics = pd.read_csv(snakemake.input.missed_contigs) 

assembly_metrics = quast_metrics.merge(
    contam_metrics, on=['combo', 'depth', 'sample', 'model'], how='inner'
).merge(
    miscontig_metrics, on=['combo', 'depth', 'sample', 'model'], how='inner'
)
assembly_metrics.to_csv(snakemake.output.master, index=False)
