import pandas as pd

# === CONFIGURATION ===
input_file = "result/raw_output.csv"
output_file = "result/bcell_linear_epitopes.csv"
threshold = 0.1512  # cutoff for epitope prediction
# =====================

# Load CSV (handles tabs, commas, or mixed delimiters)
df = pd.read_csv(input_file, sep=None, engine="python")

# Clean up column names
df.columns = [c.strip() for c in df.columns]

# Try to automatically identify columns
res_col = None
score_col = None
accession_col = None

for c in df.columns:
    if "Residue" in c:
        res_col = c
    if "linear" in c and "score" in c.lower():
        score_col = c
    if "Accession" in c:
        accession_col = c

if res_col is None or score_col is None:
    raise ValueError(f"Could not find residue/score columns. Found columns: {df.columns}")

print(f"Detected columns -> Residue: {res_col}, Score: {score_col}")

# Group by accession if available
if accession_col:
    groups = df.groupby(accession_col)
else:
    groups = [(None, df)]

epitope_records = []

for acc, group in groups:
    residues = group[res_col].astype(str).tolist()
    scores = group[score_col].astype(float).tolist()

    current_epitope = ""
    start = None
    current_scores = []

    for i, (aa, score) in enumerate(zip(residues, scores), start=1):
        if score >= threshold:
            if current_epitope == "":
                start = i
            current_epitope += aa
            current_scores.append(score)
        else:
            if current_epitope != "":
                epitope_records.append({
                    "Accession": acc if acc else "unknown",
                    "Epitope": current_epitope,
                    "Start": start,
                    "End": i - 1,
                    "Length": len(current_epitope),
                    "Mean Score": round(sum(current_scores) / len(current_scores), 6)
                })
                current_epitope = ""
                current_scores = []

    # Handle last open region
    if current_epitope != "":
        epitope_records.append({
            "Accession": acc if acc else "unknown",
            "Epitope": current_epitope,
            "Start": start,
            "End": len(residues),
            "Length": len(current_epitope),
            "Mean Score": round(sum(current_scores) / len(current_scores), 6)
        })

# Save results
epi_df = pd.DataFrame(epitope_records)
epi_df.to_csv(output_file, index=False)

print(f"Done! Found {len(epi_df)} linear epitopes (score ≥ {threshold}).")
print(f"Results saved to: {output_file}\n")
print(epi_df.head(10))
