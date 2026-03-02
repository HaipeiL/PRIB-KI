import pandas as pd
from pathlib import Path
import random

# ----------------------------
# Config
# ----------------------------
INFILE = Path("data/TheraSAbDab.csv")
OUTFILE = Path("data/demo_100.csv")

N_REAL = 60
N_VAR = 40
N_TOTAL = N_REAL + N_VAR

# random seed for reproducibility (change if you want a different shuffle)
SEED = 42

# ----------------------------
# Helpers
# ----------------------------
AA = set("ACDEFGHIKLMNPQRSTVWY")

def is_protein_sequence(s: str) -> bool:
    """Basic sanity check: length + mostly standard amino acids."""
    if not isinstance(s, str):
        return False
    s = s.strip().upper()
    if len(s) < 80:
        return False
    return sum(ch in AA for ch in s) / len(s) > 0.95

AA_HYDRO = list("AILMFWVY")  # hydrophobic aa (rough)
AA_POLAR = list("STNQ")      # polar aa
AA_CHARGED = list("DEKRH")   # charged aa

def make_stress_variant(seq: str, n_mut: int = 12) -> str:
    """
    Create a 'stress-test' variant to enrich higher developability risk profiles.
    This is NOT a claim of experimental truth; it is for demonstrator stratification.
    Heuristics:
      - increase local hydrophobicity (aggregation/viscosity proxy)
      - disturb charge/polar pattern (solubility/self-interaction proxy)
    """
    s = list(seq.strip().upper())
    L = len(s)
    if L == 0:
        return seq

    idxs = random.sample(range(L), k=min(n_mut, L))
    for i in idxs:
        r = random.random()
        if r < 0.60:
            # push toward hydrophobic residues
            s[i] = random.choice(AA_HYDRO)
        elif r < 0.85:
            # perturb with polar residues
            s[i] = random.choice(AA_POLAR)
        else:
            # if charged, neutralize; else increase hydrophobicity
            if s[i] in AA_CHARGED:
                s[i] = random.choice(AA_POLAR)
            else:
                s[i] = random.choice(AA_HYDRO)
    return "".join(s)

# ----------------------------
# Main
# ----------------------------
def main():
    random.seed(SEED)

    if not INFILE.exists():
        raise FileNotFoundError(f"Input file not found: {INFILE.resolve()}")

    df = pd.read_csv(INFILE)

    print("Loaded:", INFILE.resolve())
    print("Rows:", len(df))
    print("\nAll columns:")
    for c in df.columns:
        print(" -", c)

    # Find sequence-like columns and choose the best one
    candidate_columns = []
    for c in df.columns:
        sample = df[c].dropna().astype(str).head(300)
        score = sum(is_protein_sequence(x) for x in sample)
        if score >= 5:
            candidate_columns.append((c, score))

    if not candidate_columns:
        raise RuntimeError("No suitable sequence column found automatically. Please inspect columns printed above.")

    candidate_columns.sort(key=lambda x: x[1], reverse=True)
    seq_col = candidate_columns[0][0]
    print(f"\nUsing sequence column: {seq_col}  (score={candidate_columns[0][1]})")

    # Extract and clean sequences
    tmp = df[[seq_col]].dropna().copy()
    tmp["sequence"] = tmp[seq_col].astype(str).str.strip().str.upper()
    tmp = tmp[tmp["sequence"].apply(is_protein_sequence)]
    tmp = tmp.drop_duplicates("sequence").reset_index(drop=True)

    if len(tmp) < N_REAL:
        raise RuntimeError(f"Not enough valid unique sequences found ({len(tmp)}) to sample {N_REAL} real candidates.")

    # Sample real sequences
    real = tmp.sample(n=N_REAL, random_state=SEED).reset_index(drop=True)
    real["source"] = "public_therasabdab"

    # Create stress variants from the real set (or from tmp pool)
    variants = []
    for j in range(N_VAR):
        base_seq = real.loc[j % N_REAL, "sequence"]
        vseq = make_stress_variant(base_seq, n_mut=12)
        variants.append({"sequence": vseq, "source": "stress_variant"})

    var_df = pd.DataFrame(variants)

    # Combine and shuffle (important: do not expose category in ID order)
    out_df = pd.concat([real[["sequence", "source"]], var_df[["sequence", "source"]]], ignore_index=True)

    # Optional: ensure variants are not accidentally identical to some real sequence
    # (rare, but let's be safe)
    out_df = out_df.drop_duplicates("sequence").reset_index(drop=True)
    if len(out_df) < N_TOTAL:
        # If duplicates reduced count, top up by generating more variants
        needed = N_TOTAL - len(out_df)
        print(f"Warning: duplicates removed; topping up with {needed} additional variants.")
        extra = []
        base_pool = real["sequence"].tolist()
        for _ in range(needed):
            base_seq = random.choice(base_pool)
            extra.append({"sequence": make_stress_variant(base_seq, n_mut=14), "source": "stress_variant"})
        out_df = pd.concat([out_df, pd.DataFrame(extra)], ignore_index=True)
        out_df = out_df.drop_duplicates("sequence").reset_index(drop=True)

    # Final shuffle
    out_df = out_df.sample(frac=1.0, random_state=SEED).reset_index(drop=True)

    # Assign neutral IDs (no group leakage)
    out_df.insert(0, "id", [f"Cand_{i:03d}" for i in range(1, len(out_df) + 1)])

    # Write output
    OUTFILE.parent.mkdir(exist_ok=True)
    # Keep 'source' column for your internal checks, but you can hide it in the UI/video if you want
    out_df[["id", "sequence", "source"]].to_csv(OUTFILE, index=False)

    print(f"\nWrote {len(out_df)} candidates to {OUTFILE.resolve()}")
    print(out_df["source"].value_counts(dropna=False))

if __name__ == "__main__":
    main()