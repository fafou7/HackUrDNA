#!/usr/bin/env python3
from pathlib import Path
import ebony_classifier as ec  # import your existing classifier module

ALIGNMENT = Path("ebony_training_plus_ref_aligned.fasta")
MODEL = Path("ebony_model_v2.json")

def main():
    if not ALIGNMENT.exists():
        raise SystemExit(f"Alignment file not found: {ALIGNMENT}")
    if not MODEL.exists():
        raise SystemExit(f"Model file not found: {MODEL}")

    # reuse functions from ebony_classifier.py
    seqs = ec.read_fasta(ALIGNMENT)
    model = ec.load_model(MODEL)

    if not seqs:
        raise SystemExit("No sequences found in alignment.")

    # take the LAST sequence in the alignment = your reference
    seq_id, seq = list(seqs.items())[-1]

    result = ec.score_sequence(seq, model)

    print(f"Scoring sequence (from alignment): {seq_id}")
    print(f"Informative sites used : {result['usable_sites']}")
    print(f"Matches dark alleles   : {result['matches_dark']}")
    print(f"Matches light alleles  : {result['matches_light']}")
    if result["darkness_score"] is None:
        print("Darkness score         : N/A (no usable sites)")
    else:
        print(f"Darkness score (0=light,1=dark): {result['darkness_score']:.3f}")

if __name__ == "__main__":
    main()
