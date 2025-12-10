#!/usr/bin/env python3
import argparse
import json
from collections import Counter

# ---------------------
# Basic FASTA / CSV I/O
# ---------------------

def read_fasta(path):
    """
    Simple FASTA reader.
    Returns dict: {id: sequence}

    Uses the first whitespace-separated token after '>' as ID
    and strips any version suffix like '.1':
      '>EF114371.1 Drosophila...' -> ID = 'EF114371'
    """
    seqs = {}
    current_id = None
    current_seq = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # save previous
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq).upper()
                header = line[1:].strip()
                raw_id = header.split()[0]
                base_id = raw_id.split(".")[0]  # strip version
                current_id = base_id
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            seqs[current_id] = "".join(current_seq).upper()

    return seqs


def read_labels(path):
    """
    VERY ROBUST label reader.

    Expects a CSV-like file where:
      - first line is header (ignored)
      - each subsequent line has at least two columns separated by commas
      - first column = ID (e.g. 'EF114371')
      - last column  = phenotype ('dark' or 'light')

    Works for:
      id,phenotype
      EF114371,dark

    and also:
      accession,isolate,population,phenotype
      EF114371,U76,Uganda,dark
    """
    labels = {}
    with open(path, "r") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("labels.csv is empty")

    # skip header
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",") if p.strip() != ""]
        if len(parts) < 2:
            continue
        rid = parts[0]
        phenotype = parts[-1].lower()  # last column
        if phenotype not in {"dark", "light"}:
            raise ValueError(f"Phenotype must be 'dark' or 'light', got '{phenotype}' in line: {line}")
        labels[rid] = phenotype

    if not labels:
        raise ValueError("No labels found in labels.csv")
    return labels


def split_dark_light(seqs, labels):
    dark = {}
    light = {}
    missing = []

    for sid, seq in seqs.items():
        if sid in labels:
            if labels[sid] == "dark":
                dark[sid] = seq
            else:
                light[sid] = seq
        else:
            missing.append(sid)

    if missing:
        print("Warning: these sequences had no label and were ignored:")
        for m in missing:
            print("  ", m)

    return dark, light


def check_same_length(seqs):
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(f"Sequences have different lengths: {lengths}. "
                         "Align them first (e.g. with MAFFT).")
    return next(iter(lengths))

# ---------------------
# Model building
# ---------------------

def build_model(dark_seqs, light_seqs, min_freq=0.6, ignore_gaps=True):
    """
    Build a simple model:
      - For each alignment column, look at alleles in dark vs light.
      - If major allele in dark != major allele in light,
        and each has freq >= min_freq, mark it as informative.

    Returns:
      dict with:
        'length': alignment length
        'positions': list of {pos, dark_allele, light_allele}
    """
    if not dark_seqs or not light_seqs:
        raise ValueError("Need at least one dark and one light sequence.")

    all_seqs = {}
    all_seqs.update(dark_seqs)
    all_seqs.update(light_seqs)
    align_len = check_same_length(all_seqs)

    dark_list = list(dark_seqs.values())
    light_list = list(light_seqs.values())

    informative = []

    for i in range(align_len):
        col_dark = [s[i] for s in dark_list]
        col_light = [s[i] for s in light_list]

        if ignore_gaps:
            col_dark = [b for b in col_dark if b not in "-Nn"]
            col_light = [b for b in col_light if b not in "-Nn"]

        if not col_dark or not col_light:
            continue

        d_counts = Counter(col_dark)
        l_counts = Counter(col_light)

        d_allele, d_count = d_counts.most_common(1)[0]
        l_allele, l_count = l_counts.most_common(1)[0]

        d_freq = d_count / len(col_dark)
        l_freq = l_count / len(col_light)

        if d_allele != l_allele and d_freq >= min_freq and l_freq >= min_freq:
            informative.append({
                "pos": i,
                "dark_allele": d_allele,
                "light_allele": l_allele,
            })

    return {
        "length": align_len,
        "positions": informative,
    }


def save_model(model, path):
    with open(path, "w") as f:
        json.dump(model, f, indent=2)


def load_model(path):
    with open(path, "r") as f:
        return json.load(f)

# ---------------------
# Scoring
# ---------------------

def score_sequence(seq, model):
    seq = seq.upper()
    if len(seq) != model["length"]:
        raise ValueError(
            f"Sequence length {len(seq)} != model length {model['length']}. "
            "Align to the same reference as training."
        )

    usable = 0
    match_dark = 0
    match_light = 0

    for site in model["positions"]:
        pos = site["pos"]
        base = seq[pos]
        if base in "-Nn":
            continue
        usable += 1
        if base == site["dark_allele"]:
            match_dark += 1
        elif base == site["light_allele"]:
            match_light += 1

    if usable == 0:
        return {
            "usable_sites": 0,
            "matches_dark": 0,
            "matches_light": 0,
            "darkness_score": None,
        }

    darkness_score = match_dark / usable

    return {
        "usable_sites": usable,
        "matches_dark": match_dark,
        "matches_light": match_light,
        "darkness_score": darkness_score,
    }

# ---------------------
# CLI
# ---------------------

def main():
    parser = argparse.ArgumentParser(
        description="Train ebony dark/light model and score sequences."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # train
    p_train = subparsers.add_parser("train", help="Train model from aligned FASTA + labels.csv")
    p_train.add_argument("--fasta", required=True,
                         help="Aligned multi-FASTA (e.g. ebony_training.fasta)")
    p_train.add_argument("--labels", required=True,
                         help="labels.csv file")
    p_train.add_argument("--out", default="ebony_model.json",
                         help="Output JSON model file (default ebony_model.json)")
    p_train.add_argument("--min-freq", type=float, default=0.6,
                         help="Min major-allele frequency in each group (default 0.6)")

    # score
    p_score = subparsers.add_parser("score", help="Score a new aligned sequence")
    p_score.add_argument("--seq-fasta", required=True,
                         help="FASTA with ONE aligned sequence to score")
    p_score.add_argument("--model", required=True,
                         help="Model JSON from train step")

    args = parser.parse_args()

    if args.command == "train":
        seqs = read_fasta(args.fasta)
        labels = read_labels(args.labels)
        dark, light = split_dark_light(seqs, labels)
        print(f"Loaded {len(dark)} dark and {len(light)} light sequences.")

        model = build_model(dark, light, min_freq=args.min_freq)
        print(f"Found {len(model['positions'])} informative positions.")
        save_model(model, args.out)
        print(f"Model written to {args.out}")

    elif args.command == "score":
        model = load_model(args.model)
        seqs = read_fasta(args.seq_fasta)
        if len(seqs) != 1:
            raise ValueError("score mode expects exactly ONE sequence in --seq-fasta")

        sid, seq = next(iter(seqs.items()))
        res = score_sequence(seq, model)

        print(f"Scoring sequence: {sid}")
        print(f"Informative sites used : {res['usable_sites']}")
        print(f"Matches dark alleles   : {res['matches_dark']}")
        print(f"Matches light alleles  : {res['matches_light']}")
        if res["darkness_score"] is None:
            print("Darkness score         : N/A (no usable sites)")
        else:
            print(f"Darkness score (0=light,1=dark): {res['darkness_score']:.3f}")


if __name__ == "__main__":
    main()
