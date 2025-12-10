from pathlib import Path

# ---- CONFIG ----
genome_path = Path("dmel-all-chromosome-r6.65.fasta")  # << use this one
chrom_hint = "3R"      # we want the 3R arm; we'll search for a header containing this
start = 8390000        # 1-based inclusive
end   = 8423000        # 1-based inclusive
output_fasta = Path("ref_ebony.fasta")
# -----------------

def find_chrom_name(genome_path, hint):
    """
    Scan FASTA headers and return the first header name that matches the hint.
    Matching rules:
      - exact match
      - header ends with hint
      - hint ends with header
      - header contains hint as substring
    """
    with open(genome_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if (
                    name == hint
                    or name.endswith(hint)
                    or hint.endswith(name)
                    or hint in name
                ):
                    print(f"Using chromosome header: {name}")
                    return name
    return None

def main():
    if not genome_path.exists():
        raise FileNotFoundError(genome_path)

    chrom_name = find_chrom_name(genome_path, chrom_hint)
    if chrom_name is None:
        raise SystemExit(
            f"Could not find any FASTA header matching hint '{chrom_hint}'. "
            "Open the FASTA and check the exact chromosome names."
        )

    # collect sequence for that chromosome
    seq_chunks = []
    current = None
    with open(genome_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
            elif current == chrom_name:
                seq_chunks.append(line)

    genome_seq = "".join(seq_chunks).upper()
    print(f"Chromosome {chrom_name} length = {len(genome_seq)} bp")

    if len(genome_seq) == 0:
        raise SystemExit("No sequence collected for that chromosome. Check genome file.")

    # clip requested window (handle if end goes past sequence length)
    if start < 1 or start > len(genome_seq):
        raise SystemExit(f"Start {start} outside chromosome length {len(genome_seq)}")

    end_clipped = min(end, len(genome_seq))
    subseq = genome_seq[start - 1 : end_clipped]

    with open(output_fasta, "w") as out:
        out.write(f">{chrom_name}:{start}-{end_clipped}\n")
        for i in range(0, len(subseq), 80):
            out.write(subseq[i : i + 80] + "\n")

    print(f"Extracted {len(subseq)} bp to {output_fasta}")

if __name__ == "__main__":
    main()
