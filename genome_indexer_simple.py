#!/usr/bin/env python3
"""
Simple genome/FASTA indexer into SQLite.

Supports:
- VCF / VCF.gz
- 23andMe raw text
- FASTA (.fa / .fasta / gzipped)

Usage example (Windows):
    python genome_indexer_simple.py "C:\\Users\\Fafou\\Desktop\\dmel-all-exon-r6.65.fasta" --db dmel_exon.db
"""

import argparse
import gzip
import sqlite3
from pathlib import Path


# -----------------------
# Helpers
# -----------------------

def open_maybe_gzip(path: Path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="ignore")
    return open(p, "rt", encoding="utf-8", errors="ignore")


def guess_format(path: Path) -> str:
    """
    Very small heuristic:
      - .vcf / .vcf.gz          -> vcf
      - .fa / .fasta / .fa.gz   -> fasta
      - if any first lines start with '>' -> fasta
      - if '23andMe' or '# rsid' in header -> 23andme
      - else -> 23andme (default)
    """
    suffix = "".join(path.suffixes).lower()

    if suffix.endswith(".vcf") or suffix.endswith(".vcf.gz"):
        return "vcf"

    if (
        suffix.endswith(".fa")
        or suffix.endswith(".fasta")
        or suffix.endswith(".fa.gz")
        or suffix.endswith(".fasta.gz")
    ):
        return "fasta"

    # Look at first lines
    with open_maybe_gzip(path) as f:
        for _ in range(20):
            line = f.readline()
            if not line:
                break
            if line.startswith(">"):
                return "fasta"
            if "23andMe" in line or line.startswith("# rsid"):
                return "23andme"

    # Default
    return "23andme"


# -----------------------
# Parsers
# -----------------------

def parse_vcf(path: Path, source_name: str):
    with open_maybe_gzip(path) as f:
        for line in f:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue

            rsid = fields[2] if fields[2] not in (".", "") else None
            ref = fields[3]
            alt = fields[4]

            qual_str = fields[5]
            quality = None
            if qual_str not in (".", ""):
                try:
                    quality = float(qual_str)
                except ValueError:
                    quality = None

            genotype = ""
            if len(fields) >= 10:
                fmt = fields[8]
                sample = fields[9]
                fmt_keys = fmt.split(":")
                sample_vals = sample.split(":")
                if "GT" in fmt_keys:
                    gt_idx = fmt_keys.index("GT")
                    if gt_idx < len(sample_vals):
                        genotype = sample_vals[gt_idx]

            yield {
                "chrom": str(chrom),
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "genotype": genotype,
                "rsid": rsid,
                "source_file": source_name,
                "quality": quality,
            }


def parse_23andme(path: Path, source_name: str):
    with open_maybe_gzip(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split()
            if len(fields) < 4:
                continue

            rsid, chrom, pos_str, genotype = fields[:4]
            try:
                pos = int(pos_str)
            except ValueError:
                continue

            yield {
                "chrom": str(chrom),
                "pos": pos,
                "ref": "",
                "alt": "",
                "genotype": genotype,
                "rsid": rsid,
                "source_file": source_name,
                "quality": None,
            }


def parse_fasta(path: Path, source_name: str):
    """
    For each base in each sequence:
      chrom = header without '>'
      pos   = 1-based index
      ref   = base
    """
    with open_maybe_gzip(path) as f:
        chrom = None
        pos = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                chrom = line[1:].split()[0]
                pos = 0
                continue
            if chrom is None:
                continue

            for base in line:
                pos += 1
                yield {
                    "chrom": str(chrom),
                    "pos": pos,
                    "ref": base,
                    "alt": "",
                    "genotype": "",
                    "rsid": None,
                    "source_file": source_name,
                    "quality": None,
                }


# -----------------------
# SQLite
# -----------------------

def init_db(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            genotype TEXT,
            rsid TEXT,
            source_file TEXT,
            quality REAL
        );
        """
    )
    cur.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON variants(rsid);")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_pos ON variants(chrom, pos);")
    conn.commit()
    return conn


# -----------------------
# Indexing logic
# -----------------------

def index_file(input_path: Path, db_path: Path):
    input_path = input_path.expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    db = init_db(db_path)

    fmt = guess_format(input_path)
    print(f"Detected format: {fmt}")

    if fmt == "vcf":
        parser = parse_vcf
    elif fmt == "fasta":
        parser = parse_fasta
    else:
        parser = parse_23andme

    batch = []
    batch_size = 1000
    processed = 0
    cur = db.cursor()

    for rec in parser(input_path, source_name=input_path.name):
        batch.append(
            (
                rec["chrom"],
                rec["pos"],
                rec["ref"],
                rec["alt"],
                rec["genotype"],
                rec["rsid"],
                rec["source_file"],
                rec["quality"],
            )
        )
        processed += 1

        if len(batch) >= batch_size:
            cur.executemany(
                """
                INSERT INTO variants (chrom, pos, ref, alt, genotype, rsid, source_file, quality)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?);
                """,
                batch,
            )
            db.commit()
            batch.clear()
            if processed % 1000000 == 0:
                print(f"Inserted {processed} rows...")

    if batch:
        cur.executemany(
            """
            INSERT INTO variants (chrom, pos, ref, alt, genotype, rsid, source_file, quality)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?);
            """,
            batch,
        )
        db.commit()

    print(f"Done. Total rows inserted: {processed}")
    db.close()


# -----------------------
# CLI
# -----------------------

def main():
    parser = argparse.ArgumentParser(description="Simple genome/FASTA indexer into SQLite.")
    parser.add_argument("input", help="Path to VCF, 23andMe raw, or FASTA file")
    parser.add_argument("--db", default="genome.db", help="Output SQLite DB (default: genome.db)")
    args = parser.parse_args()

    index_file(Path(args.input), Path(args.db))


if __name__ == "__main__":
    main()
