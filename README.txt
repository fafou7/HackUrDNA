# ebony-pigmentation-classifier

Small toy project inspired by Pool & Aquadro (2007) to classify *Drosophila melanogaster* **ebony (e)** haplotypes as
more "dark" or "light" abdominal pigmentation based on natural variation.

## Idea

- Use 21 ebony-region sequences (≈20 kb each) from Uganda (dark abdomen) and Kenya (light abdomen).
- Align them, plus any new query sequence (e.g. reference genome region).
- Identify positions where allele frequencies differ between dark vs light groups.
- For a new haplotype, compute a **darkness score** between 0 (light) and 1 (dark) based on matches at those informative sites.

## Files

- `ebony_classifier.py` – core model:
  - `train` command: build JSON model from aligned FASTA + labels
  - `score` command: score a single aligned sequence
- `score_ref_from_alignment.py` – convenience script:
  - reads the last sequence from an aligned multi-FASTA
  - scores it with a given model
- `extract_region.py` – helper to extract the ebony window from `dmel-all-chromosome-r6.65.fasta`
- `labels.csv` – accession → phenotype (`dark` / `light`) for EF114370–EF114390.

## Data

Training sequences are the ebony-region haplotypes from:

> Pool JE & Aquadro CF (2007)  
> *The genetic basis of adaptive pigmentation variation in Drosophila melanogaster.*  
> Molecular Ecology 16(14):2844–2851.

GenBank accessions: **EF114370–EF114390**.

You can download them from NCBI (`Send to → File → FASTA`) into `ebony_training.fasta`.

## Usage

### 1. Align training + reference

Combine training sequences + your reference window:

```bash
cat ebony_training.fasta ref_ebony.fasta > ebony_training_plus_ref.fasta
