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

## FASTA files in this repo

This repo includes several small FASTA files as example inputs and intermediates.  
They are all subsets or alignments around the **ebony** region from *Drosophila melanogaster*.

### Core example data

- **`ebony_training.fasta`**  
  Raw, unaligned ebony-region sequences from **21 natural lines** used as training data.  
  These are the EF114370–EF114390 sequences from:  
  Pool & Aquadro (2007) *The genetic basis of adaptive pigmentation variation in Drosophila melanogaster*.

- **`ref_ebony.fasta`**  
  Example ebony-region sequence extracted from the **reference D. melanogaster genome**  
  (chromosome 3R, ~8.39–8.42 Mb).  
  This is treated as a “query” sequence when we ask,  
  *“Does the reference look more like the dark or light ebony haplotypes?”*

### Alignments

- **`ebony_training_aligned.fasta`**  
  MAFFT alignment of the 21 training sequences only  
  (Uganda = dark, Kenya = light).  
  Useful if you just want to re-run or inspect the alignment of the training set.

- **`ebony_training_plus_ref.fasta`**  
  Unaligned FASTA containing:
  - the 21 training sequences, plus  
  - the reference ebony sequence from `ref_ebony.fasta`.  
  This is the input that was sent to MAFFT to produce the joint alignment below.

- **`ebony_training_plus_ref_aligned.fasta`**  
  MAFFT alignment of the 21 training sequences **plus** the reference ebony sequence.  
  This is the main aligned file used by `ebony_classifier.py train` to:
  - identify informative positions that differ between dark and light lines, and  
  - later score the reference sequence.

- **`ebony_plus_ref_aligned.fasta`**  
  Earlier/alternate alignment of training + reference sequences kept as an intermediate.  
  Not required for the basic pipeline, but left here as a reference/example.

### Single-sequence alignment snippets

- **`K60_aligned.fasta`**, **`U70_aligned.fasta`**  
  Individual aligned sequences (extracted from the multi-FASTA) for specific lines  
  used during testing and debugging. They can be used with the `score` command as  
  minimal examples of how to score one aligned sequence at a time.

You should get something like:

K60 → darkness_score close to 0
(more matches to the “light” alleles, consistent with provided gene)

U70 → darkness_score close to 1
(more matches to the “dark” alleles, consistent with provided gene)

