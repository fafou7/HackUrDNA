EBONY PIGMENTATION CLASSIFIER – QUICK GUIDE
===========================================

Purpose
- Small demo project: given an aligned ebony gene sequence from Drosophila melanogaster,
  estimate whether it looks more like a “dark” (high-altitude Uganda) or “light”
  (low-altitude Kenya) haplotype.
- Output: a darkness score between 0 (light-like) and 1 (dark-like).
- Training data: ebony haplotypes EF114370–EF114390 from:

  Pool JE & Aquadro CF (2007)
  “The genetic basis of adaptive pigmentation variation in Drosophila melanogaster”
  Molecular Ecology 16(14):2844–2851
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2650379/


1. REQUIREMENTS
----------------
- Python 3 installed
- Git only if you clone from GitHub (not needed if you download ZIP)

All small FASTA files and labels for the demo are included in the repo.
You do NOT need a full Drosophila genome to run the basic test.


2. GET THE CODE
----------------
Option 1 – Using git:
  git clone https://github.com/YOUR_USERNAME/HackUrDNA.git
  cd HackUrDNA

Option 2 – Using ZIP:
  - Download the repo as ZIP from GitHub
  - Unzip it
  - Open a terminal in that folder

Example on Windows (PowerShell):
  cd "C:\Users\YOURNAME\Desktop\HackUrDNA"


3. QUICK DEMO: TRAIN AND SCORE THE REFERENCE SEQUENCE
-----------------------------------------------------

This is the shortest path to see the classifier working end-to-end.

Step 3.1 – Train the model

From inside the repo folder, run:
  python ebony_classifier.py train --fasta ebony_training_plus_ref_aligned.fasta --labels labels.csv --out ebony_model_v2.json

This:
- reads the aligned multi-FASTA “ebony_training_plus_ref_aligned.fasta”
  (21 natural ebony haplotypes + 1 reference sequence),
- uses “labels.csv” to know which sequences are “dark” (Uganda) and “light” (Kenya),
- finds informative positions where dark and light alleles differ,
- writes the model to “ebony_model_v2.json”.

You should see messages like:
  Loaded 11 dark and 10 light sequences.
  Found XX informative positions.
  Model written to ebony_model_v2.json

Step 3.2 – Score the reference ebony haplotype

Run:
  python score_ref_from_alignment.py

This:
- reads “ebony_training_plus_ref_aligned.fasta”,
- takes the last sequence (header: “3R:8390000-8423000” = reference ebony window),
- loads “ebony_model_v2.json”,
- prints a darkness score.

Example output:
  Scoring sequence (from alignment): 3R:8390000-8423000
  Informative sites used : 90
  Matches dark alleles   : 25
  Matches light alleles  : 40
  Darkness score (0=light,1=dark): 0.278

Interpretation:
The reference ebony haplotype looks more like the “light” Kenya lines than the “dark” Uganda lines.
If you see a darkness score printed, the basic setup is working.


4. OPTIONAL: SCORE SPECIFIC TRAINING LINES (K60 AND U70)
--------------------------------------------------------

You can also score individual training haplotypes directly from the alignment.

From the FASTA headers in “ebony_training_plus_ref_aligned.fasta”:
- K60 (Kenya, light) has header containing:
    EF114384.1 Drosophila melanogaster isolate K60 ...
- U70 (Uganda, dark) has header containing:
    EF114375.1 Drosophila melanogaster isolate U70 ...

Internally, IDs are simplified to the first token without “.1”:
- K60 -> EF114384
- U70 -> EF114375

Use “score_from_alignment.py” to look up sequences by these accessions.

Step 4.1 – Score K60 (light line)
  python score_from_alignment.py --sample EF114384

Expected:
Darkness score close to 0 (light-like).

Step 4.2 – Score U70 (dark line)
  python score_from_alignment.py --sample EF114375

Expected:
Darkness score close to 1 (dark-like).

If K60 scores light and U70 scores dark, the classifier is consistent with the training labels.


5. OPTIONAL: REBUILD THE REFERENCE SEQUENCE FROM A FULL GENOME
--------------------------------------------------------------

Only needed if you want to redo the reference ebony extraction from a full genome.

Step 5.1 – Download a Drosophila genome (e.g. from FlyBase)
Example filename:
  dmel-all-chromosome-r6.65.fasta

Place this file in the repo folder.

Step 5.2 – Extract the ebony window
Run:
  python extract_region.py --genome dmel-all-chromosome-r6.65.fasta --chrom 3R --start 8390000 --end 8423000 --out ref_ebony.fasta

Step 5.3 – Rebuild the training + reference FASTA

On Linux/macOS:
  cat ebony_training.fasta ref_ebony.fasta > ebony_training_plus_ref.fasta

On Windows (PowerShell):
  type ebony_training.fasta ref_ebony.fasta > ebony_training_plus_ref.fasta

Step 5.4 – Realign with MAFFT
Align “ebony_training_plus_ref.fasta” with MAFFT (locally or via the MAFFT web server)
and save the result as:
  ebony_training_plus_ref_aligned.fasta

Step 5.5 – Retrain and rescore
Repeat:
  python ebony_classifier.py train --fasta ebony_training_plus_ref_aligned.fasta --labels labels.csv --out ebony_model_v2.json
  python score_ref_from_alignment.py


6. NOTES
--------

- This is a small proof-of-concept, not a production-ready predictor.
- All sequences used for scoring must be in the same multiple-sequence alignment
  (same length and gaps) as the one used for training, otherwise the script will
  refuse to score and report a length mismatch.
