import subprocess
from pathlib import Path

# config
mafft_path = Path("mafft.bat")  # same folder; change if it's elsewhere
input_fasta = Path("ebony_training_plus_ref.fasta")
output_fasta = Path("ebony_training_plus_ref_aligned.fasta")

if not mafft_path.exists():
    raise SystemExit(f"Cannot find {mafft_path}. Put mafft.bat in this folder or update mafft_path.")

if not input_fasta.exists():
    raise SystemExit(f"Cannot find {input_fasta}.")

print(f"Running MAFFT on {input_fasta} ...")

with open(output_fasta, "w") as out:
    # run: mafft.bat --auto ebony_training_plus_ref.fasta
    subprocess.run(
        [str(mafft_path), "--auto", str(input_fasta)],
        stdout=out,
        check=True,
        shell=True,  # needed for .bat on Windows
    )

print(f"Aligned FASTA written to {output_fasta}")
