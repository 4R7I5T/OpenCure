import os
import subprocess
import sys


def run_step(cmd, description):
    print(f"\n=== {description} ===")
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    python_executable = sys.executable

    # Step 1: Build prompt from annotated VCF
    run_step([python_executable, "Step1.py"], "Step1: Build prompt from VCF")

    # Step 3: Call OpenAI to generate CRISPR edit CSV
    run_step([python_executable, "Step3.py"], "Step3: Generate CRISPR edit CSV with OpenAI")

    # Step 4a: Fetch sequences and run FlashFry
    run_step([python_executable, "Step4a.py"], "Step4a: Fetch sequences and run FlashFry")

    # Step 4b: Run Augustus gene prediction
    run_step(["bash", "Step4b.sh"], "Step4b: Run Augustus gene prediction")

    # Step 5: Run CROPSR on sequences and predictions
    run_step(["bash", "step5.sh"], "Step5: Run CROPSR on sequences and predictions")

    # Step 9: Annotate guides with Cas-OFFinder off-target risk
    run_step([python_executable, "Step9.py"], "Step9: Annotate guides with off-target risk")

    # Step 6: Sample top guides per chromosome
    run_step([python_executable, "Step6.py"], "Step6: Sample top guides per chromosome")


if __name__ == "__main__":
    main()
