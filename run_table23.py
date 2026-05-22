#!/usr/bin/env python3
import glob, os, shutil, subprocess, sys, time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
CODE_DIR  = REPO_ROOT / "code"
OUT_DIR   = REPO_ROOT / "results" / "intermediate"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def find_rscript():
    """Return Rscript command as a list, handling PATH, apptainer (PACE), and Windows."""
    # 1. On PATH (macOS, standard Linux)
    if shutil.which("Rscript"):
        return ["Rscript"]

    # 2. PACE / apptainer container
    if shutil.which("apptainer"):
        sifs = sorted(glob.glob("/usr/local/pace-apps/manual/packages/r/*/r-*.sif"))
        if sifs:
            print(f"Found R via apptainer: {sifs[-1]}")
            return ["apptainer", "exec", sifs[-1], "Rscript"]

    # 3. Windows user-level install
    if sys.platform == "win32":
        for base in [
            Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "R",
            Path(os.environ.get("PROGRAMFILES", r"C:\Program Files")) / "R",
        ]:
            if base.is_dir():
                for r_dir in sorted(base.iterdir(), reverse=True):
                    exe = r_dir / "bin" / "Rscript.exe"
                    if exe.exists():
                        print(f"Found Rscript at: {exe}")
                        return [str(exe)]

    sys.exit("ERROR: Rscript not found. Please install R.")

RSCRIPT = find_rscript()
print("Running Tables 2-3...")
print(f"Output directory: {OUT_DIR}\n")

def run(cmd):
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        sys.exit(f"\nFailed: {' '.join(str(c) for c in cmd)}")

def timed(label, out_file, cmd):
    t0 = time.time()
    run(cmd)
    elapsed = int(time.time() - t0)
    (OUT_DIR / out_file).write_text(str(elapsed))
    print(f"{label} runtime (sec): {elapsed}")

script_start = time.time()

timed("STGP",           "stgp_runtime.txt",   [*RSCRIPT, str(CODE_DIR / "Table2-Table3(STGP).R")])
timed("twinGP+Binning", "twingp_runtime.txt",  [*RSCRIPT, str(CODE_DIR / "Table2_Table3(twinGP+Binning).R")])
timed("NN",             "nn_runtime.txt",      [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", str(CODE_DIR / "Table2-Table3(NN).ipynb"), "--inplace"])
timed("BNN",            "bnn_runtime.txt",     [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", str(CODE_DIR / "Table2-Table3(BNN).ipynb"), "--inplace"])
timed("XGBoost",        "xgboost_runtime.txt", [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", str(CODE_DIR / "Table2-Table3(XGBoost).ipynb"), "--inplace"])
timed("Update results", "python_runtime.txt",  [sys.executable, str(CODE_DIR / "update_final_results.py")])

total = int(time.time() - script_start)
(OUT_DIR / "table23_total_runtime.txt").write_text(str(total))
print(f"\nTotal runtime (sec): {total}")
print("\nTables 2-3 finished successfully.")
