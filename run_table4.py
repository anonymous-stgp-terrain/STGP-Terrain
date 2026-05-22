#!/usr/bin/env python3
import glob, os, shutil, subprocess, sys, time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
CODE_DIR  = REPO_ROOT / "code"
OUT_DIR   = REPO_ROOT / "results" / "intermediate"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def find_rscript():
    if shutil.which("Rscript"):
        return ["Rscript"]
    if shutil.which("apptainer"):
        sifs = sorted(glob.glob("/usr/local/pace-apps/manual/packages/r/*/r-*.sif"))
        if sifs:
            print(f"Found R via apptainer: {sifs[-1]}")
            return ["apptainer", "exec", sifs[-1], "Rscript"]
    if sys.platform == "win32":
        for base in [
            Path(os.environ.get("LOCALAPPDATA","")) / "Programs" / "R",
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
print("Running Table 4...")
print(f"Output directory: {OUT_DIR}\n")

def run(cmd, cwd=REPO_ROOT):
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        sys.exit(f"\nFailed: {' '.join(str(c) for c in cmd)}")

def timed(label, out_file, cmd, cwd=REPO_ROOT):
    t0 = time.time()
    run(cmd, cwd=cwd)
    elapsed = int(time.time()-t0)
    (OUT_DIR/out_file).write_text(str(elapsed))
    print(f"{label} runtime (sec): {elapsed}")

script_start = time.time()

timed("STGP",           "table4_stgp_runtime.txt",          [*RSCRIPT, str(CODE_DIR / "Table4(STGP).R")])
timed("NN",             "table4_nn_runtime.txt",            [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", str(CODE_DIR / "Table4(NN).ipynb"),      "--inplace"])
timed("BNN",            "table4_bnn_runtime.txt",           [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", str(CODE_DIR / "Table4(BNN).ipynb"),     "--inplace"])
timed("XGBoost",        "table4_xgboost_runtime.txt",       [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", str(CODE_DIR / "Table4(XGBoost).ipynb"), "--inplace"])
timed("twinGP+Binning", "table4_twinGP_binning_runtime.txt",[*RSCRIPT, str(CODE_DIR / "Table4(twinGP+Binning).R")])
timed("Binning-hetGP",  "table4_hetgp_runtime.txt",         [*RSCRIPT, str(CODE_DIR / "Table4(Binning-hetGP).R")])
timed("Update results", "table4_python_update_runtime.txt", [sys.executable, str(CODE_DIR / "update_final_results.py")])

total = int(time.time()-script_start)
(OUT_DIR/"table4_total_runtime.txt").write_text(str(total))
print(f"\nTotal runtime (sec): {total}")
print("\nTable 4 finished successfully.")
