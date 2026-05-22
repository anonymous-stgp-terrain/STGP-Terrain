import pandas as pd
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

def update_final_results(method, table_id, rmse=np.nan, nlpd=np.nan, runtime=np.nan,
                         version="x", final_csv=None):
    """
    version : "x"  -> writes to RMSE(x)   column (Table 4 only)
              "xs" -> writes to RMSE(x+s)  column (Table 4 only)
              ignored for Table 2 and Table 3
    """
    if final_csv is None:
        final_csv = REPO_ROOT / "results" / "final results.csv"
    else:
        final_csv = Path(final_csv)

    if not final_csv.exists():
        raise FileNotFoundError(f"{final_csv} does not exist. Please create the template first.")

    df = pd.read_csv(final_csv, header=None, dtype=str).fillna("")

    def find_row(table_label, method_name):
        table_rows = df.index[df[0] == table_label].tolist()
        if not table_rows:
            raise ValueError(f"Could not find {table_label}")
        table_row = table_rows[0]
        next_table_rows = df.index[
            (df[0].isin(["Table 2", "Table 3", "Table 4"])) & (df.index > table_row)
        ].tolist()
        end_row = len(df) - 1 if not next_table_rows else next_table_rows[0] - 1
        method_rows = df.index[df[0] == method_name].tolist()
        method_rows = [r for r in method_rows if table_row < r <= end_row]
        if not method_rows:
            raise ValueError(f"Could not find method '{method_name}' under {table_label}")
        return method_rows[0]

    r = find_row(table_id, method)
    fmt = lambda x: "" if pd.isna(x) else f"{x:.2f}"

    if table_id in ["Table 2", "Table 3"]:
        df.iat[r, 1] = fmt(rmse)
        df.iat[r, 2] = fmt(nlpd)
        df.iat[r, 3] = fmt(runtime)
    elif table_id == "Table 4":
        if version == "x":
            df.iat[r, 1] = fmt(rmse)
        elif version == "xs":
            df.iat[r, 2] = fmt(rmse)
        else:
            raise ValueError("version must be 'x' or 'xs'")
    else:
        raise ValueError("table_id must be one of: Table 2, Table 3, Table 4")

    df.to_csv(final_csv, header=False, index=False)
