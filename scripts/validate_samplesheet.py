#!/usr/bin/env python3
"""
fastq2tracks v3.0.4 — Sample sheet validator (patched: tech-replicate aware)
Usage: python3 scripts/validate_samplesheet.py [--check-files] samplesheet.csv

Changes vs v3.0.2:
  - Duplicate sample_id is ALLOWED when tech_replicate values differ.
    The enforced unique key is (sample_id, replicate, tech_replicate).
  - Within a tech-replicate group (same sample_id + replicate) all metadata
    columns except fastq_1, fastq_2, and tech_replicate must be identical.
  - True duplicate rows (all fields identical) still fail.
  - Duplicate composite keys still fail.
  - control_id cross-reference uses the unique_sids set (unchanged semantics).
"""
import csv, sys, os, argparse

REQUIRED_COLS = [
    "sample_id","fastq_1","fastq_2","layout","genome",
    "assay","factor","condition","treatment","cell_type",
    "replicate","tech_replicate","is_control","control_id",
    "macs2_mode","blacklist","chipqc_annotation","output_prefix"
]

CONSISTENT_COLS = [
    "layout","genome","assay","factor","condition","treatment",
    "cell_type","is_control","control_id","macs2_mode","blacklist",
    "chipqc_annotation","output_prefix"
]

VALID_LAYOUTS = {"PE","SE"}
VALID_GENOMES = {"hg38","mm39"}
VALID_MACS2   = {"narrow","broad","both","none"}
VALID_ISCTR   = {"TRUE","FALSE","true","false","1","0","yes","no"}

def err(row, msg):  print(f"  [ERROR] row {row}: {msg}", file=sys.stderr); return 1
def warn(row, msg): print(f"  [WARN]  row {row}: {msg}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("samplesheet")
    ap.add_argument("--check-files", action="store_true",
                    help="Verify FASTQ/BED/RDS files exist on disk")
    args = ap.parse_args()
    errors = 0
    print(f"Validating: {args.samplesheet}")

    with open(args.samplesheet, newline="") as fh:
        reader = csv.DictReader(fh)
        missing_cols = [c for c in REQUIRED_COLS if c not in (reader.fieldnames or [])]
        if missing_cols:
            print(f"[FATAL] Missing columns: {missing_cols}", file=sys.stderr)
            sys.exit(1)
        rows = list(reader)

    composite_keys       = {}   # (sid, rep, tech_rep) -> row-number
    seen_raw_rows        = {}   # frozenset(row.items()) -> row-number
    unique_sids          = set()
    control_ids_declared = set()
    ip_to_control        = {}   # sid -> ctrl_id (first occurrence)
    group_meta           = {}   # (sid, rep) -> {col: val, "first_row": N}

    for i, row in enumerate(rows, start=2):
        sid      = row["sample_id"].strip()
        rep      = row["replicate"].strip()
        tech_rep = row["tech_replicate"].strip()

        # True duplicate row check
        raw_key = frozenset((k, v.strip()) for k, v in row.items())
        if raw_key in seen_raw_rows:
            errors += err(i, f"Exact duplicate row (identical to row {seen_raw_rows[raw_key]})")
        else:
            seen_raw_rows[raw_key] = i

        # Composite key uniqueness (replaces bare sample_id check)
        comp = (sid, rep, tech_rep)
        if comp in composite_keys:
            errors += err(i, f"Duplicate composite key (sample_id='{sid}', "
                             f"replicate='{rep}', tech_replicate='{tech_rep}') "
                             f"— first seen at row {composite_keys[comp]}")
        else:
            composite_keys[comp] = i

        unique_sids.add(sid)

        layout = row["layout"].strip().upper()
        if layout not in VALID_LAYOUTS:
            errors += err(i, f"Invalid layout '{layout}' (PE or SE)")
        if not row["fastq_1"].strip():
            errors += err(i, "fastq_1 is empty")
        if layout == "PE" and not row["fastq_2"].strip():
            errors += err(i, "PE sample missing fastq_2")
        if layout == "SE" and row["fastq_2"].strip():
            warn(i, "SE sample has fastq_2 — will be ignored")
        if row["genome"].strip() not in VALID_GENOMES:
            errors += err(i, f"Invalid genome '{row['genome']}'")
        if row["macs2_mode"].strip().lower() not in VALID_MACS2:
            errors += err(i, f"Invalid macs2_mode '{row['macs2_mode']}'")

        is_ctr = row["is_control"].strip()
        if is_ctr not in VALID_ISCTR:
            errors += err(i, f"Invalid is_control '{is_ctr}'")
        is_ctr_bool = is_ctr.lower() in {"true","1","yes"}
        if is_ctr_bool:
            control_ids_declared.add(sid)
            if row["control_id"].strip():
                warn(i, "Control sample should have empty control_id")
            if row["macs2_mode"].strip().lower() != "none":
                warn(i, "Control sample: macs2_mode should be 'none'")
        else:
            ctrl = row["control_id"].strip()
            if not ctrl:
                warn(i, "IP sample has no control_id — MACS2 will run without control")
            else:
                ip_to_control[sid] = ctrl

        try:
            int(rep)
        except ValueError:
            errors += err(i, f"replicate must be integer, got '{rep}'")

        try:
            int(tech_rep)
        except ValueError:
            errors += err(i, f"tech_replicate must be integer, got '{tech_rep}'")

        if args.check_files:
            for col in ["fastq_1","fastq_2","blacklist","chipqc_annotation"]:
                p = row[col].strip()
                if p and not os.path.exists(p):
                    errors += err(i, f"File not found: {col}={p}")

        # Metadata consistency within tech-replicate group
        grp = (sid, rep)
        meta_snap = {c: row[c].strip() for c in CONSISTENT_COLS}
        if grp not in group_meta:
            group_meta[grp] = {"first_row": i, **meta_snap}
        else:
            ref = group_meta[grp]
            mismatches = [c for c in CONSISTENT_COLS if meta_snap[c] != ref[c]]
            if mismatches:
                errors += err(i,
                    f"Inconsistent metadata within tech-replicate group "
                    f"(sample_id='{sid}', replicate='{rep}'): "
                    f"columns {mismatches} differ from row {ref['first_row']}")

    for sid, ctrl in ip_to_control.items():
        if ctrl not in unique_sids:
            errors += err("?", f"Sample '{sid}' references control_id '{ctrl}' "
                               f"not found in samplesheet")

    n_ip  = sum(1 for r in rows if r["is_control"].strip().lower() not in ("true","1","yes"))
    n_grp = len(group_meta)
    print(f"\nRows parsed          : {len(rows)}")
    print(f"Unique sample groups : {n_grp}  (sample_id + replicate)")
    print(f"IP sample rows       : {n_ip}")
    print(f"Control sample IDs   : {len(control_ids_declared)}")

    if errors:
        print(f"\n[FAIL] {errors} error(s) — fix before running workflow.",
              file=sys.stderr)
        sys.exit(1)
    print("\n[OK] Sample sheet is valid.")


if __name__ == "__main__":
    main()
