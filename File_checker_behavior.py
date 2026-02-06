# -*- coding: utf-8 -*-
"""
This short script checks all the files in raw_data of a list of mice to check their integrity

@author: Xander
"""
import os
import re
#from collections import defaultdict

# Root paths
ROOTS = {
    "experiment_file": r"Z:\raw_data\experiment_file",
    "FLIR":            r"Z:\raw_data\FLIR",
    "TDMS":            r"Z:\raw_data\TDMS",
    "pupil_cam":       r"Z:\raw_data\pupil_cam",
}

# Session IDs
sess_ids = ["250814_XT001","250814_XT002","250814_XT003","250814_XT004","250814_XT005",
            "250815_XT001","250815_XT002","250815_XT003","250815_XT004","250815_XT005",
            "250816_XT001","250816_XT002","250816_XT003","250816_XT004","250816_XT005",
            "250817_XT001","250817_XT002","250817_XT003","250817_XT004","250817_XT005",
            "250818_XT001","250818_XT002","250818_XT003","250818_XT004","250818_XT005",
            "250819_XT001","250819_XT002","250819_XT003","250819_XT004","250819_XT005",
            "250820_XT001","250820_XT002","250820_XT003","250820_XT004","250820_XT005",
            "250821_XT001","250821_XT002","250821_XT003","250821_XT004","250821_XT005",
            "250822_XT001","250822_XT002","250822_XT003","250822_XT004","250822_XT005",
            "250823_XT001","250823_XT002","250823_XT003","250823_XT004","250823_XT005"
            ]  

sess_ids = ["250929_XT006","250929_XT007","250929_XT008"]
sess_ids = ["250930_XT006","250930_XT007","250930_XT008"]
sess_ids = ["251001_XT006","251001_XT007","251001_XT008"]
sess_ids = ["251002_XT006","251002_XT007","251002_XT008"]
sess_ids = ["251003_XT006","251003_XT007","251003_XT008"]
sess_ids = ["251004_XT006","251004_XT007","251004_XT008"]
sess_ids = ["251005_XT006","251005_XT007","251005_XT008"]
sess_ids = ["251006_XT006","251006_XT007","251006_XT008"]
sess_ids = ["251007_XT006","251007_XT007","251007_XT008"]
#sess_ids = ["251008_XT006","251008_XT007","251008_XT008"]
#sess_ids = ["251009_XT006","251009_XT007","251009_XT008"]
sess_ids = ["251003_XT006","251003_XT007","251003_XT008","250817_XT001","250817_XT002","250817_XT003","250817_XT004","250817_XT005"]

 
#The above are all done and OK

# Your session IDs

expected_flir_seq = 130  # <-- set this (optional; used for an extra check)

# -------------------- helpers --------------------

def get_expected_seq_count(sess_id):
    if isinstance(expected_flir_seq, dict):
        return expected_flir_seq.get(sess_id, None)
    return expected_flir_seq

def list_dir(path):
    try:
        return os.listdir(path)
    except FileNotFoundError:
        return None

def status_dict():
    return {"ok": True, "warns": [], "fails": []}

def merge_status(base, new):
    base["ok"] = base["ok"] and new["ok"]
    base["warns"].extend(new["warns"])
    base["fails"].extend(new["fails"])
    return base

def print_and_record(level, sess_id, msg, stash):
    tag = {"ok":"[OK]  ", "warn":"[WARN]", "fail":"[FAIL]"}[level]
    print(f"{tag}: {msg}")
    if level == "warn":
        stash["warns"].append(msg)
    elif level == "fail":
        stash["ok"] = False
        stash["fails"].append(msg)

# -------------------- checks --------------------

def find_exp_pair(exp_dir, sess_id):
    """
    In experiment_file/<sess_id>:
      - exactly one .txt and one .mat
      - both filenames end with _<suffix> (suffix = part after '_' in sess_id)
      - both share the same base (e.g., Day4_XT005.txt & Day4_XT005.mat)
    Returns (status, exp_prefix or None)
    """
    st = status_dict()
    suffix = sess_id.split("_", 1)[1] if "_" in sess_id else sess_id

    files = list_dir(exp_dir)
    if files is None:
        print_and_record("fail", sess_id, f"experiment_file folder missing: {exp_dir}", st)
        return st, None

    pat = re.compile(rf".*_{re.escape(suffix)}\.(txt|mat)$", re.IGNORECASE)
    matches = [f for f in files if pat.fullmatch(f)]
    txts = [f for f in matches if f.lower().endswith(".txt")]
    mats = [f for f in matches if f.lower().endswith(".mat")]

    if len(txts) == 1 and len(mats) == 1:
        base_txt = os.path.splitext(txts[0])[0]
        base_mat = os.path.splitext(mats[0])[0]
        if base_txt == base_mat:
            print_and_record("ok", sess_id, f"experiment_file: matching pair found: {txts[0]}, {mats[0]}", st)
            return st, base_txt  # exp_prefix (e.g., Day4_XT005)
        else:
            print_and_record("fail", sess_id, f"experiment_file: base names differ: {txts[0]} vs {mats[0]}", st)
            return st, None
    else:
        print_and_record("fail", sess_id, f"experiment_file: need exactly one .txt and one .mat ending with _{suffix}; found: {matches}", st)
        return st, None

def check_flir(flir_dir, sess_id):
    """
    In FLIR/<sess_id>:
      - exclusively .seq files
      - count matches expected (if provided)
    Returns (status, n_trials)
    """
    st = status_dict()
    files = list_dir(flir_dir)
    if files is None:
        print_and_record("fail", sess_id, f"FLIR folder missing: {flir_dir}", st)
        return st, None

    non_dirs = [f for f in files if os.path.isfile(os.path.join(flir_dir, f))]
    seqs = [f for f in non_dirs if f.lower().endswith(".seq")]
    others = [f for f in non_dirs if not f.lower().endswith(".seq")]

    n_trials = len(seqs)
    expected = get_expected_seq_count(sess_id)

    if expected is None:
        print_and_record("warn", sess_id, f"FLIR: no expected .seq count provided; found {n_trials} .seq (n_trials set)", st)
    else:
        if n_trials != expected:
            print_and_record("warn", sess_id, f"FLIR: found {n_trials} .seq, expected {expected}", st)

    if others:
        print_and_record("warn", sess_id, f"FLIR: non-.seq present: {others}", st)

    if st["ok"]:
        print_and_record("ok", sess_id, f"FLIR: exactly {n_trials} .seq files and nothing else (n_trials={n_trials})", st)

    return st, n_trials

def check_tdms(tdms_dir, sess_id, exp_prefix):
    """
    In TDMS/<sess_id>:
      - total files >= 4
      - includes ≥1 .tdms (or .tdsm)
      - includes ≥1 .tdms_index
      - includes ≥2 .txt starting with f"{exp_prefix}_"
      - any files beyond the minimal core are reported as warnings (extras)
    """
    st = status_dict()
    files = list_dir(tdms_dir)
    if files is None:
        print_and_record("fail", sess_id, f"TDMS folder missing: {tdms_dir}", st)
        return st

    non_dirs = [f for f in files if os.path.isfile(os.path.join(tdms_dir, f))]
    total = len(non_dirs)

    tdms_like = [f for f in non_dirs if os.path.splitext(f)[1].lower() in (".tdms", ".tdsm")]
    tdms_index = [f for f in non_dirs if f.lower().endswith(".tdms_index")]
    txts = [f for f in non_dirs if f.lower().endswith(".txt")]
    txts_good = [f for f in txts if f.startswith(exp_prefix + "_")]

    # Build minimal required set to identify extras
    minimal = set()
    if tdms_like:  minimal.add(tdms_like[0])
    if tdms_index: minimal.add(tdms_index[0])
    minimal.update(txts_good[:2])
    extras = [f for f in non_dirs if f not in minimal] if total > 4 else []

    if total < 4:
        print_and_record("fail", sess_id, f"TDMS: only {total} files (need ≥4)", st)
    if not tdms_like:
        print_and_record("fail", sess_id, "TDMS: missing .tdms file", st)
    if not tdms_index:
        print_and_record("fail", sess_id, "TDMS: missing .tdms_index file", st)
    if len(txts_good) < 2:
        print_and_record("fail", sess_id, f"TDMS: {len(txts_good)} matching .txt (need ≥2 starting with '{exp_prefix}_')", st)

    extra_txts_wrong_prefix = [f for f in txts if f not in txts_good]
    if extra_txts_wrong_prefix:
        print_and_record("warn", sess_id, f"TDMS: TXT not matching prefix '{exp_prefix}_': {extra_txts_wrong_prefix}", st)

    if extras:
        print_and_record("warn", sess_id, f"TDMS: extra/error files detected: {extras}", st)

    if st["ok"]:
        print_and_record("ok", sess_id, "TDMS: core files present (≥4, tdms(+index), 2 prefixed .txt).", st)
    return st

def check_pupil_cam(pupil_dir, sess_id, n_trials):
    """
    In pupil_cam/<sess_id>:
      - exactly n_trials .tif files (where n_trials = number of FLIR .seq files)
      - exactly two subfolders total:
            * one named '.camlog' (case-insensitive)
            * one additional folder (name arbitrary)
      - warn if any non-.tif files are present
    """
    st = status_dict()
    if n_trials is None:
        print_and_record("fail", sess_id, "pupil_cam check skipped: n_trials unknown (FLIR missing).", st)
        return st

    files = list_dir(pupil_dir)
    if files is None:
        print_and_record("fail", sess_id, f"pupil_cam folder missing: {pupil_dir}", st)
        return st

    non_dirs = [f for f in files if os.path.isfile(os.path.join(pupil_dir, f))]
    dirs     = [f for f in files if os.path.isdir(os.path.join(pupil_dir, f))]

    tif_files = [f for f in non_dirs if f.lower().endswith(".tif")]
    other_files = [f for f in non_dirs if not f.lower().endswith(".tif")]

    # .camlog folder + exactly one more folder
    dir_lower = [d.lower() for d in dirs]
    camlog_present = ".camlog" in dir_lower
    if not camlog_present:
        print_and_record("warn", sess_id, "pupil_cam: missing '.camlog' folder", st)

    if len(dirs) != 2:
        print_and_record("warn", sess_id, f"pupil_cam: expected exactly 2 subfolders ('.camlog' + 1 other), found {len(dirs)}: {dirs}", st)

    # tif count must match n_trials
    if not(len(tif_files) == n_trials or len(tif_files) == n_trials+1) :
        print_and_record("fail", sess_id, f"pupil_cam: found {len(tif_files)} .tif, expected n_trials={n_trials}", st)

    if other_files:
        print_and_record("warn", sess_id, f"pupil_cam: non-.tif files present: {other_files}", st)

    if st["ok"]:
        print_and_record("ok", sess_id, f"pupil_cam: {len(tif_files)} .tif files match FLIR files and required folders present.", st)

    return st

# -------------------- main --------------------

def main():
    session_summaries = {}

    for sess_id in sess_ids:
        print("\n" + "="*72)
        print(f"Checking {sess_id}")
        print("="*72)

        exp_dir     = os.path.join(ROOTS["experiment_file"], sess_id)
        flir_dir    = os.path.join(ROOTS["FLIR"], sess_id)
        tdms_dir    = os.path.join(ROOTS["TDMS"], sess_id)
        pupil_dir   = os.path.join(ROOTS["pupil_cam"], sess_id)

        overall = status_dict()

        # 1) experiment_file
        st_exp, exp_prefix = find_exp_pair(exp_dir, sess_id)
        merge_status(overall, st_exp)

        # 2) FLIR
        st_flir, n_trials = check_flir(flir_dir, sess_id)
        merge_status(overall, st_flir)

        # 3) TDMS (only if we have exp_prefix to validate TXT naming)
        if exp_prefix:
            st_tdms = check_tdms(tdms_dir, sess_id, exp_prefix)
            merge_status(overall, st_tdms)
        else:
            print_and_record("fail", sess_id, "TDMS check skipped: missing experiment_file prefix (fix experiment_file first).", overall)

        # 4) pupil_cam (uses n_trials from FLIR)
        st_pupil = check_pupil_cam(pupil_dir, sess_id, n_trials)
        merge_status(overall, st_pupil)

        # --- Final status per session ---
        if overall["fails"]:
            final = "FAIL"
        elif overall["warns"]:
            final = "OK with warnings"
        else:
            final = "OK"

        print("-" * 72)
        print(f"FINAL STATUS [{sess_id}]: {final}")
        print(f"  n_trials (from FLIR): {n_trials if n_trials is not None else 'unknown'}")
        if overall["warns"]:
            print(f"  Warnings ({len(overall['warns'])}):")
            for w in overall["warns"]:
                print(f"   - {w}")
        if overall["fails"]:
            print(f"  Failures ({len(overall['fails'])}):")
            for f in overall["fails"]:
                print(f"   - {f}")

        session_summaries[sess_id] = {"final": final, "n_trials": n_trials}

    # --- Overall summary across sessions ---
    print("\n" + "="*72)
    print("OVERALL SUMMARY")
    print("="*72)
    counts = {"OK":0, "OK with warnings":0, "FAIL":0}
    for sess_id, info in session_summaries.items():
        final = info["final"]
        counts[final] = counts.get(final, 0) + 1
        print(f"{sess_id}: {final} (n_trials={info['n_trials']})")

if __name__ == "__main__":
    main()
