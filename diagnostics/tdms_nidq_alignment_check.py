"""
TDMS ↔ NIDQ alignment diagnostic
==================================
Compares valve-open TTL timestamps extracted from the LabVIEW TDMS file
against the same events extracted from the SpikeGLX NIDQ file.

Both recordings capture the same physical TTL wire:
  - TDMS  : behavioural NI-DAQ, usually AI0 (0-5 V, 1 kHz)
  - NIDQ  : neuropixels auxiliary board, digital word 0

If the two clocks drift, TTL timestamps will diverge over the recording.
A systematic offset or a ramp in the difference plot indicates a problem.

Usage
-----
Set the paths in CONFIG and run:
    python tdms_nidq_alignment_check.py

Outputs (saved to SAVE_PATH):
  alignment_<session_id>.png  — per-session diagnostic figure
  alignment_summary.txt       — offset / drift table across all sessions
"""

import os, sys, glob, pickle
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from nptdms import TdmsFile

sys.path.append(r"C:\Users\xaand\Documents\PhD\Analysis\Informatics\Python\Electrophysiology\helpers")
import readSGLX as sglx

# ─── CONFIG ───────────────────────────────────────────────────────────────────
# Paths on Z: drive (mount before running)
NIDQ_ROOT  = r"Z:\raw_data\neuropixels"   # each session subfolder contains *.nidq.bin
TDMS_ROOT  = r"Z:\raw_data\TDMS"          # each session subfolder contains *.tdms

SAVE_PATH  = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\Analysis\tdms_nidq_alignment"

# Which session IDs to check — leave empty [] to process all found subfolders
SESSION_IDS = []

# TDMS channel that carries the valve/trial TTL
TDMS_AI_CHANNEL = "AI0"
TDMS_SR         = 1000.0   # Hz (wf_increment = 0.001 s)

# NIDQ: digital word index and line index that carry the trial TTL
NIDQ_DW         = 0        # digital word 0
NIDQ_LINE       = 0        # line 0 (adjust if needed)

# TTL detection threshold (fraction of signal range)
TTL_THRESHOLD   = 0.5

# For OFC sessions: ep_onsets = TTL_TS[ses][1::4] (valve openings only)
# Set the same stride here so we compare the same subset of TTLs
NIDQ_TTL_STRIDE  = 4       # every 4th event
NIDQ_TTL_OFFSET  = 1       # start index within stride
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(SAVE_PATH, exist_ok=True)


def rising_edges(signal, sr, threshold=TTL_THRESHOLD):
    """Return rising-edge timestamps (seconds) from a 1-D TTL signal."""
    lo, hi = signal.min(), signal.max()
    mid = lo + (hi - lo) * threshold
    above = (signal > mid).astype(np.int8)
    edges = np.where(np.diff(above) == 1)[0] + 1
    return edges / sr


def load_tdms_ttl(tdms_path):
    """Extract rising-edge timestamps from TDMS AI channel."""
    tf  = TdmsFile.read(tdms_path)
    sig = tf["Measured Data"][TDMS_AI_CHANNEL].data
    return rising_edges(sig, TDMS_SR)


def load_nidq_ttl(nidq_bin_path):
    """Extract rising-edge timestamps from neuropixels NIDQ digital channel."""
    meta   = sglx.readMeta(Path(nidq_bin_path))
    raw    = sglx.makeMemMapRaw(Path(nidq_bin_path), meta)
    nsr    = sglx.SampRate(meta)
    dig    = sglx.ExtractDigital(raw, 0, raw.shape[1] - 1,
                                  NIDQ_DW, [NIDQ_LINE], meta)
    sig    = dig[0].astype(float)
    return rising_edges(sig, nsr)


def find_sessions():
    if SESSION_IDS:
        return SESSION_IDS
    subs = [d for d in os.listdir(TDMS_ROOT)
            if os.path.isdir(os.path.join(TDMS_ROOT, d))]
    return sorted(subs)


def check_session(ses_id):
    tdms_dir = os.path.join(TDMS_ROOT, ses_id)
    nidq_dir = os.path.join(NIDQ_ROOT, ses_id)

    tdms_files = glob.glob(os.path.join(tdms_dir, "*.tdms"))
    nidq_files = glob.glob(os.path.join(nidq_dir, "*.nidq.bin"))

    if not tdms_files:
        print(f"  [SKIP] no .tdms found in {tdms_dir}")
        return None
    if not nidq_files:
        print(f"  [SKIP] no .nidq.bin found in {nidq_dir}")
        return None

    tdms_path = tdms_files[0]
    nidq_path = nidq_files[0]

    print(f"  TDMS : {os.path.basename(tdms_path)}")
    print(f"  NIDQ : {os.path.basename(nidq_path)}")

    tdms_ts = load_tdms_ttl(tdms_path)
    nidq_ts = load_nidq_ttl(nidq_path)

    # Select valve-open subset from NIDQ (matches ep_preprocess stride)
    nidq_ts_ep = nidq_ts[NIDQ_TTL_OFFSET::NIDQ_TTL_STRIDE]

    n_tdms = len(tdms_ts)
    n_nidq = len(nidq_ts_ep)
    n_match = min(n_tdms, n_nidq)
    print(f"  TDMS events: {n_tdms}  |  NIDQ valve events: {n_nidq}  |  comparing {n_match}")

    if n_match < 2:
        print("  [SKIP] too few matched events")
        return None

    tdms_m = tdms_ts[:n_match]
    nidq_m = nidq_ts_ep[:n_match]
    diffs  = tdms_m - nidq_m   # TDMS minus NIDQ, in seconds

    median_offset = np.median(diffs)
    drift_ms_per_min = (diffs[-1] - diffs[0]) / (nidq_m[-1] / 60) * 1000

    print(f"  Median TDMS–NIDQ offset : {median_offset*1000:+.1f} ms")
    print(f"  Drift across session    : {drift_ms_per_min:+.2f} ms/min")

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle(f'{ses_id}  —  TDMS ↔ NIDQ alignment', fontsize=11)

    t_min = nidq_m / 60

    # Panel 0: raw TTL onset times side by side
    axes[0].scatter(t_min, tdms_m, s=4, color='steelblue', label='TDMS', alpha=0.7)
    axes[0].scatter(t_min, nidq_m, s=4, color='tomato',    label='NIDQ', alpha=0.7)
    axes[0].set_ylabel('Event time [s]', fontsize=9)
    axes[0].legend(fontsize=8, loc='upper left')
    axes[0].tick_params(direction='in', labelsize=8)

    # Panel 1: difference (TDMS − NIDQ) over time — should be flat if no drift
    axes[1].axhline(0, color='gray', lw=0.8, ls='--')
    axes[1].axhline(median_offset, color='k', lw=0.8, ls=':',
                    label=f'median {median_offset*1000:+.1f} ms')
    axes[1].scatter(t_min, diffs * 1000, s=4, color='darkorange', alpha=0.8)
    axes[1].set_ylabel('TDMS − NIDQ [ms]', fontsize=9)
    axes[1].set_xlabel('Time in session [min]', fontsize=9)
    axes[1].legend(fontsize=8)
    axes[1].tick_params(direction='in', labelsize=8)

    # Annotate drift
    axes[1].text(0.98, 0.05,
                 f'drift: {drift_ms_per_min:+.2f} ms/min',
                 transform=axes[1].transAxes, ha='right', va='bottom', fontsize=8,
                 color='firebrick' if abs(drift_ms_per_min) > 1 else 'k')

    fig.tight_layout()
    fname = os.path.join(SAVE_PATH, f'alignment_{ses_id}.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    print(f"  Saved: {fname}")

    return {
        'ses_id':        ses_id,
        'n_events':      n_match,
        'median_offset': median_offset,
        'drift_ms_per_min': drift_ms_per_min,
        'max_abs_diff':  np.max(np.abs(diffs)),
    }


# ─── Main ─────────────────────────────────────────────────────────────────────
sessions = find_sessions()
results  = []

for ses_id in sessions:
    print(f"\n{'='*60}")
    print(f"Session: {ses_id}")
    r = check_session(ses_id)
    if r:
        results.append(r)

if results:
    summary_path = os.path.join(SAVE_PATH, 'alignment_summary.txt')
    with open(summary_path, 'w') as fout:
        header = f"{'session':<25}  {'n_events':>9}  {'offset_ms':>10}  {'drift_ms/min':>13}  {'max_abs_ms':>10}"
        fout.write(header + '\n')
        fout.write('-' * len(header) + '\n')
        for r in results:
            line = (f"{r['ses_id']:<25}  {r['n_events']:>9}  "
                    f"{r['median_offset']*1000:>+10.2f}  "
                    f"{r['drift_ms_per_min']:>+13.3f}  "
                    f"{r['max_abs_diff']*1000:>10.2f}")
            fout.write(line + '\n')
            print(line)
    print(f"\nSummary saved to {summary_path}")
