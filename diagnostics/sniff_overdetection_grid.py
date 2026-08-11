"""
Sniff detection — two-sided grid-search sweep
=============================================
Sweeps a full pre-processing parameter grid and scores each combination on BOTH
failure modes of peak-based inhalation detection, with no reliance on a "good"
reference session to define ground truth:

  - OVER-detection  → short_isi_frac : fraction of ISIs < 1/MAX_BR_HZ, i.e.
                      physiologically-impossible intervals (spurious double-counts).
  - UNDER-detection → gap_frac       : a missed sniff leaves an ISI ≈ k× the LOCAL
                      breathing period, so estimated misses per interval =
                      round(ISI / local_median) − 1. This is the mirror image of the
                      short-ISI metric and is computed from the data itself.

The best operating point minimises  cost = gap_frac + OD_WEIGHT × short_isi_frac
— it fills in missed sniffs without introducing impossible ones. The residual
gap floor (gaps that no parameter removes) is the recall ceiling of threshold
detection on this signal; if it is too high the fix is preprocessing, not params.

The grid session ({GRID_SES}) stores, under `ml_grid`, the 390 parameter combos
(SmoothSpan × MinPeakDistance × MinPeakHeightFactor × MinPeakPromFactor), each
with its own per-trial `ml_inh_onsets`. A REFERENCE session ({REF_SES}) is shown
for CONTEXT ONLY (heatmap colour centring / green dots) — it does not set the
recommended combo; optimal params are signal-dependent and need not transfer.

Always uses ml_* fields (ft_* never used for breathing analysis).

Usage
-----
Edit the CONFIG block, then run. Outputs go to SAVE_DIR:
  grid_gap_frac_<stamp>.png       — UNDER-detection (missed-sniff) heatmap
  grid_short_isi_frac_<stamp>.png — OVER-detection heatmap
  grid_mean_rate_<stamp>.png      — detection-rate heatmap
  grid_summary_<stamp>.txt        — best operating point + top-20 by two-sided cost
@author: Xander
"""
#%%
import sys
sys.path.append(r"C:\Users\xaand\Documents\PhD\Analysis\Informatics\Python\Electrophysiology\helpers")
import os
import glob
from datetime import datetime
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ============================================================
# Paths
# ============================================================
MAT_ROOT = r"C:\Users\xaand\Documents\PhD\Experiments\postprocessed_files"
SAVE_DIR = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\Analysis\sniff_overdetection_grid"

# ============================================================
# CONFIGURATION
# ============================================================
GRID_SES  = '250118_KK152_grid'   # session whose sniffdata holds ml_grid / ml_grid_params
REF_SES   = '260512_XT026'        # trusted reference session (detector='ml')

PRE_EVENT     = 4.0           # s before odor onset (assumed if no pre_event in mat)
ANALYSIS_WIN  = (-4.0, 8.0)   # s rel. to odor onset over which sniffs are counted
MAX_BR_HZ     = 14.0          # max physiologically plausible inhalation rate (Hz)
MIN_BR_HZ     = 0.3           # below this an ISI is a long pause, excluded from rate stats
FLAG_RATIO    = 1.25          # combo flagged as over-detecting if metric > FLAG_RATIO × ref

# Gap (under-detection) metric — a missed sniff leaves an ISI ≈ k× the LOCAL breathing
# period. Estimated misses per gap = round(ISI / local_median) − 1.
GAP_LOCAL_K   = 3             # local median window = ±GAP_LOCAL_K neighbouring ISIs
GAP_MIN_RATIO = 1.5           # ISI must exceed this × local median to count as a gap
GAP_MAX_ISI   = 1.0 / MIN_BR_HZ   # ISIs longer than this are real pauses, not misses

# Two-sided operating point: choose the combo that minimises a combined cost of
# under-detection (gaps) and over-detection (short ISIs). Weight balances the two.
OD_WEIGHT     = 1.0           # cost = gap_frac + OD_WEIGHT × short_isi_frac

# parameter axes (used for the heatmap layout)
PARAM_KEYS = ['SmoothSpan', 'MinPeakDistance', 'MinPeakHeightFactor', 'MinPeakPromFactor']
FACET_KEYS = ['SmoothSpan', 'MinPeakDistance']   # → one heatmap panel per combination
AX_X       = 'MinPeakPromFactor'                 # heatmap x-axis
AX_Y       = 'MinPeakHeightFactor'               # heatmap y-axis

# ============================================================
# Metric core (shared logic with sniff_overdetection_check.py)
# ============================================================
WIN_DUR = ANALYSIS_WIN[1] - ANALYSIS_WIN[0]
MIN_ISI = 1.0 / MAX_BR_HZ
MAX_ISI = 1.0 / MIN_BR_HZ


def _trial_onsets_s(inh_frames, sf, pre):
    if inh_frames is None or len(inh_frames) == 0:
        return np.array([])
    s = np.sort(np.asarray(inh_frames, float) / sf - pre)
    return s[(s >= ANALYSIS_WIN[0]) & (s <= ANALYSIS_WIN[1])]


def _local_median(isis, k=GAP_LOCAL_K):
    """Median ISI in a ±k-neighbour window around each ISI (local breathing period)."""
    out = np.empty_like(isis, dtype=float)
    for i in range(len(isis)):
        out[i] = np.median(isis[max(0, i - k):min(len(isis), i + k + 1)])
    return out


def metrics_from_onsets(onsets_list, sf, pre):
    """Per-trial onset arrays → two-sided detection-quality metric dict."""
    counts, rates, all_isis = [], [], []
    n_detected, n_missed = 0, 0          # for the gap (under-detection) metric
    for t in range(len(onsets_list)):
        s = _trial_onsets_s(onsets_list[t], sf, pre)
        counts.append(len(s))
        rates.append(len(s) / WIN_DUR)
        if len(s) >= 2:
            isis = np.diff(s)
            all_isis.append(isis)
        if len(s) >= 2 * GAP_LOCAL_K:    # need enough ISIs for a stable local median
            isis = np.diff(s)
            lm   = _local_median(isis)
            # estimated missed sniffs per interval; only gaps that are not real pauses
            est  = np.maximum(0.0, np.round(isis / lm) - 1.0)
            est[(isis < GAP_MIN_RATIO * lm) | (isis > GAP_MAX_ISI)] = 0.0
            n_detected += len(isis)
            n_missed   += float(est.sum())
    all_isis = np.concatenate(all_isis) if all_isis else np.array([])
    valid    = all_isis[(all_isis >= MIN_ISI) & (all_isis <= MAX_ISI)]
    expected = n_detected + n_missed
    return dict(
        n_trials       = len(onsets_list),
        mean_count     = float(np.mean(counts)),
        mean_rate      = float(np.mean(rates)),
        short_isi_frac = float(np.mean(all_isis < MIN_ISI)) if all_isis.size else np.nan,
        gap_frac       = float(n_missed / expected) if expected else np.nan,
        median_isi     = float(np.median(valid)) if valid.size else np.nan,
        inst_rate_med  = float(np.median(1.0 / valid)) if valid.size else np.nan,
    )


def _load_mat(folder):
    """Load the *_1_sniffdata_center.mat inside a session folder.
    The file is prefixed by the unique_id, which may differ from the folder name
    (e.g. the grid folder '250118_KK152_grid' holds '250118_KK152_1_...')."""
    d = os.path.join(MAT_ROOT, folder)
    hits = glob.glob(os.path.join(d, '*_1_sniffdata_center.mat'))
    if not hits:
        raise FileNotFoundError(f'No *_1_sniffdata_center.mat in {d}')
    return sio.loadmat(hits[0], simplify_cells=True)


# ============================================================
# Reference metrics
# ============================================================
#%%
os.makedirs(SAVE_DIR, exist_ok=True)
RUN_STAMP = datetime.now().strftime('%Y%m%d_%H%M%S')

ref_mat = _load_mat(REF_SES)
ref = metrics_from_onsets(ref_mat['ml_inh_onsets'],
                          int(ref_mat['samp_freq']),
                          float(ref_mat.get('pre_event', PRE_EVENT)))
print(f"Reference {REF_SES}: rate={ref['mean_rate']:.3f} Hz, "
      f"inh/tr={ref['mean_count']:.2f}, shortISI={ref['short_isi_frac']*100:.2f}%")

# ============================================================
# Grid metrics
# ============================================================
grid_mat = _load_mat(GRID_SES)
sf_grid  = int(grid_mat['samp_freq'])
pre_grid = float(grid_mat.get('pre_event', PRE_EVENT))
grid     = grid_mat['ml_grid']           # list of dicts: {'params':..., 'ml_inh_onsets':...}

rows = []
for entry in grid:
    p = entry['params']
    m = metrics_from_onsets(entry['ml_inh_onsets'], sf_grid, pre_grid)
    rows.append({**{k: p[k] for k in PARAM_KEYS}, **m})
print(f"Computed metrics for {len(rows)} grid combinations.")

# axis value sets
facet_vals = {k: sorted({r[k] for r in rows}) for k in FACET_KEYS}
xs = sorted({r[AX_X] for r in rows})
ys = sorted({r[AX_Y] for r in rows})
lookup = {(r[FACET_KEYS[0]], r[FACET_KEYS[1]], r[AX_X], r[AX_Y]): r for r in rows}


def _grid_array(metric, f0, f1):
    """2-D array [y, x] of a metric for one facet (f0=FACET_KEYS[0], f1=FACET_KEYS[1])."""
    a = np.full((len(ys), len(xs)), np.nan)
    for yi, yv in enumerate(ys):
        for xi, xv in enumerate(xs):
            r = lookup.get((f0, f1, xv, yv))
            if r is not None:
                a[yi, xi] = r[metric]
    return a


# ============================================================
# Heatmaps — one figure per metric, faceted by SmoothSpan × MinPeakDistance
# ============================================================
def heatmap_figure(metric, ref_val, title, cmap='RdBu_r', center_on_ref=True):
    f0s, f1s = facet_vals[FACET_KEYS[0]], facet_vals[FACET_KEYS[1]]
    nr, nc = len(f1s), len(f0s)
    fig, axes = plt.subplots(nr, nc, figsize=(3.6 * nc, 3.2 * nr), squeeze=False)

    allvals = np.array([r[metric] for r in rows], float)
    vmax = np.nanmax(allvals)
    vmin = np.nanmin(allvals)
    if center_on_ref and not np.isnan(ref_val) and vmin < ref_val < vmax:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=ref_val, vmax=vmax)
    else:
        norm = None

    for ri, f1 in enumerate(f1s):          # MinPeakDistance → rows
        for ci, f0 in enumerate(f0s):      # SmoothSpan → cols
            ax = axes[ri, ci]
            A  = _grid_array(metric, f0, f1)
            im = ax.imshow(A, origin='lower', aspect='auto', cmap=cmap, norm=norm)
            ax.set_xticks(range(len(xs))); ax.set_xticklabels(xs, fontsize=6, rotation=90)
            ax.set_yticks(range(len(ys))); ax.set_yticklabels(ys, fontsize=7)
            ax.set_title(f'{FACET_KEYS[0]}={f0}, {FACET_KEYS[1]}={f1}', fontsize=8)
            if ri == nr - 1:
                ax.set_xlabel(AX_X, fontsize=7)
            if ci == 0:
                ax.set_ylabel(AX_Y, fontsize=7)
            # annotate cells close to reference
            if not np.isnan(ref_val):
                for yi in range(len(ys)):
                    for xi in range(len(xs)):
                        v = A[yi, xi]
                        if not np.isnan(v) and ref_val > 0 and v <= FLAG_RATIO * ref_val:
                            ax.text(xi, yi, '•', ha='center', va='center',
                                    fontsize=8, color='lime')

    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.ax.tick_params(labelsize=7)
    if not np.isnan(ref_val):
        cbar.ax.axhline(cbar.norm(ref_val), color='k', lw=1.5)
        cbar.set_label(f'{metric}  (black line = {REF_SES} ref = {ref_val:.4g})', fontsize=8)
    fig.suptitle(title + f'   (green • = within {FLAG_RATIO:g}× the {REF_SES} reference value)',
                 fontsize=12, fontweight='bold')
    f = os.path.join(SAVE_DIR, f'grid_{metric}_{RUN_STAMP}.png')
    fig.savefig(f, dpi=140, bbox_inches='tight')
    plt.show()
    return f


fpaths = []
fpaths.append(heatmap_figure('gap_frac', ref['gap_frac'],
              'Gap fraction across grid  (UNDER-detection / missed-sniff fingerprint)'))
fpaths.append(heatmap_figure('short_isi_frac', ref['short_isi_frac'],
              'Short-ISI fraction across grid  (OVER-detection fingerprint)'))
fpaths.append(heatmap_figure('mean_rate', ref['mean_rate'],
              'Mean inhalation rate (Hz) across grid'))

# ============================================================
# Two-sided operating point + summary
# ============================================================
# Self-contained cost: under-detection (gaps) + weighted over-detection (short ISIs).
# This needs NO reference — the best combo is the one that fills in missed sniffs
# without introducing physiologically-impossible ones.
for r in rows:
    g = r['gap_frac'] if not np.isnan(r['gap_frac']) else 1.0
    s = r['short_isi_frac'] if not np.isnan(r['short_isi_frac']) else 1.0
    r['cost'] = g + OD_WEIGHT * s
    # reference closeness (secondary, context only)
    rr = r['mean_rate'] / ref['mean_rate'] if ref['mean_rate'] else np.inf
    r['rate_ratio'] = rr

ranked = sorted(rows, key=lambda r: r['cost'])      # best operating point first
best   = ranked[0]

lines = []
lines.append(f"GRID two-sided detection sweep — {GRID_SES}")
lines.append(f"Window {ANALYSIS_WIN} s | max plausible {MAX_BR_HZ:g} Hz "
             f"(min ISI {MIN_ISI*1000:.0f} ms)")
lines.append(f"Gap metric: ISI > {GAP_MIN_RATIO:g}× local median (±{GAP_LOCAL_K} ISIs) "
             f"= missed sniff(s); cost = gap_frac + {OD_WEIGHT:g}×short_isi_frac")
lines.append(f"Reference {REF_SES} (context only): rate={ref['mean_rate']:.3f} Hz, "
             f"shortISI={ref['short_isi_frac']*100:.2f}%, gap={ref['gap_frac']*100:.2f}%")
lines.append('')
lines.append("=> BEST OPERATING POINT (lowest combined miss+false-positive cost):")
lines.append(f"   SmoothSpan={best['SmoothSpan']}, MinPeakDistance={best['MinPeakDistance']}, "
             f"MinPeakHeightFactor={best['MinPeakHeightFactor']}, "
             f"MinPeakPromFactor={best['MinPeakPromFactor']}")
lines.append(f"   rate={best['mean_rate']:.3f} Hz | gap(miss)={best['gap_frac']*100:.2f}% | "
             f"shortISI(false+)={best['short_isi_frac']*100:.2f}% | cost={best['cost']:.4f}")
lines.append('')

hdr = (f"{'Smooth':>6} {'MinDist':>7} {'HeightF':>7} {'PromF':>6}  "
       f"{'rate_Hz':>8} {'gap%':>7} {'shortISI%':>9} {'cost':>7} {'rate/ref':>8}")
lines.append('TOP 20 combinations by two-sided cost (lowest = best):')
lines.append(hdr)
lines.append('-' * len(hdr))
for r in ranked[:20]:
    lines.append(
        f"{r['SmoothSpan']:>6} {r['MinPeakDistance']:>7} {r['MinPeakHeightFactor']:>7} "
        f"{r['MinPeakPromFactor']:>6}  {r['mean_rate']:>8.3f} {r['gap_frac']*100:>6.2f}% "
        f"{r['short_isi_frac']*100:>8.2f}% {r['cost']:>7.4f} {r['rate_ratio']:>8.2f}")

summary = '\n'.join(lines)
print('\n' + summary)
sp = os.path.join(SAVE_DIR, f'grid_summary_{RUN_STAMP}.txt')
with open(sp, 'w', encoding='utf-8') as fout:
    fout.write(summary + '\n')

print('\nSaved:')
for f in fpaths:
    print('  ' + f)
print('  ' + sp)
print('Done.')
#%%
