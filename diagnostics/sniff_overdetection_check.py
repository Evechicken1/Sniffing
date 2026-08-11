"""
Sniff over-detection check
==========================
Quantifies whether a new inhalation-detection setting OVER-detects sniffs
compared to a reference session, using physiological plausibility metrics
rather than trial-by-trial eyeballing.

Motivation
----------
A more sensitive peak detector (e.g. lowering the prominence to 0.70 sd) can
pick up spurious inhalation onsets. The tell-tale signatures of over-detection
are:
  1. Inter-sniff intervals (ISIs) that are physiologically impossible — mice
     cannot inhale faster than ~MAX_BR Hz, so a spike of very short ISIs is the
     fingerprint of double-counted / noise peaks.
  2. An inflated inhalation count per trial.
  3. An inflated mean / median breathing rate.

This script loads each session's `ml_inh_onsets` (always the ml_* fields — the
ft_* fields are never used for breathing analysis), restricts to an analysis
window, and compares the above metrics across configurable session groups.

Usage
-----
Edit the CONFIG block, then run. Outputs go to SAVE_DIR:
  overdetection_hist_<stamp>.png   — ISI + breathing-rate distributions
  overdetection_perTrial_<stamp>.png — per-trial count / rate distributions
  overdetection_summary_<stamp>.txt  — metric table
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

# ============================================================
# Paths
# ============================================================
MAT_ROOT = r"C:\Users\xaand\Documents\PhD\Experiments\postprocessed_files"
SAVE_DIR = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\Analysis\sniff_overdetection_check"

# ============================================================
# CONFIGURATION — edit these blocks
# ============================================================
# A group is either a whole session or one grid-search parameter combo:
#   source='mat'  (default) → loads {MAT_ROOT}/{session_id}/{session_id}_1_sniffdata_center.mat
#                             detector='ml' → ml_inh_onsets (default; always ml_* for breathing)
#                             detector='ft' → ft_inh_onsets (only for explicit detector compares)
#   source='grid'           → one combo from a grid session's ml_grid; provide
#                             grid_session=<folder> and params=dict(...)  (always ml onsets)
#
# GRID parameter values available in 250118_KK152_grid (390 combos):
#   SmoothSpan 1,2,3 | MinPeakDistance 3,5 | MinPeakHeightFactor 0.1,0.15,0.2,0.25,0.3 |
#   MinPeakPromFactor 0.4,0.45,...,0.95,1.0 (step 0.05)
#
# Put the CANDIDATE (new) setting first and the REFERENCE second — the summary
# flags the candidate as over-detecting if its short-ISI fraction or mean rate
# is materially above the reference.

GROUPS = [
    # CANDIDATE: highest detection rate that stays close to standard and clean.
    dict(label='h0.2, prom0.8  (highest rate)', color='darkorange',
         source='grid', grid_session='250118_KK152_grid',
         params=dict(SmoothSpan=3, MinPeakDistance=5,
                     MinPeakHeightFactor=0.2, MinPeakPromFactor=0.8)),
    # REFERENCE: the standard values (smooth=3, height=0.25, prom=1.0, MinDist=5).
    dict(label='STANDARD (h0.25, prom1.0)', color='black',
         source='grid', grid_session='250118_KK152_grid',
         params=dict(SmoothSpan=3, MinPeakDistance=5,
                     MinPeakHeightFactor=0.25, MinPeakPromFactor=1.0)),
    # Intermediate: only prominence lowered from standard.
    dict(label='h0.25, prom0.8', color='steelblue',
         source='grid', grid_session='250118_KK152_grid',
         params=dict(SmoothSpan=3, MinPeakDistance=5,
                     MinPeakHeightFactor=0.25, MinPeakPromFactor=0.8)),
]

PRE_EVENT       = 4.0          # s before odor onset (assumed if no pre_event in mat)
ANALYSIS_WIN    = (-4.0, 8.0)  # s rel. to odor onset over which sniffs are counted
MAX_BR_HZ       = 14.0         # max physiologically plausible inhalation rate (Hz)
#   → any ISI shorter than 1/MAX_BR_HZ is treated as "implausible" (over-detection)
MIN_BR_HZ       = 0.3          # below this an ISI is a long pause, excluded from rate stats
ISI_HIST_MAX    = 1.0          # s — x-axis cap for the ISI histogram
N_ISI_BINS      = 60
FLAG_RATIO      = 1.25         # candidate flagged if metric > FLAG_RATIO × reference

DETECTOR_FIELDS = {'ml': 'ml_inh_onsets', 'ft': 'ft_inh_onsets'}

# ============================================================
# Helpers
# ============================================================
def _load_mat(uid):
    p = os.path.join(MAT_ROOT, uid, f'{uid}_1_sniffdata_center.mat')
    return sio.loadmat(p, simplify_cells=True)


_GRID_CACHE = {}
def _load_grid_mat(grid_session):
    """Load (and cache) a grid session's sniffdata mat. The file is prefixed by the
    unique_id, which can differ from the folder name, so glob for it."""
    if grid_session not in _GRID_CACHE:
        d    = os.path.join(MAT_ROOT, grid_session)
        hits = glob.glob(os.path.join(d, '*_1_sniffdata_center.mat'))
        if not hits:
            raise FileNotFoundError(f'No *_1_sniffdata_center.mat in {d}')
        _GRID_CACHE[grid_session] = sio.loadmat(hits[0], simplify_cells=True)
    return _GRID_CACHE[grid_session]


def _grid_combo(grid_session, params):
    """Return (per-trial onsets, sf, pre) for one grid parameter combo."""
    m = _load_grid_mat(grid_session)
    for entry in m['ml_grid']:
        p = entry['params']
        if all(abs(float(p[k]) - float(v)) < 1e-9 for k, v in params.items()):
            return (entry['ml_inh_onsets'], int(m['samp_freq']),
                    float(m.get('pre_event', PRE_EVENT)))
    raise ValueError(f"No grid combo matches {params} in {grid_session}. "
                     f"Example available combo: {m['ml_grid'][0]['params']}")


def _group_onsets(group):
    """Resolve a group to (per-trial onsets, sf, pre, id_string), for either source."""
    if group.get('source') == 'grid':
        onsets, sf, pre = _grid_combo(group['grid_session'], group['params'])
        pstr = ','.join(f'{k}={v}' for k, v in group['params'].items())
        return onsets, sf, pre, f"{group['grid_session']} [{pstr}]"
    m   = _load_mat(group['session_id'])
    fld = DETECTOR_FIELDS[group.get('detector', 'ml')]
    return (m[fld], int(m['samp_freq']),
            float(m.get('pre_event', PRE_EVENT)), group['session_id'])


def _trial_onsets_s(inh_frames, sf, pre):
    """Frames → seconds relative to odor onset, restricted to ANALYSIS_WIN, sorted."""
    if inh_frames is None or len(inh_frames) == 0:
        return np.array([])
    s = np.sort(np.asarray(inh_frames, float) / sf - pre)
    return s[(s >= ANALYSIS_WIN[0]) & (s <= ANALYSIS_WIN[1])]


def collect(group):
    """Return per-trial metrics + pooled ISIs for one session/grid group."""
    onsets, sf, pre, id_str = _group_onsets(group)

    win_dur          = ANALYSIS_WIN[1] - ANALYSIS_WIN[0]
    min_isi          = 1.0 / MAX_BR_HZ
    max_isi          = 1.0 / MIN_BR_HZ

    per_trial_count  = []   # inhalations in window
    per_trial_rate   = []   # mean inhalation rate in window (count / win_dur)
    per_trial_shortf = []   # fraction of that trial's ISIs below min_isi
    all_isis         = []   # pooled ISIs (s)

    for t in range(len(onsets)):
        s = _trial_onsets_s(onsets[t], sf, pre)
        per_trial_count.append(len(s))
        per_trial_rate.append(len(s) / win_dur)
        if len(s) >= 2:
            isis = np.diff(s)
            all_isis.append(isis)
            per_trial_shortf.append(np.mean(isis < min_isi))
        else:
            per_trial_shortf.append(np.nan)

    all_isis = np.concatenate(all_isis) if all_isis else np.array([])
    valid    = all_isis[(all_isis >= min_isi) & (all_isis <= max_isi)]

    return dict(
        label   = group['label'],
        color   = group['color'],
        ses     = id_str,
        sf      = sf,
        n_trials= len(onsets),
        count   = np.array(per_trial_count, float),
        rate    = np.array(per_trial_rate, float),
        shortf  = np.array(per_trial_shortf, float),
        isis    = all_isis,
        # over-detection metrics
        short_isi_frac = float(np.mean(all_isis < min_isi)) if all_isis.size else np.nan,
        mean_count     = float(np.mean(per_trial_count)),
        mean_rate      = float(np.mean(per_trial_rate)),
        median_isi     = float(np.median(valid)) if valid.size else np.nan,
        inst_rate_med  = float(np.median(1.0 / valid)) if valid.size else np.nan,
    )


# ============================================================
# Run
# ============================================================
#%%
os.makedirs(SAVE_DIR, exist_ok=True)
RUN_STAMP = datetime.now().strftime('%Y%m%d_%H%M%S')

stats   = [collect(g) for g in GROUPS]
min_isi = 1.0 / MAX_BR_HZ

# ── Figure 1: distribution-level comparison ─────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

# (a) ISI histogram — over-detection shows as a spike left of the min-ISI line
for s in stats:
    if s['isis'].size:
        axes[0].hist(s['isis'], bins=np.linspace(0, ISI_HIST_MAX, N_ISI_BINS),
                     histtype='step', lw=1.8, color=s['color'], density=True,
                     label=f"{s['label']}  (n={s['isis'].size} ISIs)")
axes[0].axvline(min_isi, color='red', ls='--', lw=1.2,
                label=f'physiological min ISI = {min_isi*1000:.0f} ms ({MAX_BR_HZ:.0f} Hz)')
axes[0].set_xlabel('Inter-sniff interval (s)')
axes[0].set_ylabel('density')
axes[0].set_title('ISI distribution\n(mass left of red line = over-detection)', fontsize=10)
axes[0].legend(fontsize=7)
axes[0].tick_params(direction='in')

# (b) instantaneous breathing-rate distribution
for s in stats:
    if s['isis'].size:
        rate = 1.0 / s['isis'][s['isis'] > 0]
        axes[1].hist(rate, bins=np.linspace(0, 25, 50), histtype='step', lw=1.8,
                     color=s['color'], density=True, label=s['label'])
axes[1].axvline(MAX_BR_HZ, color='red', ls='--', lw=1.2,
                label=f'max plausible = {MAX_BR_HZ:.0f} Hz')
axes[1].set_xlabel('Instantaneous breathing rate (Hz)')
axes[1].set_ylabel('density')
axes[1].set_title('Instantaneous rate distribution\n(mass right of red line = over-detection)',
                  fontsize=10)
axes[1].legend(fontsize=7)
axes[1].tick_params(direction='in')

fig.suptitle('Sniff over-detection check — distributions', fontsize=13, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.95])
f1 = os.path.join(SAVE_DIR, f'overdetection_hist_{RUN_STAMP}.png')
fig.savefig(f1, dpi=150, bbox_inches='tight')
plt.show()

# ── Figure 2: per-trial spread ──────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
pos     = np.arange(len(stats))
labels  = [s['label'].split(':')[0] for s in stats]
colors  = [s['color'] for s in stats]

def _box(ax, data, title, ylabel):
    bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True,
                    showmeans=True, meanline=True)
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c); patch.set_alpha(0.4)
    for i, d in enumerate(data):
        jit = np.random.default_rng(0).normal(0, 0.04, len(d))
        ax.scatter(pos[i] + jit, d, s=8, color=colors[i], alpha=0.6, zorder=3)
    ax.set_xticks(pos); ax.set_xticklabels(labels, fontsize=8)
    ax.set_title(title, fontsize=10); ax.set_ylabel(ylabel)
    ax.tick_params(direction='in')

_box(axes[0], [s['count'] for s in stats],
     f'Inhalations per trial\n({ANALYSIS_WIN[0]:.0f}–{ANALYSIS_WIN[1]:.0f} s)', 'count')
_box(axes[1], [s['rate'] for s in stats], 'Mean inhalation rate per trial', 'rate (Hz)')
_box(axes[2], [s['shortf'][~np.isnan(s['shortf'])] for s in stats],
     f'Short-ISI fraction per trial\n(ISI < {min_isi*1000:.0f} ms)', 'fraction')

fig.suptitle('Sniff over-detection check — per-trial spread', fontsize=13, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.95])
f2 = os.path.join(SAVE_DIR, f'overdetection_perTrial_{RUN_STAMP}.png')
fig.savefig(f2, dpi=150, bbox_inches='tight')
plt.show()

# ── Summary table + verdict ─────────────────────────────────────────────────
lines = []
hdr = (f"{'session / setting':<28}  {'n_tr':>5}  {'inh/tr':>7}  {'rate_Hz':>8}  "
       f"{'med_ISI_ms':>10}  {'med_instHz':>10}  {'shortISI%':>9}")
lines.append(hdr)
lines.append('-' * len(hdr))
for s in stats:
    lines.append(
        f"{s['label']:<28}  {s['n_trials']:>5}  {s['mean_count']:>7.2f}  "
        f"{s['mean_rate']:>8.3f}  {s['median_isi']*1000:>10.1f}  "
        f"{s['inst_rate_med']:>10.2f}  {s['short_isi_frac']*100:>8.2f}%")

# verdict: candidate (first) vs reference (second)
if len(stats) >= 2:
    cand, ref = stats[0], stats[1]
    lines.append('')
    lines.append(f"Verdict (candidate '{cand['label']}' vs reference '{ref['label']}'):")
    checks = [
        ('mean inhalation rate', cand['mean_rate'], ref['mean_rate']),
        ('inhalations per trial', cand['mean_count'], ref['mean_count']),
        ('short-ISI fraction',
         cand['short_isi_frac'] if not np.isnan(cand['short_isi_frac']) else 0.0,
         ref['short_isi_frac']  if not np.isnan(ref['short_isi_frac'])  else 0.0),
    ]
    over = False
    for name, c, r in checks:
        # Elevation requires the candidate to actually EXCEED the reference; a
        # candidate at or below reference is never over-detecting (0-vs-0 is fine).
        if c <= r:
            ratio, flag = (c / r if r else 1.0), False
        elif r == 0:
            ratio, flag = float('inf'), True          # ref clean but candidate has spurious ISIs
        else:
            ratio = c / r
            flag  = ratio > FLAG_RATIO
        over |= flag
        lines.append(f"  {name:<24}: {c:.4g} vs {r:.4g}  "
                     f"(×{ratio:.2f}) {'<-- ELEVATED' if flag else 'ok'}")
    lines.append('')
    lines.append("  => OVER-DETECTION LIKELY" if over else
                 "  => no over-detection signal (candidate comparable to reference)")

summary = '\n'.join(lines)
print('\n' + summary)
sp = os.path.join(SAVE_DIR, f'overdetection_summary_{RUN_STAMP}.txt')
with open(sp, 'w', encoding='utf-8') as fout:
    fout.write(summary + '\n')

print(f"\nSaved:\n  {f1}\n  {f2}\n  {sp}")
print("Done.")
#%%
