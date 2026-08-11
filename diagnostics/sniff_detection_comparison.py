"""
Sniff detection comparison — individual trial visualisation
Compares inhalation-onset detection quality across configurable session groups,
trial by trial. Each group picks a detector (ml_inh_onsets or ft_inh_onsets), so
you can compare detectors against each other — including on the SAME session for a
like-for-like comparison (the fixed RNG seed makes equal sessions sample equal trials).
@author: Xander
"""
#%%
import sys
sys.path.append(r"C:\Users\xaand\Documents\PhD\Analysis\Informatics\Python\Electrophysiology\helpers")
import os
import glob
import pickle
from datetime import datetime
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ============================================================
# Paths
# ============================================================
DATA_PATH = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\data\raw_data\ofc_spks_newpreproc"
MAT_ROOT  = r"C:\Users\xaand\Documents\PhD\Experiments\postprocessed_files"
SAVE_DIR  = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\Analysis\sniff_detection_comparison"

# ============================================================
# CONFIGURATION — edit these blocks to change session groups
# ============================================================
# source='sniffs_pkl' → sessions already loaded from sniffs.pkl; provide session_indices (0-based)
# source='mat'        → loads {MAT_ROOT}/{id}/*_sniffdata_*.mat (any sniffdata_xxxx;
#                       prefers center when present); provide session_ids
#                       Works for ANY session in postprocessed_files: YYMMDD_KKxxx or YYMMDD_XTxxx
# source='grid'       → pulls one parameter combo from a grid-search session's ml_grid.
#                       provide grid_session=<folder> and params=dict(...) (see below).
#                       The gray cam_resp trace is the grid session's ml_cam_resp (same raw
#                       signal for every combo — only the inh-onset overlay changes).
#
# detector='ml'  → uses ml_inh_onsets (+ ml_cam_resp gray trace)   [default if omitted]
# detector='ft'  → uses ft_inh_onsets (+ ft_cam_resp gray trace)
#   To compare detectors on one session, make two groups with the same session and
#   different detectors — identical RNG_SEED + n_per → identical trials per group.
#   (detector is ignored for source='grid' — the grid stores only ml onsets.)
#
# GRID parameter values available in 250118_KK152_grid (390 combos):
#   SmoothSpan          : 1, 2, 3
#   MinPeakDistance     : 3, 5
#   MinPeakHeightFactor : 0.1, 0.15, 0.2, 0.25, 0.3
#   MinPeakPromFactor   : 0.4, 0.45, 0.5, ... 0.95, 1.0  (step 0.05)
#   Same RNG_SEED + n_per on two grid groups → identical trials → like-for-like overlay.
#
# Available sessions in postprocessed_files (copy IDs directly into session_ids):
#
#   KK064: 231020_KK064  231024_KK064
#   KK065: 231020_KK065  231024_KK065
#   KK066: 231020_KK066  231024_KK066
#   KK082: 231215_KK082  231219_KK082
#   KK084: 231215_KK084  231219_KK084
#   KK085: 231215_KK085  231219_KK085
#   KK087: 240201_KK087  240205_KK087
#   KK088: 240201_KK088  240205_KK088
#   KK089: 240201_KK089  240205_KK089
#   KK092: 240211_KK092  240215_KK092
#   KK093: 240211_KK093  240215_KK093
#   KK099: 240301_KK099  240305_KK099
#   KK100: 240301_KK100  240305_KK100
#
#   XT009: 260302_XT009  260303_XT009  260306_XT009  260307_XT009
#   XT010: 260302_XT010  260303_XT010  260306_XT010  260307_XT010
#   XT011: 260302_XT011  260303_XT011  260306_XT011  260307_XT011
#   XT012: 260302_XT012  260303_XT012  260306_XT012  260307_XT012
#   XT013: 260324_XT013  260325_XT013  260328_XT013  260329_XT013  260402_XT013
#   XT014: 260324_XT014  260325_XT014  260328_XT014  260329_XT014  260402_XT014
#   XT015: 260324_XT015  260325_XT015  260328_XT015  260329_XT015  260402_XT015
#   XT016: 260324_XT016  260325_XT016  260328_XT016  260329_XT016  260402_XT016
#   XT017: 260417_XT017  260418_XT017  260423_XT017  260424_XT017  260430_XT017
#   XT018: 260417_XT018  260418_XT018  260423_XT018  260424_XT018  260430_XT018
#   XT019: 260417_XT019  260418_XT019  260423_XT019  260424_XT019  260430_XT019
#

PRE_EVENT       = 4.0   # s before odor onset (assumed equal across all datasets)
N_PER_SES       = 20     # trials sampled per session for the figure
RNG_SEED        = 42
TRIALS_PER_PAGE = 40
XLIM_TRIAL      = (-PRE_EVENT, 12.0 - PRE_EVENT)   # full 12 s window (one panel per trial)

# detector → (inhalation-onset field, raw camera-response field)
DETECTOR_FIELDS = {
    'ml': ('ml_inh_onsets', 'ml_cam_resp'),
    'ft': ('ft_inh_onsets', 'ft_cam_resp'),
}

GRID_SESSION = '250118_KK152_grid'   # default grid folder for source='grid' groups

GROUPS = [
    # Three-way: STANDARD vs the two floor picks that tie at 8.75% miss / 0% impossible.
    # Eyeball whether the extra marks land on real cam_resp deflections, and whether
    # the most-sensitive pick (orange) adds anything over the conservative one (blue).
    dict(
        label='XT024',
        color='black',
        source='sess_ids',
        session_ids=['260619_XT024']
        #grid_session=GRID_SESSION,
        #params=dict(SmoothSpan=3, MinPeakDistance=5,
        #            MinPeakHeightFactor=0.25, MinPeakPromFactor=1.0),
    ),
    dict(
        label='XT025',
        color='black',
        source='sess_ids', #grid, sess_ids, or sniffs_pkl
        session_ids=['260619_XT025'],
        #grid_session=GRID_SESSION,
        #params=dict(SmoothSpan=3, MinPeakDistance=5,
        #            MinPeakHeightFactor=0.25, MinPeakPromFactor=1.0),
    ),
]

# ============================================================
# Load sniff data (sniffs.pkl needed only for source='sniffs_pkl' groups)
# ============================================================
with open(os.path.join(DATA_PATH, 'sniffs.pkl'), 'rb') as f:
    sniffs = pickle.load(f)

os.makedirs(SAVE_DIR, exist_ok=True)

SF    = int(sniffs[0]['samp_freq'])
T_VEC = np.arange(-PRE_EVENT, 12.0 - PRE_EVENT, 1 / SF)

# ============================================================
# Helpers
# ============================================================
def _inst_br(inh_frames, sf, pre):
    if inh_frames is None or len(inh_frames) < 2:
        return np.full(len(T_VEC), np.nan)
    inh_s = np.asarray(inh_frames, float) / sf - pre
    isis  = np.diff(inh_s)
    ok    = (isis > 0.05) & (isis < 2.0)
    if ok.sum() < 2:
        return np.full(len(T_VEC), np.nan)
    mid  = inh_s[:-1][ok] + isis[ok] / 2
    rate = 1.0 / isis[ok]
    return interp1d(mid, rate, kind='linear', bounds_error=False,
                    fill_value=np.nan)(T_VEC)


def _find_sniffdata_mat(d):
    """Return the path to a sniffdata mat in folder `d`, matching any
    *_sniffdata_*.mat (center, left, right, ...). The file is prefixed by the
    unique_id, which can differ from the folder name, so glob rather than assume.
    Prefer 'center' if present, else fall back to the first match."""
    hits = sorted(glob.glob(os.path.join(d, '*_sniffdata_*.mat')))
    if not hits:
        raise FileNotFoundError(f'No *_sniffdata_*.mat in {d}')
    for h in hits:
        if 'sniffdata_center' in os.path.basename(h):
            return h
    return hits[0]


def _load_mat(uid):
    p = _find_sniffdata_mat(os.path.join(MAT_ROOT, uid))
    return sio.loadmat(p, simplify_cells=True)


_GRID_CACHE = {}
def _load_grid_mat(grid_session):
    """Load (and cache) a grid session's sniffdata mat. The file is prefixed by the
    unique_id, which can differ from the folder name (e.g. '250118_KK152_grid' holds
    '250118_KK152_1_...'), so glob for it rather than assuming folder==prefix."""
    if grid_session not in _GRID_CACHE:
        p = _find_sniffdata_mat(os.path.join(MAT_ROOT, grid_session))
        _GRID_CACHE[grid_session] = sio.loadmat(p, simplify_cells=True)
    return _GRID_CACHE[grid_session]


def _grid_combo(grid_session, params):
    """Return (per-trial onsets, sf, pre, cam_resp list) for one grid parameter combo."""
    m    = _load_grid_mat(grid_session)
    for entry in m['ml_grid']:
        p = entry['params']
        if all(abs(float(p[k]) - float(v)) < 1e-9 for k, v in params.items()):
            sf   = int(m['samp_freq'])
            pre  = float(m.get('pre_event', PRE_EVENT))
            cam  = m['ml_cam_resp'] if 'ml_cam_resp' in m else None
            return entry['ml_inh_onsets'], sf, pre, cam
    raise ValueError(f"No grid combo matches {params} in {grid_session}. "
                     f"Example available combo: {m['ml_grid'][0]['params']}")


def _sample_group(group, n_per):
    """Sample n_per trials per session; returns list of (inh, raw, label, sf, pre)."""
    onset_fld, raw_fld = DETECTOR_FIELDS[group.get('detector', 'ml')]
    rng    = np.random.default_rng(RNG_SEED)
    sample = []
    if group['source'] == 'grid':
        onsets, sf, pre, cam = _grid_combo(group['grid_session'], group['params'])
        n_t    = len(onsets)
        chosen = rng.choice(n_t, size=min(n_per, n_t), replace=False)
        raws   = cam if cam is not None else [None] * n_t
        pstr   = ','.join(f'{k}={v}' for k, v in group['params'].items())
        for t in chosen:
            sample.append((onsets[t], raws[t], f'{group["grid_session"]} [{pstr}] t{t}',
                           sf, pre))
    elif group['source'] == 'sniffs_pkl':
        for sidx in group['session_indices']:
            ses    = sniffs[sidx]
            n_t    = len(ses[onset_fld])
            chosen = rng.choice(n_t, size=min(n_per, n_t), replace=False)
            uid    = ses.get('unique_id', ses.get('folder_identifier', f'ses{sidx}'))
            raws   = ses[raw_fld] if raw_fld in ses else [None] * n_t
            sf     = int(ses['samp_freq'])
            for t in chosen:
                sample.append((ses[onset_fld][t], raws[t], f'{uid} t{t}', sf, PRE_EVENT))
    else:
        for uid in group['session_ids']:
            m      = _load_mat(uid)
            n_t    = len(m[onset_fld])
            chosen = rng.choice(n_t, size=min(n_per, n_t), replace=False)
            sf     = int(m['samp_freq'])
            pre    = float(m.get('pre_event', PRE_EVENT))
            raws   = m[raw_fld] if raw_fld in m else [None] * n_t
            for t in chosen:
                sample.append((m[onset_fld][t], raws[t], f'{uid} t{t}', sf, pre))
    return sample


def _plot_panel(ax, inh, raw, sf, pre, xlim, color, row_lbl='', show_n_inh=True):
    if raw is not None and len(raw) > 0:
        raw_a = np.asarray(raw, float)
        t_raw = np.arange(len(raw_a)) / sf - pre
        bl_m  = (t_raw >= -3.5) & (t_raw < -0.5)
        mu = raw_a[bl_m].mean() if bl_m.any() else raw_a.mean()
        sd = raw_a[bl_m].std()  if bl_m.any() else raw_a.std()
        ax.plot(t_raw, (raw_a - mu) / (sd + 1e-9), color='gray', lw=0.7, alpha=0.8)

    inh_s = (np.asarray(inh, float) / sf - pre
             if (inh is not None and len(inh) > 0) else np.array([]))
    vis = inh_s[(inh_s >= xlim[0]) & (inh_s <= xlim[1])] if len(inh_s) else []
    for onset in vis:
        ax.axvline(onset, color=color, lw=0.8, alpha=0.5)

    ax2 = ax.twinx()
    br  = (_inst_br(inh, sf, pre)
           if (inh is not None and len(inh) > 0) else np.full(len(T_VEC), np.nan))
    ax2.plot(T_VEC, br, color='tomato', lw=1.5)
    ax2.set_ylim(0, 16)
    ax2.tick_params(axis='y', labelcolor='tomato', labelsize=5)
    ax2.set_ylabel('BR (Hz)', color='tomato', fontsize=5)

    ax.axvline(0, color='k', ls='--', lw=0.9)
    ax.set_xlim(xlim)
    ax.set_ylim(-8, 8)
    if row_lbl:
        ax.set_ylabel(row_lbl, fontsize=7)
    if show_n_inh:
        ax.set_title(f'{len(vis)} inh', fontsize=7, color='dimgray', pad=2)


# ============================================================
# Individual trial comparison — paginated
# ============================================================
#%%
n_g     = len(GROUPS)
n_cols  = n_g
samples = [_sample_group(g, N_PER_SES) for g in GROUPS]
n_rows  = min(len(s) for s in samples)

col_headers = [g['label'] for g in GROUPS]

RUN_STAMP = datetime.now().strftime('%Y%m%d_%H%M%S')  # unique per run → no overwrites

n_pages = max(1, int(np.ceil(n_rows / TRIALS_PER_PAGE)))

for page in range(n_pages):
    slc   = slice(page * TRIALS_PER_PAGE, (page + 1) * TRIALS_PER_PAGE)
    batch = [s[slc] for s in samples]
    n_r   = min(len(b) for b in batch)
    if n_r == 0:
        break

    fig, axes = plt.subplots(n_r, n_cols, figsize=(9.0 * n_cols, n_r * 2.2),
                             squeeze=False)

    for row in range(n_r):
        for gi, g in enumerate(GROUPS):
            if row >= len(batch[gi]):
                continue
            inh, raw, lbl, sf, pre = batch[gi][row]
            col = gi
            _plot_panel(axes[row, col], inh, raw, sf, pre, XLIM_TRIAL,
                        color=g['color'], row_lbl=lbl, show_n_inh=True)
            if row == 0:
                axes[row, col].set_title(col_headers[col], fontsize=8, fontweight='bold')
            if row == n_r - 1:
                axes[row, col].set_xlabel('Time from odor onset (s)', fontsize=8)

        for gi in range(n_g - 1):
            axes[row, gi].spines['right'].set_linewidth(2)
            axes[row, gi].spines['right'].set_color('dimgray')
            axes[row, gi + 1].spines['left'].set_linewidth(2)
            axes[row, gi + 1].spines['left'].set_color('dimgray')

    group_str = ' vs '.join(
        f"{g['label']} ({g['color']} {DETECTOR_FIELDS[g.get('detector', 'ml')][0]})"
        for g in GROUPS)
    fig.suptitle(
        f'Individual trial comparison — {group_str} — page {page+1}/{n_pages}\n'
        'Gray: cam_resp z-score  |  Red: BR  |  Colored marks: inh onsets (per group detector)'
        '  |  Panel title: # inh in window',
        fontsize=15, fontweight='bold', y=0.998, va='top')
    # Reserve top space for the 2-line suptitle so it never overlaps the first row.
    # Reserved fraction shrinks as the page gets taller (more rows).
    top_frac = 1.0 - min(0.08, 1.4 / (n_r * 2.2 + 1.4))
    fig.tight_layout(rect=[0, 0, 1, top_frac])
    fname = os.path.join(SAVE_DIR, f'cmp_trials_{RUN_STAMP}_p{page+1:02d}.png')
    fig.savefig(fname, dpi=120, bbox_inches='tight')
    plt.show()
    print(f"Saved: {fname}")

print("Done.")
#%%
