#%% OFC-LDTg axonal stim — sniffing to novel/familiar odours (CONTROL animals)
#
# Control cohort for the OFC->LDTg opto experiment. Same novel/familiar odour
# paradigm + ChR2-style stimulation trials, but in control animals where the
# stimulation is not expected to drive a behavioural (sniffing) effect.
#
# Two session types are loaded together and every plot is labelled by type:
#   - Laser (10 mW)
#   - LED   (5 mW)
#
# Some sessions are missing on disk (e.g. the LED session for XT021); the loader
# below silently drops any session whose folder is missing or contains no .mat,
# and keeps the per-group animal lists aligned with what was actually imported.
#
# NOTE: all sniff/breathing quantities use the ml_* fields (ml_inh_onsets).

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\xaand\Documents\PhD\Analysis\Informatics\Python\Sniffing')
import sniff_tools as st
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import matplotlib.colors as mcolors

#%% Session IDs (control animals)

# Laser (10 mW) control sessions
laser_sess_ids = ['260528_XT020', '260528_XT021', '260528_XT022', '260528_XT023']
#laser_sess_ids = ['260528_XT020', '260528_XT023']

# LED (5 mW) control sessions  (XT021 LED session is missing on disk)
led_sess_ids = ['260602_XT020', '260602_XT022', '260602_XT023']
#led_sess_ids = ['260602_XT020', '260602_XT023']


DATA_ROOT = r"C:\Users\xaand\Documents\PhD\Experiments\postprocessed_files"
# DATA_ROOT = r"Y:\postprocessed_data"


def load_sessions(sess_ids, data_root=DATA_ROOT):
    """Import sniff .mat structs, dropping sessions whose folder is missing or empty.

    Returns (sniffs_list, kept_ids) with the two lists index-aligned, so a missing
    session never shifts the animal it is later paired with.
    """
    existing = [s for s in sess_ids if os.path.isdir(os.path.join(data_root, s))
                and glob.glob(os.path.join(data_root, s, '*.mat'))]
    missing = [s for s in sess_ids if s not in existing]
    if missing:
        print('  ! skipped (folder missing / no .mat):', missing)

    paths = [os.path.join(data_root, s) for s in existing]
    sniffs = st.import_sniff_mat_select(paths)

    # Defensive: drop anything that still imported empty.
    kept = [(s, sn) for s, sn in zip(existing, sniffs) if len(sn) > 0]
    kept_ids = [s for s, _ in kept]
    kept_sniffs = [sn for _, sn in kept]
    return kept_sniffs, kept_ids


print('Loading laser (10 mW) controls...')
sniffs_laser, ids_laser = load_sessions(laser_sess_ids)
print('Loading LED (5 mW) controls...')
sniffs_led, ids_led = load_sessions(led_sess_ids)

# Each group: (display label, sniffs list, session ids, shading colour)
groups = [
    ('Laser (10 mW)', sniffs_laser, ids_laser, 'red'),
    ('LED (5 mW)',    sniffs_led,   ids_led,   'red'),
]

for label, sniffs, ids, _ in groups:
    print(f'{label}: {len(sniffs)} sessions ->',
          [s.get('folder_identifier', i) for s, i in zip(sniffs, ids)])

#%% Parameters
savefig = 1
savefolder = r"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\Sniffing to nov fam\Controls"
os.makedirs(savefolder, exist_ok=True)

fps      = 713 / 12
t_range  = np.array([4.3, 7.3]) * fps     # odour-response window (frames)
bl_range = np.array([0.5, 4.0]) * fps     # pre-odour baseline window (frames)
t_time   = (t_range[1]  - t_range[0])  / fps
bl_time  = (bl_range[1] - bl_range[0]) / fps

# Condition styling (shared across plots)
COND_COLORS = {'nov': 'green', 'fam': 'purple', 'opto': 'violet'}
COND_LABELS = {'nov': 'Novel', 'fam': 'Familiar', 'opto': 'Familiar + stim'}


def slug(label):
    """Filename-safe group label, e.g. 'Laser (10 mW)' -> 'Laser_10mW'."""
    return (label.replace(' (', '_').replace(' ', '').replace('mW)', 'mW')
                 .replace(')', '').replace('(', ''))


#%% Per-mouse scalar response  [nov, fam, opto]  (baseline-subtracted inh/s)
def compute_rates(sniffs_list):
    nov_list, fam_list, opto_list = [], [], []
    for s in sniffs_list:
        tidx = s['trial_idx'] - 1
        nov_mask  = (s['trial_novelty']     == 1) & (s['trial_occur'] <= 3) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
        fam_mask  = (s['trial_familiarity'] == 1) & (s['trial_occur'] <= 3) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
        opto_mask = (s['trial_opto']        == 1) & (s['trial_occur'] <= 3) & (s['trial_blank'] == 0)

        def mean_delta(mask):
            idxs = tidx[mask]
            if not len(idxs):
                return np.nan
            onsets = s['ml_inh_onsets'][idxs]
            stim = np.array([((i > t_range[0])  & (i < t_range[1])).sum()  for i in onsets], float) / t_time
            bl   = np.array([((i > bl_range[0]) & (i < bl_range[1])).sum() for i in onsets], float) / bl_time
            return float((stim - bl).mean())

        nov_list.append(mean_delta(nov_mask))
        fam_list.append(mean_delta(fam_mask))
        opto_list.append(mean_delta(opto_mask))
    return np.array(nov_list), np.array(fam_list), np.array(opto_list)


#%% Mean sniffing change bar plot (Familiar / Fam+stim / Novel) — per group
def p_to_stars(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 5e-2 else 'ns'


def add_sig_bar(ax, x1, x2, y, h, stars):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c='k', clip_on=False, zorder=5)
    ax.text((x1 + x2) / 2, y + h * 1.05, stars, ha='center', va='bottom',
            clip_on=False, zorder=6)


for label, sniffs, ids, _ in groups:
    if not len(sniffs):
        continue
    nov, fam, opto = compute_rates(sniffs)
    col_data   = [fam, opto, nov]
    col_labels = ['Familiar', 'Familiar\n+ stim', 'Novel']
    bar_colors = ['purple', 'violet', 'green']
    n_mice = len(sniffs)

    fig, ax = plt.subplots(figsize=(4.5, 5), dpi=300)

    # per-mouse connecting traces + dots
    for m in range(n_mice):
        ax.plot([0, 1, 2], [fam[m], opto[m], nov[m]], color='gray', lw=1, alpha=0.75)
    for x, ys, c in [(0, fam, 'purple'), (1, opto, 'violet'), (2, nov, 'green')]:
        ax.scatter(np.full(n_mice, x), ys, color=c, s=28, zorder=3, alpha=0.85)

    # bars (mean ± SEM)
    for i, (vals, c) in enumerate(zip(col_data, bar_colors)):
        v = vals[np.isfinite(vals)]
        if not len(v):
            continue
        m_val = v.mean()
        sem   = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
        rgba  = mcolors.to_rgba(c)
        ax.bar(i, m_val, width=0.6, color=(*rgba[:3], 0.25),
               edgecolor=(0, 0, 0, 0.6), linewidth=1, zorder=0)
        ax.errorbar(i, m_val, yerr=sem, capsize=2, color='black', lw=1, zorder=2, alpha=0.7)

    # pairwise paired Wilcoxon, Bonferroni-corrected
    sig_pairs = [(fam, opto, 0, 1), (fam, nov, 0, 2), (opto, nov, 1, 2)]
    raw_ps = []
    for a, b, xi, xj in sig_pairs:
        mask = np.isfinite(a) & np.isfinite(b)
        raw_ps.append(stats.wilcoxon(a[mask], b[mask])[1] if mask.sum() >= 2 else np.nan)
    n_valid = sum(np.isfinite(p) for p in raw_ps) or 1

    finite = np.concatenate([v[np.isfinite(v)] for v in col_data if np.isfinite(v).any()])
    ymin, ymax = finite.min(), finite.max()
    yrng = max(ymax - ymin, 1.0)
    base_y, step_h, line_h = ymax + 0.03 * yrng, 0.08 * yrng, 0.015 * yrng

    sig_bars = []
    for (a, b, xi, xj), p in zip(sig_pairs, raw_ps):
        if np.isfinite(p):
            stars = p_to_stars(min(p * n_valid, 1.0))
            if stars != 'ns':
                sig_bars.append((xi, xj, stars))
    for k, (xi, xj, stars) in enumerate(sig_bars):
        add_sig_bar(ax, xi, xj, base_y + k * step_h, line_h, stars)

    ax.set_xticks(range(3))
    ax.set_xticklabels(col_labels, fontsize=14)
    ax.set_ylabel("Δ avg inhalations/sec", fontsize=14)
    ax.set_title(f"Mean sniffing change\nControls — {label}", pad=30, fontsize=14, weight='bold')
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_ylim(ymin - 0.05 * yrng, base_y + max(len(sig_bars), 1) * step_h + 0.1 * yrng)

    if savefig:
        plt.savefig(f"{savefolder}/Mean sniffing change_controls_{slug(label)}.png",
                    dpi=300, bbox_inches='tight')
    plt.show()


#%% Binned sniffing time trace (baseline-subtracted counts/bin) — across mice, per group
_BINS_FR  = np.linspace(0, 720, 13)        # 12 bins of ~60 frames
_BL_BINS  = list(range(0, 4))               # first 4 bins = pre-odour baseline
_BIN_CTRS = np.linspace(-3.5, 7.5, num=12)


def compute_time_traces_binned(sniffs_list):
    nov_tr, fam_tr, opto_tr = [], [], []
    for s in sniffs_list:
        tidx      = s['trial_idx'] - 1
        nov_mask  = (s['trial_novelty']     == 1) & (s['trial_occur'] <= 3) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
        fam_mask  = (s['trial_familiarity'] == 1) & (s['trial_occur'] <= 3) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
        opto_mask = (s['trial_opto']        == 1) & (s['trial_occur'] <= 3) & (s['trial_blank'] == 0)

        def get_trace(mask):
            idxs = tidx[mask]
            if not len(idxs):
                return np.full(12, np.nan)
            trials = np.array([np.histogram(s['ml_inh_onsets'][i], bins=_BINS_FR)[0].astype(float)
                               for i in idxs])
            bl = trials[:, _BL_BINS].mean(axis=1, keepdims=True)
            return (trials - bl).mean(axis=0)

        nov_tr.append(get_trace(nov_mask))
        fam_tr.append(get_trace(fam_mask))
        opto_tr.append(get_trace(opto_mask))
    return np.array(nov_tr), np.array(fam_tr), np.array(opto_tr)


for label, sniffs, ids, shade_c in groups:
    if not len(sniffs):
        continue
    nov_tr, fam_tr, opto_tr = compute_time_traces_binned(sniffs)
    specs = [
        (fam_tr,  'purple', 'Familiar',         'solid'),
        (opto_tr, 'violet', 'Familiar + stim',  'dashed'),
        (nov_tr,  'green',  'Novel',            'solid'),
    ]

    fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)
    for data, color, leg, ls in specs:
        valid = np.all(np.isfinite(data), axis=1)
        d = data[valid]
        if not len(d):
            continue
        mean = d.mean(axis=0)
        sem  = d.std(axis=0, ddof=1) / np.sqrt(len(d)) if len(d) > 1 else np.zeros_like(mean)
        ax.plot(_BIN_CTRS, mean, color=color, label=leg, ls=ls)
        ax.errorbar(_BIN_CTRS, mean, yerr=sem, fmt='o', color=color,
                    ecolor=color, elinewidth=1, capsize=3, ls=ls)

    ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
    ax.axvline(0, color='black', lw=1, alpha=0.7)
    ax.axvspan(0.3, 2.3, color=shade_c, alpha=0.3, linewidth=0)
    ax.set_xlabel("Time from odor onset (s)", fontsize=13)
    ax.set_ylabel("Δ avg inhalations (inh/bin)", fontsize=13)
    ax.set_title(f"Mean sniffing (first 3 presentations)\nControls — {label}",
                 pad=10, fontsize=13, weight='bold')
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()

    if savefig:
        plt.savefig(f"{savefolder}/Sniffing time trace binned_controls_{slug(label)}.png",
                    dpi=300, bbox_inches='tight')
    plt.show()


#%% Per-mouse binned time trace — per group
for label, sniffs, ids, shade_c in groups:
    if not len(sniffs):
        continue
    nov_tr, fam_tr, opto_tr = compute_time_traces_binned(sniffs)
    mouse_ids = [s.get('folder_identifier', i) for s, i in zip(sniffs, ids)]

    for m in range(len(sniffs)):
        fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)
        for data, color, leg, ls in [(fam_tr[m],  'purple', 'Familiar',        'solid'),
                                      (opto_tr[m], 'violet', 'Familiar + stim', 'dashed'),
                                      (nov_tr[m],  'green',  'Novel',           'solid')]:
            if np.isfinite(data).any():
                ax.plot(_BIN_CTRS, data, color=color, label=leg, ls=ls)
        ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
        ax.axvline(0, color='black', lw=1, alpha=0.7)
        ax.axvspan(0.3, 2.3, color=shade_c, alpha=0.3, linewidth=0)
        ax.set_xlabel("Time from odor onset (s)", fontsize=13)
        ax.set_ylabel("Δ avg inhalations (inh/bin)", fontsize=13)
        ax.set_title(f"Mean sniffing — {mouse_ids[m]} ({label})",
                     pad=10, fontsize=13, weight='bold')
        ax.spines[['right', 'top']].set_visible(False)
        ax.legend()
        if savefig:
            plt.savefig(f"{savefolder}/Sniffing time trace binned_{mouse_ids[m]}_{slug(label)}.png",
                        dpi=300, bbox_inches='tight')
        plt.show()


#%% Habituation: sniffing change per presentation number — per group
def compute_habituation(sniffs_list, n_occur=8):
    nov_hab, fam_hab, opto_hab = [], [], []
    for s in sniffs_list:
        tidx = s['trial_idx'] - 1
        nov_row, fam_row, opto_row = [], [], []
        for occ in range(1, n_occur + 1):
            nov_mask  = (s['trial_novelty']     == 1) & (s['trial_occur'] == occ) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
            fam_mask  = (s['trial_familiarity'] == 1) & (s['trial_occur'] == occ) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
            opto_mask = (s['trial_opto']        == 1) & (s['trial_occur'] == occ) & (s['trial_blank'] == 0)

            def mean_delta(mask):
                idxs = tidx[mask]
                if not len(idxs):
                    return np.nan
                onsets = s['ml_inh_onsets'][idxs]
                stim = np.array([((i > t_range[0]) & (i < t_range[1])).sum() for i in onsets], float) / t_time
                bl   = np.array([((i > bl_range[0]) & (i < bl_range[1])).sum() for i in onsets], float) / bl_time
                return float((stim - bl).mean())

            nov_row.append(mean_delta(nov_mask))
            fam_row.append(mean_delta(fam_mask))
            opto_row.append(mean_delta(opto_mask))
        nov_hab.append(nov_row)
        fam_hab.append(fam_row)
        opto_hab.append(opto_row)
    return np.array(nov_hab), np.array(fam_hab), np.array(opto_hab)


for label, sniffs, ids, _ in groups:
    if not len(sniffs):
        continue
    nov_hab, fam_hab, opto_hab = compute_habituation(sniffs)
    x_ticks = np.arange(1, 9)
    specs = [
        (fam_hab,  'purple', 'Familiar',        'solid'),
        (opto_hab, 'violet', 'Familiar + stim', 'dashed'),
        (nov_hab,  'green',  'Novel',           'solid'),
    ]

    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
    for data, color, leg, ls in specs:
        n_valid = np.sum(np.isfinite(data), axis=0).astype(float)
        n_valid[n_valid < 2] = np.nan
        mean = np.nanmean(data, axis=0)
        sem  = np.nanstd(data, axis=0, ddof=1) / np.sqrt(n_valid)
        ax.plot(x_ticks, mean, color=color, label=leg, ls=ls, lw=1.5, ms=5)
        ax.fill_between(x_ticks, mean - sem, mean + sem, color=color, alpha=0.2, linewidth=0)

    ax.axhline(0, color='black', lw=1, ls='--', alpha=0.6)
    ax.set_xlabel("Presentation #", fontsize=13)
    ax.set_ylabel("Δ sniffing rate (inh/s)", fontsize=13)
    ax.set_title(f"Habituation\nControls — {label}", pad=10, fontsize=13, weight='bold')
    ax.set_xticks(x_ticks)
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()

    if savefig:
        plt.savefig(f"{savefolder}/Habituation_controls_{slug(label)}.png",
                    dpi=300, bbox_inches='tight')
    plt.show()


#%% Sniffing to opto blanks vs non-opto blanks — per group
def compute_blank_rates(sniffs_list):
    """Per-mouse baseline-subtracted inh/s for blank+stim vs blank trials."""
    opto_bl_vals, nonopto_bl_vals = [], []
    for s in sniffs_list:
        tidx = s['trial_idx'] - 1
        opto_blanks     = (s['trial_blank'] == 1) & (s['trial_opto'] == 1)
        non_opto_blanks = (s['trial_blank'] == 1) & (s['trial_opto'] == 0)

        def mean_delta(mask):
            idxs = tidx[mask]
            if not len(idxs):
                return np.nan
            onsets = s['ml_inh_onsets'][idxs]
            stim = np.array([((i > t_range[0])  & (i < t_range[1])).sum()  for i in onsets], float) / t_time
            bl   = np.array([((i > bl_range[0]) & (i < bl_range[1])).sum() for i in onsets], float) / bl_time
            return float((stim - bl).mean())

        nonopto_bl_vals.append(mean_delta(non_opto_blanks))
        opto_bl_vals.append(mean_delta(opto_blanks))
    return np.array(nonopto_bl_vals), np.array(opto_bl_vals)


for label, sniffs, ids, _ in groups:
    if not len(sniffs):
        continue
    nonopto, opto = compute_blank_rates(sniffs)
    col_data   = [nonopto, opto]
    col_labels = ['Blank', 'Blank\n+ stim']
    bar_colors = ['gray', 'violet']
    n_mice = len(sniffs)

    fig, ax = plt.subplots(figsize=(4, 5), dpi=300)
    for m in range(n_mice):
        ax.plot([0, 1], [nonopto[m], opto[m]], color='gray', lw=1, alpha=0.75)
    for x, ys, c in [(0, nonopto, 'gray'), (1, opto, 'violet')]:
        ax.scatter(np.full(n_mice, x), ys, color=c, s=28, zorder=3, alpha=0.85)

    for i, (vals, c) in enumerate(zip(col_data, bar_colors)):
        v = vals[np.isfinite(vals)]
        if not len(v):
            continue
        m_val = v.mean()
        sem   = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
        rgba  = mcolors.to_rgba(c)
        ax.bar(i, m_val, width=0.6, color=(*rgba[:3], 0.25),
               edgecolor=(0, 0, 0, 0.6), linewidth=1, zorder=0)
        ax.errorbar(i, m_val, yerr=sem, capsize=2, color='black', lw=1, zorder=2, alpha=0.7)

    a, b = nonopto, opto
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() >= 2:
        stars = p_to_stars(stats.wilcoxon(a[mask], b[mask])[1])
        finite = np.concatenate([v[np.isfinite(v)] for v in col_data if np.isfinite(v).any()])
        ymin, ymax = finite.min(), finite.max()
        yrng = max(ymax - ymin, 1.0)
        add_sig_bar(ax, 0, 1, ymax + 0.03 * yrng, 0.015 * yrng, stars)
        ax.set_ylim(ymin - 0.05 * yrng, ymax + 0.18 * yrng)

    ax.set_xticks(range(2))
    ax.set_xticklabels(col_labels, fontsize=14)
    ax.set_ylabel("Δ avg inhalations/sec", fontsize=14)
    ax.set_title(f"Blank trials\nControls — {label}", pad=30, fontsize=14, weight='bold')
    ax.spines[['right', 'top']].set_visible(False)

    if savefig:
        plt.savefig(f"{savefolder}/Mean sniffing change_blank_controls_{slug(label)}.png",
                    dpi=300, bbox_inches='tight')
    plt.show()

# %%
