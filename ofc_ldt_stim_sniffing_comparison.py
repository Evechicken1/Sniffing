

#%%
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\xaand\Documents\PhD\Analysis\Informatics\Python\Sniffing')
import sniff_tools as st
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import matplotlib.colors as mcolors

#%% Session IDs

# LED (5 mW) sessions – on target
led_sess_ids = ['260302_XT009','260302_XT010','260302_XT011','260302_XT012',\
                '260324_XT013','260324_XT014','260324_XT015','260324_XT016',
                '260417_XT017','260417_XT018','260417_XT019',
                '260619_XT024','260619_XT025'   
                ]
# all LED: ['260302_XT009','260302_XT010','260302_XT011','260302_XT012',
#            '260324_XT013','260324_XT014','260324_XT015','260324_XT016']

# Laser (10 mW) sessions – on target
laser_sess_ids = ['260306_XT009','260306_XT010','260306_XT011','260306_XT012',\
                  '260328_XT013','260328_XT014','260328_XT015','260328_XT016',\
                  '260423_XT017','260423_XT018','260423_XT019',
                  '260611_XT024','260611_XT025']
# all laser: ['260306_XT009','260306_XT010','260306_XT011','260306_XT012',
#              '260328_XT013','260328_XT014','260328_XT015','260328_XT016']
# alternative 4-mouse set: ['260328_XT013','260328_XT014','260328_XT015','260328_XT016']

led_paths   = [rf"C:\Users\xaand\Documents\PhD\Experiments\postprocessed_files\\{s}" for s in led_sess_ids]
laser_paths = [fr"C:\Users\xaand\Documents\PhD\Experiments\postprocessed_files\\{s}" for s in laser_sess_ids]

sniffs_led   = st.import_sniff_mat_select(led_paths)
sniffs_laser = st.import_sniff_mat_select(laser_paths)

print('LED sessions:', len(sniffs_led), 'mice', [s['folder_identifier'] for s in sniffs_led])
print('Laser sessions:', len(sniffs_laser), 'mice', [s['folder_identifier'] for s in sniffs_laser])

#%% Parameters
savefig    = 1
savefolder = r"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\Sniffing to nov fam"

fps      = 713 / 12
t_range  = np.array([4.3, 7.3]) * fps
bl_range = np.array([0.5, 4.0]) * fps
t_time   = (t_range[1]  - t_range[0])  / fps
bl_time  = (bl_range[1] - bl_range[0]) / fps

#%% Compute per-mouse [nov, fam, opto] from a sniffs list
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

nov_led,   fam_led,   opto_led   = compute_rates(sniffs_led)
nov_laser, fam_laser, opto_laser = compute_rates(sniffs_laser)

#%% Bar plot of change in breathing rate between conditions

def p_to_stars(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 5e-2 else 'ns'

def add_sig_bar(ax, x1, x2, y, h, stars):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c='k', clip_on=False, zorder=5)
    ax.text((x1+x2)/2, y+h*1.05, stars, ha='center', va='bottom', clip_on=False, zorder=6)

graph      = ["Familiar", "Fam+5mW", "Fam+10mW", "Novel"]
bar_colors = ['purple', 'violet', "#e77295", 'green']

n_mice = len(nov_led)   # same animals in both sessions (paired)

# Average fam and nov across both days per mouse → one value per animal per column
fam_avg = (fam_led + fam_laser) / 2
nov_avg = (nov_led + nov_laser) / 2

# One row per mouse: [fam_avg, opto_led, opto_laser, nov_avg]
col_data = [fam_avg, opto_led, opto_laser, nov_avg]

fig, ax = plt.subplots(figsize=(6, 5), dpi=300)


# One connecting line per mouse through all 4 columns
for i in range(n_mice):
    ax.plot([0, 1, 2, 3],
            [fam_avg[i], opto_led[i], opto_laser[i], nov_avg[i]],
            color='gray', lw=1, alpha=0.75)

# One dot per mouse per column
for x, ys, c in [(0, fam_avg,   'purple'),
                 (1, opto_led,  'violet'),
                 (2, opto_laser,"#e77295"),
                 (3, nov_avg,   'green')]:
    ax.scatter(np.full(n_mice, x), ys, color=c, s=28, zorder=3, alpha=0.85)

# Bars (mean ± SEM)
for i, (vals, c) in enumerate(zip(col_data, bar_colors)):
    v = vals[np.isfinite(vals)]
    if not len(v):
        continue
    m_val = v.mean()
    sem   = v.std(ddof=1) / np.sqrt(len(v))
    rgba  = mcolors.to_rgba(c)
    ax.bar(i, m_val, width=0.6, color=(*rgba[:3], 0.25), zorder=0)
    ax.errorbar(i, m_val, yerr=sem, capsize=2, color='black', lw=1, zorder=2, alpha=0.7)

# Significance: all pairwise paired Wilcoxon (same n_mice animals across all columns)
sig_pairs = [
    (fam_avg,   opto_led,   0, 1),
    (fam_avg,   opto_laser, 0, 2),
    (fam_avg,   nov_avg,    0, 3),
    (opto_led,  opto_laser, 1, 2),
    (opto_led,  nov_avg,    1, 3),
    (opto_laser,nov_avg,    2, 3),
]

raw_ps = []
for a, b, xi, xj in sig_pairs:
    a, b = np.asarray(a), np.asarray(b)
    mask = np.isfinite(a) & np.isfinite(b)
    p = stats.wilcoxon(a[mask], b[mask])[1] if mask.sum() >= 2 else np.nan
    raw_ps.append(p)

# Bonferroni correction
n_valid = sum(np.isfinite(p) for p in raw_ps) or 1
ymin   = min(np.nanmin(v) for v in col_data if np.isfinite(v).any())
ymax   = max(np.nanmax(v) for v in col_data if np.isfinite(v).any())
yrng   = max(ymax - ymin, 1.0)
base_y = ymax  + 0.03 * yrng
step_h = 0.08  * yrng
line_h = 0.015 * yrng

sig_bars = []
for (a, b, xi, xj), p in zip(sig_pairs, raw_ps):
    p_adj = min(p * n_valid, 1.0) if np.isfinite(p) else np.nan
    stars = p_to_stars(p_adj) if np.isfinite(p_adj) else 'n/a'
    if stars not in ('ns', 'n/a'):
        sig_bars.append((xi, xj, stars))

for k, (xi, xj, stars) in enumerate(sig_bars):
    add_sig_bar(ax, xi, xj, base_y + k * step_h, line_h, stars)

ax.set_xticks(range(4))
ax.set_xticklabels(graph, fontsize=14)
ax.set_ylabel("Δ avg inhalations/sec", fontsize=14)
ax.set_title("Mean sniffing change", pad=40, fontsize=14, weight='bold')
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylim(ymin - 0.05*yrng, base_y + max(len(sig_bars), 1)*step_h + 0.1*yrng)

if savefig:
    plt.savefig(f"{savefolder}/Mean sniffing change_LED_vs_laser.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Mean sniffing change_LED_vs_laser.svg", bbox_inches='tight')
plt.show()


#%% Time trace: instantaneous sniffing rate via Gaussian kernel smoothing

# Time axis: bins are 60 frames wide, so frame/60 - 4.0 maps frame 240 → t=0 (odor onset)
_N_FR   = 720
_T_ARR  = np.arange(_N_FR) / 60.0 - 4.0   # seconds relative to odor onset
_BL_MSK = _T_ARR < 0                        # pre-odor baseline frames

def compute_time_traces(sniffs_list, sigma_s=0.15):
    """
    Convolve each trial's spike train with a Gaussian (sigma=sigma_s seconds),
    scale to Hz, baseline-subtract (pre-odor mean), then average across trials
    and return one (n_mice, _N_FR) array per condition.
    """
    sigma_fr = sigma_s * fps
    nov_tr, fam_tr, opto_tr = [], [], []
    for s in sniffs_list:
        tidx      = s['trial_idx'] - 1
        nov_mask  = (s['trial_novelty']     == 1) & (s['trial_occur'] <= 3) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
        fam_mask  = (s['trial_familiarity'] == 1) & (s['trial_occur'] <= 3) & (s['trial_opto'] == 0) & (s['trial_blank'] == 0)
        opto_mask = (s['trial_opto']        == 1) & (s['trial_occur'] <= 3) & (s['trial_blank'] == 0)

        def get_trace(mask):
            idxs = tidx[mask]
            if not len(idxs):
                return np.full(_N_FR, np.nan)
            trial_rates = []
            for i in idxs:
                spike = np.zeros(_N_FR)
                ons   = np.asarray(s['ml_inh_onsets'][i]).astype(int)
                ons   = ons[(ons >= 0) & (ons < _N_FR)]
                spike[ons] = 1.0
                trial_rates.append(gaussian_filter1d(spike, sigma=sigma_fr) * fps)
            trial_rates = np.array(trial_rates)           # (n_trials, _N_FR)
            bl = trial_rates[:, _BL_MSK].mean(axis=1, keepdims=True)
            return (trial_rates - bl).mean(axis=0)

        nov_tr.append(get_trace(nov_mask))
        fam_tr.append(get_trace(fam_mask))
        opto_tr.append(get_trace(opto_mask))
    return np.array(nov_tr), np.array(fam_tr), np.array(opto_tr)

nov_tr_led,   fam_tr_led,   opto_tr_led   = compute_time_traces(sniffs_led)
nov_tr_laser, fam_tr_laser, opto_tr_laser = compute_time_traces(sniffs_laser)

# Average familiar and novel across both days per mouse (same animals)
fam_tr_avg = (fam_tr_led  + fam_tr_laser) / 2   # (n_mice, _N_FR)
nov_tr_avg = (nov_tr_led  + nov_tr_laser) / 2   # (n_mice, _N_FR)

trace_specs = [
    (fam_tr_avg,    'purple',     'Familiar',         'solid'),
    (opto_tr_led,   'violet',     'Fam+5mW',    'dashed'),
    (opto_tr_laser, "#e77295",   'Fam+10mW', 'dashed'),
    (nov_tr_avg,    'green',      'Novel',            'solid'),
]

fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)

for data, color, label, ls in trace_specs:
    valid = np.all(np.isfinite(data), axis=1)
    d = data[valid]
    if not len(d):
        continue
    mean = d.mean(axis=0)
    sem  = d.std(axis=0, ddof=1) / np.sqrt(len(d))
    ax.plot(_T_ARR, mean, color=color, label=label, ls=ls, lw=1.5)
    ax.fill_between(_T_ARR, mean - sem, mean + sem, color=color, alpha=0.2)

ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
ax.axvline(0, color='black', lw=1, alpha=0.7)
ax.axvspan(0.3, 2.3, 0, 0.1, color='red', alpha=0.3, linewidth=0)

ax.set_xlabel("Time from odor onset (s)", fontsize=13)
ax.set_ylabel("Δ sniffing rate (Hz)", fontsize=13)
ax.set_title("Mean sniffing (first 3 presentations)", pad=10, fontsize=13, weight='bold')
ax.spines[['right', 'top']].set_visible(False)
ax.legend()

if savefig:
    plt.savefig(f"{savefolder}/Sniffing time trace_LED_vs_laser.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Sniffing time trace_LED_vs_laser.svg", bbox_inches='tight')
plt.show()

#%% Time trace: binned sniffing (histogram, baseline-subtracted counts per bin)

_BINS_FR  = np.linspace(0, 720, 13)   # 12 bins of ~60 frames each
_BL_BINS  = list(range(0, 4))          # first 4 bins = pre-odor baseline
_BIN_CTRS = np.linspace(-3.5, 7.5, num=12)

def compute_time_traces_binned(sniffs_list):
    """Per-mouse (12,) baseline-subtracted bin-count arrays per condition."""
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

nov_btr_led,   fam_btr_led,   opto_btr_led   = compute_time_traces_binned(sniffs_led)
nov_btr_laser, fam_btr_laser, opto_btr_laser = compute_time_traces_binned(sniffs_laser)

fam_btr_avg = (fam_btr_led + fam_btr_laser) / 2
nov_btr_avg = (nov_btr_led + nov_btr_laser) / 2

btrace_specs = [
    (fam_btr_avg,    'purple',   'Familiar',         'solid'),
    (opto_btr_led,   'violet',   'Fam+5mW',    'solid'),
    (opto_btr_laser, '#e77295', 'Fam+10mW', 'dashed'),
    (nov_btr_avg,    'green',    'Novel',            'solid'),
]

fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)

for data, color, label, ls in btrace_specs:
    valid = np.all(np.isfinite(data), axis=1)
    d = data[valid]
    if not len(d):
        continue
    mean = d.mean(axis=0)
    sem  = d.std(axis=0, ddof=1) / np.sqrt(len(d))
    ax.plot(_BIN_CTRS, mean, color=color, label=label, ls=ls)
    ax.errorbar(_BIN_CTRS, mean, yerr=sem, fmt='o', color=color,
                ecolor=color, elinewidth=1, capsize=3, ls=ls)

ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
ax.axvline(0, color='black', lw=1, alpha=0.7)
ax.axvspan(0.3, 2.3, color='red', alpha=0.3, linewidth=0)

ax.set_xlabel("Time from odor onset (s)", fontsize=13)
ax.set_ylabel("Δ avg inhalations (inh/bin)", fontsize=13)
ax.set_title("Mean sniffing (first 3 presentations)", pad=10, fontsize=13, weight='bold')
ax.spines[['right', 'top']].set_visible(False)
ax.legend()

if savefig:
    plt.savefig(f"{savefolder}/Sniffing time trace binned_LED_vs_laser.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Sniffing time trace binned_LED_vs_laser.svg", bbox_inches='tight')
plt.show()

#%% Time trace: binned sniffing per mouse individually

mouse_ids = [s['folder_identifier'] for s in sniffs_led]

for m in range(fam_btr_avg.shape[0]):
    fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)

    pm_specs = [
        (fam_btr_avg[m],    'purple',   'Familiar',         'solid'),
        (opto_btr_led[m],   'violet',   'Fam+5mW',    'solid'),
        (opto_btr_laser[m], '#e77295', 'Fam+10mW', 'dashed'),
        (nov_btr_avg[m],    'green',    'Novel',            'solid'),
    ]

    for data, color, label, ls in pm_specs:
        if np.isfinite(data).any():
            ax.plot(_BIN_CTRS, data, color=color, label=label, ls=ls)

    ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
    ax.axvline(0, color='black', lw=1, alpha=0.7)
    ax.axvspan(0.3, 2.3, color='red', alpha=0.3, linewidth=0)

    ax.set_xlabel("Time from odor onset (s)", fontsize=13)
    ax.set_ylabel("Δ avg inhalations (inh/bin)", fontsize=13)
    mid = mouse_ids[m] if m < len(mouse_ids) else f"mouse {m}"
    ax.set_title(f"Mean sniffing – {mid}", pad=10, fontsize=13, weight='bold')
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()

    if savefig:
        plt.savefig(f"{savefolder}/Sniffing time trace binned_{mid}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{savefolder}/Sniffing time trace binned_{mid}.svg", bbox_inches='tight')
    plt.show()

#%% Habituation: sniffing change per presentation number

def compute_habituation(sniffs_list, n_occur=8):
    """Per-mouse (n_occur,) arrays of baseline-subtracted sniffing rate per occurrence."""
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

nov_hab_led,   fam_hab_led,   opto_hab_led   = compute_habituation(sniffs_led)
nov_hab_laser, fam_hab_laser, opto_hab_laser = compute_habituation(sniffs_laser)

fam_hab_avg = (fam_hab_led + fam_hab_laser) / 2   # (n_mice, 8)
nov_hab_avg = (nov_hab_led + nov_hab_laser) / 2   # (n_mice, 8) — NaN beyond occurrence 3

x_ticks = np.arange(1, 9)

hab_specs = [
    (fam_hab_avg,    'purple',   'Familiar',         'solid'),
    (opto_hab_led,   'violet',   'Fam+5mW',    'dashed'),
    (opto_hab_laser, '#e77295', 'Fam+10mW', 'dashed'),
    (nov_hab_avg,    'green',    'Novel',            'solid'),
]

fig, ax = plt.subplots(figsize=(6, 4), dpi=300)

for data, color, label, ls in hab_specs:
    n_valid = np.sum(np.isfinite(data), axis=0).astype(float)
    n_valid[n_valid < 2] = np.nan
    mean = np.nanmean(data, axis=0)
    sem  = np.nanstd(data, axis=0, ddof=1) / np.sqrt(n_valid)
    ax.plot(x_ticks, mean, color=color, label=label, ls=ls, lw=1.5, ms=5)
    ax.errorbar(x_ticks, mean, yerr=sem, fmt='none', ecolor=color,
                elinewidth=1, capsize=3)

ax.axhline(0, color='black', lw=1, ls='--', alpha=0.6)
ax.set_xlabel("Presentation #", fontsize=13)
ax.set_ylabel("Δ sniffing rate (inh/s)", fontsize=13)
ax.set_title("Habituation", pad=10, fontsize=13, weight='bold')
ax.set_xticks(x_ticks)
ax.spines[['right', 'top']].set_visible(False)
ax.legend()

if savefig:
    plt.savefig(f"{savefolder}/Habituation_LED_vs_laser.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Habituation_LED_vs_laser.svg", bbox_inches='tight')
plt.show()


# =====================================================================
# ============================  CONTROLS  =============================
# =====================================================================
# Everything below ADDS the control cohort. It reuses the compute_* helpers
# defined above (they are group-agnostic) and only creates NEW figures — none
# of the experimental plots above are modified.
#
# Control LED has fewer sessions than control laser (the XT021 LED session is
# missing), so the two are NOT index-aligned. All control averaging across
# laser+LED days is therefore matched by ANIMAL NAME, and missing sessions are
# silently dropped at load time.

#%% Load control sessions (laser + LED), dropping missing/empty folders
import os, glob

DATA_ROOT = r"C:\Users\xaand\Documents\PhD\Experiments\postprocessed_files"

# Control cohort session IDs
laser_ctrl_ids_req = ['260528_XT020', '260528_XT021', '260528_XT022', '260528_XT023']
#laser_ctrl_ids_req = ['260528_XT020', '260528_XT023']

led_ctrl_ids_req   = ['260602_XT020', '260602_XT022', '260602_XT023']  # XT021 LED missing
#led_ctrl_ids_req   = ['260602_XT020', '260602_XT023']  # XT021,22 LED missing


def load_ctrl_sessions(sess_ids, data_root=DATA_ROOT):
    """Import sniff structs, dropping sessions whose folder is missing/empty.
    Returns (sniffs_list, kept_ids) kept index-aligned."""
    existing = [s for s in sess_ids if os.path.isdir(os.path.join(data_root, s))
                and glob.glob(os.path.join(data_root, s, '*.mat'))]
    missing = [s for s in sess_ids if s not in existing]
    if missing:
        print('  ! controls skipped (folder missing / no .mat):', missing)
    sniffs = st.import_sniff_mat_select([os.path.join(data_root, s) for s in existing])
    kept = [(s, sn) for s, sn in zip(existing, sniffs) if len(sn) > 0]
    return [sn for _, sn in kept], [s for s, _ in kept]


sniffs_laser_ctrl, laser_ctrl_ids = load_ctrl_sessions(laser_ctrl_ids_req)
sniffs_led_ctrl,   led_ctrl_ids   = load_ctrl_sessions(led_ctrl_ids_req)

print('Control LED sessions:  ', len(sniffs_led_ctrl),  [s['folder_identifier'] for s in sniffs_led_ctrl])
print('Control laser sessions:', len(sniffs_laser_ctrl), [s['folder_identifier'] for s in sniffs_laser_ctrl])


def _animal(sess_id):
    """'260528_XT020' -> 'XT020'."""
    return sess_id.split('_')[-1]


def align_by_animal(vals, sess_ids, master):
    """Reindex per-session values onto a master animal list (NaN where absent).
    Works for scalar-per-session (n,) and vector-per-session (n, k) arrays."""
    vals = np.asarray(vals, float)
    out = np.full((len(master),) + vals.shape[1:], np.nan)
    idx = {_animal(s): i for i, s in enumerate(sess_ids)}
    for j, a in enumerate(master):
        if a in idx:
            out[j] = vals[idx[a]]
    return out


# Master animal list = union of laser + LED control animals
ctrl_animals = sorted(set(_animal(s) for s in laser_ctrl_ids) |
                      set(_animal(s) for s in led_ctrl_ids))
print('Control animals:', ctrl_animals)

#%% Control: per-mouse scalar rates, traces and habituation (reuse helpers)
nov_led_c,   fam_led_c,   opto_led_c   = compute_rates(sniffs_led_ctrl)
nov_laser_c, fam_laser_c, opto_laser_c = compute_rates(sniffs_laser_ctrl)

# Align onto the master animal list, then average fam/nov across available days
def _avg_days(led_vals, led_ids, laser_vals, laser_ids):
    a = align_by_animal(led_vals,   led_ids,   ctrl_animals)
    b = align_by_animal(laser_vals, laser_ids, ctrl_animals)
    with np.errstate(invalid='ignore'):
        return np.nanmean(np.stack([a, b]), axis=0)

fam_avg_c  = _avg_days(fam_led_c, led_ctrl_ids, fam_laser_c, laser_ctrl_ids)
nov_avg_c  = _avg_days(nov_led_c, led_ctrl_ids, nov_laser_c, laser_ctrl_ids)
opto_led_c_al   = align_by_animal(opto_led_c,   led_ctrl_ids,   ctrl_animals)
opto_laser_c_al = align_by_animal(opto_laser_c, laser_ctrl_ids, ctrl_animals)

#%% Control mirror — Mean sniffing change bar (Familiar / Fam+5mW / Fam+10mW / Novel)
graph_c      = ["Familiar", "Fam+5mW", "Fam+10mW", "Novel"]
bar_colors_c = ['purple', 'violet', "#e77295", 'green']
col_data_c   = [fam_avg_c, opto_led_c_al, opto_laser_c_al, nov_avg_c]
n_ctrl       = len(ctrl_animals)

fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

# one connecting line per animal (NaN -> gap, so missing days are skipped)
for i in range(n_ctrl):
    ax.plot([0, 1, 2, 3],
            [fam_avg_c[i], opto_led_c_al[i], opto_laser_c_al[i], nov_avg_c[i]],
            color='gray', lw=1, alpha=0.75)
for x, ys, c in [(0, fam_avg_c, 'purple'), (1, opto_led_c_al, 'violet'),
                 (2, opto_laser_c_al, "#e77295"), (3, nov_avg_c, 'green')]:
    ax.scatter(np.full(n_ctrl, x), ys, color=c, s=28, zorder=3, alpha=0.85)

for i, (vals, c) in enumerate(zip(col_data_c, bar_colors_c)):
    v = vals[np.isfinite(vals)]
    if not len(v):
        continue
    m_val = v.mean()
    sem   = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    rgba  = mcolors.to_rgba(c)
    ax.bar(i, m_val, width=0.6, color=(*rgba[:3], 0.25), zorder=0)
    ax.errorbar(i, m_val, yerr=sem, capsize=2, color='black', lw=1, zorder=2, alpha=0.7)

# pairwise paired Wilcoxon (animals present in both columns), Bonferroni-corrected
sig_pairs_c = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
raw_ps_c = []
for xi, xj in sig_pairs_c:
    a, b = col_data_c[xi], col_data_c[xj]
    mask = np.isfinite(a) & np.isfinite(b)
    raw_ps_c.append(stats.wilcoxon(a[mask], b[mask])[1] if mask.sum() >= 2 else np.nan)
n_valid_c = sum(np.isfinite(p) for p in raw_ps_c) or 1

finite_c = np.concatenate([v[np.isfinite(v)] for v in col_data_c if np.isfinite(v).any()])
ymin_c, ymax_c = finite_c.min(), finite_c.max()
yrng_c = max(ymax_c - ymin_c, 1.0)
base_y_c, step_h_c, line_h_c = ymax_c + 0.03 * yrng_c, 0.08 * yrng_c, 0.015 * yrng_c

sig_bars_c = []
for (xi, xj), p in zip(sig_pairs_c, raw_ps_c):
    if np.isfinite(p):
        stars = p_to_stars(min(p * n_valid_c, 1.0))
        if stars != 'ns':
            sig_bars_c.append((xi, xj, stars))
for k, (xi, xj, stars) in enumerate(sig_bars_c):
    add_sig_bar(ax, xi, xj, base_y_c + k * step_h_c, line_h_c, stars)

ax.set_xticks(range(4))
ax.set_xticklabels(graph_c, fontsize=14)
ax.set_ylabel("Δ avg inhalations/sec", fontsize=14)
ax.set_title("Mean sniffing change — CONTROLS", pad=40, fontsize=14, weight='bold')
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylim(ymin_c - 0.05 * yrng_c, base_y_c + max(len(sig_bars_c), 1) * step_h_c + 0.1 * yrng_c)

if savefig:
    plt.savefig(f"{savefolder}/Mean sniffing change_CONTROLS.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Mean sniffing change_CONTROLS.svg", bbox_inches='tight')
plt.show()

#%% Control mirror — Gaussian-smoothed sniffing time trace
nov_tr_led_c,   fam_tr_led_c,   opto_tr_led_c   = compute_time_traces(sniffs_led_ctrl)
nov_tr_laser_c, fam_tr_laser_c, opto_tr_laser_c = compute_time_traces(sniffs_laser_ctrl)

fam_tr_avg_c  = _avg_days(fam_tr_led_c, led_ctrl_ids, fam_tr_laser_c, laser_ctrl_ids)
nov_tr_avg_c  = _avg_days(nov_tr_led_c, led_ctrl_ids, nov_tr_laser_c, laser_ctrl_ids)
opto_tr_led_c_al   = align_by_animal(opto_tr_led_c,   led_ctrl_ids,   ctrl_animals)
opto_tr_laser_c_al = align_by_animal(opto_tr_laser_c, laser_ctrl_ids, ctrl_animals)

trace_specs_c = [
    (fam_tr_avg_c,        'purple',  'Familiar',          'solid'),
    (opto_tr_led_c_al,    'violet',  'Fam+5mW',     'dashed'),
    (opto_tr_laser_c_al,  "#e77295", 'Fam+10mW',  'dashed'),
    (nov_tr_avg_c,        'green',   'Novel',             'solid'),
]

fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)
for data, color, label, ls in trace_specs_c:
    valid = np.all(np.isfinite(data), axis=1)
    d = data[valid]
    if not len(d):
        continue
    mean = d.mean(axis=0)
    sem  = d.std(axis=0, ddof=1) / np.sqrt(len(d)) if len(d) > 1 else np.zeros_like(mean)
    ax.plot(_T_ARR, mean, color=color, label=label, ls=ls, lw=1.5)
    ax.fill_between(_T_ARR, mean - sem, mean + sem, color=color, alpha=0.2)

ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
ax.axvline(0, color='black', lw=1, alpha=0.7)
ax.axvspan(0.3, 2.3, 0, 0.1, color='red', alpha=0.3, linewidth=0)
ax.set_xlabel("Time from odor onset (s)", fontsize=13)
ax.set_ylabel("Δ sniffing rate (Hz)", fontsize=13)
ax.set_title("Mean sniffing (first 3 presentations) — CONTROLS", pad=10, fontsize=13, weight='bold')
ax.spines[['right', 'top']].set_visible(False)
ax.legend()

if savefig:
    plt.savefig(f"{savefolder}/Sniffing time trace_CONTROLS.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Sniffing time trace_CONTROLS.svg", bbox_inches='tight')
plt.show()

#%% Control mirror — Binned sniffing time trace
nov_btr_led_c,   fam_btr_led_c,   opto_btr_led_c   = compute_time_traces_binned(sniffs_led_ctrl)
nov_btr_laser_c, fam_btr_laser_c, opto_btr_laser_c = compute_time_traces_binned(sniffs_laser_ctrl)

fam_btr_avg_c = _avg_days(fam_btr_led_c, led_ctrl_ids, fam_btr_laser_c, laser_ctrl_ids)
nov_btr_avg_c = _avg_days(nov_btr_led_c, led_ctrl_ids, nov_btr_laser_c, laser_ctrl_ids)
opto_btr_led_c_al   = align_by_animal(opto_btr_led_c,   led_ctrl_ids,   ctrl_animals)
opto_btr_laser_c_al = align_by_animal(opto_btr_laser_c, laser_ctrl_ids, ctrl_animals)

btrace_specs_c = [
    (fam_btr_avg_c,       'purple',  'Familiar',          'solid'),
    (opto_btr_led_c_al,   'violet',  'Fam+5mW',     'solid'),
    (opto_btr_laser_c_al, '#e77295', 'Fam+10mW',  'dashed'),
    (nov_btr_avg_c,       'green',   'Novel',             'solid'),
]

fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)
for data, color, label, ls in btrace_specs_c:
    valid = np.all(np.isfinite(data), axis=1)
    d = data[valid]
    if not len(d):
        continue
    mean = d.mean(axis=0)
    sem  = d.std(axis=0, ddof=1) / np.sqrt(len(d)) if len(d) > 1 else np.zeros_like(mean)
    ax.plot(_BIN_CTRS, mean, color=color, label=label, ls=ls)
    ax.errorbar(_BIN_CTRS, mean, yerr=sem, fmt='o', color=color,
                ecolor=color, elinewidth=1, capsize=3, ls=ls)

ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
ax.axvline(0, color='black', lw=1, alpha=0.7)
ax.axvspan(0.3, 2.3, color='red', alpha=0.3, linewidth=0)
ax.set_xlabel("Time from odor onset (s)", fontsize=13)
ax.set_ylabel("Δ avg inhalations (inh/bin)", fontsize=13)
ax.set_title("Mean sniffing (first 3 presentations) — CONTROLS", pad=10, fontsize=13, weight='bold')
ax.spines[['right', 'top']].set_visible(False)
ax.legend()

if savefig:
    plt.savefig(f"{savefolder}/Sniffing time trace binned_CONTROLS.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Sniffing time trace binned_CONTROLS.svg", bbox_inches='tight')
plt.show()

#%% Control mirror — Habituation
nov_hab_led_c,   fam_hab_led_c,   opto_hab_led_c   = compute_habituation(sniffs_led_ctrl)
nov_hab_laser_c, fam_hab_laser_c, opto_hab_laser_c = compute_habituation(sniffs_laser_ctrl)

fam_hab_avg_c = _avg_days(fam_hab_led_c, led_ctrl_ids, fam_hab_laser_c, laser_ctrl_ids)
nov_hab_avg_c = _avg_days(nov_hab_led_c, led_ctrl_ids, nov_hab_laser_c, laser_ctrl_ids)
opto_hab_led_c_al   = align_by_animal(opto_hab_led_c,   led_ctrl_ids,   ctrl_animals)
opto_hab_laser_c_al = align_by_animal(opto_hab_laser_c, laser_ctrl_ids, ctrl_animals)

x_ticks = np.arange(1, 9)
hab_specs_c = [
    (fam_hab_avg_c,       'purple',  'Familiar',          'solid'),
    (opto_hab_led_c_al,   'violet',  'Fam+5mW',     'dashed'),
    (opto_hab_laser_c_al, '#e77295', 'Fam+10mW',  'dashed'),
    (nov_hab_avg_c,       'green',   'Novel',             'solid'),
]

fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
for data, color, label, ls in hab_specs_c:
    n_valid = np.sum(np.isfinite(data), axis=0).astype(float)
    n_valid[n_valid < 2] = np.nan
    mean = np.nanmean(data, axis=0)
    sem  = np.nanstd(data, axis=0, ddof=1) / np.sqrt(n_valid)
    ax.plot(x_ticks, mean, color=color, label=label, ls=ls, lw=1.5, ms=5)
    ax.errorbar(x_ticks, mean, yerr=sem, fmt='none', ecolor=color,
                elinewidth=1, capsize=3)

ax.axhline(0, color='black', lw=1, ls='--', alpha=0.6)
ax.set_xlabel("Presentation #", fontsize=13)
ax.set_ylabel("Δ sniffing rate (inh/s)", fontsize=13)
ax.set_title("Habituation — CONTROLS", pad=10, fontsize=13, weight='bold')
ax.set_xticks(x_ticks)
ax.spines[['right', 'top']].set_visible(False)
ax.legend()

if savefig:
    plt.savefig(f"{savefolder}/Habituation_CONTROLS.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Habituation_CONTROLS.svg", bbox_inches='tight')
plt.show()

#%% Direct contrast — opto (Fam+stim) Δ sniffing: EXPERIMENTAL vs CONTROL
# Unpaired (different animals): experimental cohort vs control cohort, per stim type.
# Uses experimental opto_led / opto_laser computed near the top of this script.

def _u_stars(a, b):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) >= 1 and len(b) >= 1 and (len(a) + len(b)) >= 3:
        return p_to_stars(stats.mannwhitneyu(a, b, alternative='two-sided')[1])
    return 'n/a'

contrast = [
    # (condition label, experimental vals, control vals, exp colour)
    ("Fam+5mW",   opto_led,   opto_led_c,   'violet'),
    ("Fam+10mW", opto_laser, opto_laser_c, '#e77295'),
]

fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
width = 0.34
xc = np.arange(len(contrast))

for k, (clabel, exp_vals, ctrl_vals, exp_c) in enumerate(contrast):
    for off, vals, c, lbl in [(-width / 2, exp_vals, exp_c, 'Experimental'),
                              (+width / 2, ctrl_vals, 'gray', 'Control')]:
        v = np.asarray(vals, float)
        v = v[np.isfinite(v)]
        rgba = mcolors.to_rgba(c)
        if len(v):
            m_val = v.mean()
            sem   = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
            ax.bar(xc[k] + off, m_val, width=width, color=(*rgba[:3], 0.30),
                   edgecolor=(0, 0, 0, 0.6), linewidth=1, zorder=0,
                   label=lbl if k == 0 else None)
            ax.errorbar(xc[k] + off, m_val, yerr=sem, capsize=2, color='black',
                        lw=1, zorder=2, alpha=0.7)
        # individual animals
        ax.scatter(np.full(len(v), xc[k] + off), v, color=c, s=22, zorder=3,
                   alpha=0.85, edgecolor='white', linewidth=0.4)

    # exp vs control significance per condition
    stars = _u_stars(np.asarray(exp_vals, float), np.asarray(ctrl_vals, float))
    allv = np.concatenate([np.asarray(exp_vals, float), np.asarray(ctrl_vals, float)])
    allv = allv[np.isfinite(allv)]
    if len(allv):
        ytop = allv.max() + 0.08 * (allv.max() - allv.min() + 1e-9)
        add_sig_bar(ax, xc[k] - width / 2, xc[k] + width / 2, ytop, 0.02, stars)

ax.axhline(0, color='black', lw=1, ls='--', alpha=0.5)
ax.set_xticks(xc)
ax.set_xticklabels([c[0] for c in contrast], fontsize=13)
ax.set_ylabel("Δ avg inhalations/sec", fontsize=14)
ax.set_title("Opto sniffing effect: experimental vs control", pad=15, fontsize=14, weight='bold')
ax.spines[['right', 'top']].set_visible(False)
ax.legend(frameon=False)

if savefig:
    plt.savefig(f"{savefolder}/Opto effect_experimental_vs_control.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{savefolder}/Opto effect_experimental_vs_control.svg", bbox_inches='tight')
plt.show()

#%% Smoothed + edge-trimmed time traces (correct for Gaussian boundary artifact)
# Additional versions of the Gaussian-smoothed sniffing time traces above.
# The per-trial traces are already convolved with a Gaussian; here we apply a
# further light Gaussian smooth to the across-mouse mean/SEM to clean up the
# line, then crop the first and last second — where the kernel runs off the
# edge of the recording window (_T_ARR spans -4 s to +8 s) and produces
# boundary artifacts.

_SMOOTH_S   = 0.30                                   # extra smoothing sigma (s)
_TRIM_S     = 1.0                                    # seconds trimmed each end
_smooth_sig = _SMOOTH_S * fps
_keep_mask  = (_T_ARR >= _T_ARR[0] + _TRIM_S) & (_T_ARR <= _T_ARR[-1] - _TRIM_S)


def plot_smoothed_trimmed(specs, title, fname):
    """Plot across-mouse mean±SEM traces with an extra Gaussian smooth applied,
    then cropped to remove the 1 s edge (smoothing) artifact on each side."""
    fig, ax = plt.subplots(figsize=(6, 3.5), dpi=300)
    t = _T_ARR[_keep_mask]
    for data, color, label, ls in specs:
        valid = np.all(np.isfinite(data), axis=1)
        d = data[valid]
        if not len(d):
            continue
        mean = d.mean(axis=0)
        sem  = d.std(axis=0, ddof=1) / np.sqrt(len(d)) if len(d) > 1 else np.zeros_like(mean)
        mean_s = gaussian_filter1d(mean, sigma=_smooth_sig)[_keep_mask]
        sem_s  = gaussian_filter1d(sem,  sigma=_smooth_sig)[_keep_mask]
        ax.plot(t, mean_s, color=color, label=label, ls=ls, lw=1.5)
        ax.fill_between(t, mean_s - sem_s, mean_s + sem_s, color=color, alpha=0.2)

    ax.axhline(0, color='black', lw=1, alpha=0.1, ls='dotted')
    ax.axvline(0, color='black', lw=1, alpha=0.7)
    ax.axvspan(0.3, 2.3, 0, 0.1, color='red', alpha=0.3, linewidth=0)
    ax.set_xlabel("Time from odor onset (s)", fontsize=13)
    ax.set_ylabel("Δ sniffing rate (Hz)", fontsize=13)
    ax.set_title(title, pad=10, fontsize=13, weight='bold')
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()

    if savefig:
        plt.savefig(f"{savefolder}/{fname}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{savefolder}/{fname}.svg", bbox_inches='tight')
    plt.show()


plot_smoothed_trimmed(
    trace_specs,
    "Mean sniffing (smoothed, edge-trimmed)",
    "Sniffing time trace_smoothed_LED_vs_laser",
)
plot_smoothed_trimmed(
    trace_specs_c,
    "Mean sniffing (smoothed, edge-trimmed) — CONTROLS",
    "Sniffing time trace_smoothed_CONTROLS",
)

# %%
