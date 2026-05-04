

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
                '260324_XT013','260324_XT014','260324_XT015','260324_XT016']
# all LED: ['260302_XT009','260302_XT010','260302_XT011','260302_XT012',
#            '260324_XT013','260324_XT014','260324_XT015','260324_XT016']

# Laser (10 mW) sessions – on target
laser_sess_ids = ['260306_XT009','260306_XT010','260306_XT011','260306_XT012',\
                  '260328_XT013','260328_XT014','260328_XT015','260328_XT016']
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
t_range  = np.array([4.3, 7.6]) * fps
bl_range = np.array([0.0, 4.0]) * fps
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

graph      = ["Familiar", "Fam+5mW\n(LED)", "Fam+10mW\n(Laser)", "Novel"]
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
    (opto_tr_led,   'violet',     'Fam+5mW (LED)',    'dashed'),
    (opto_tr_laser, "#e77295",   'Fam+10mW (Laser)', 'dashed'),
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
    (opto_hab_led,   'violet',   'Fam+5mW (LED)',    'dashed'),
    (opto_hab_laser, '#e77295', 'Fam+10mW (Laser)', 'dashed'),
    (nov_hab_avg,    'green',    'Novel',            'solid'),
]

fig, ax = plt.subplots(figsize=(6, 4), dpi=300)

for data, color, label, ls in hab_specs:
    n_valid = np.sum(np.isfinite(data), axis=0).astype(float)
    n_valid[n_valid < 2] = np.nan
    mean = np.nanmean(data, axis=0)
    sem  = np.nanstd(data, axis=0, ddof=1) / np.sqrt(n_valid)
    ax.plot(x_ticks, mean, color=color, label=label, ls=ls, lw=1.5, marker='o', ms=5)
    ax.fill_between(x_ticks, mean - sem, mean + sem, color=color, alpha=0.2, linewidth=0)

ax.axhline(0, color='black', lw=1, ls='--', alpha=0.6)
ax.set_xlabel("Presentation #", fontsize=13)
ax.set_ylabel("Δ sniffing rate (inh/s)", fontsize=13)
ax.set_title("Habituation", pad=10, fontsize=13, weight='bold')
ax.set_xticks(x_ticks)
ax.spines[['right', 'top']].set_visible(False)
ax.legend()

if savefig:
    plt.savefig(f"{savefolder}/Habituation_LED_vs_laser.png", dpi=300, bbox_inches='tight')
plt.show()

# %%
