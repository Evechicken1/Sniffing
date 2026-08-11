"""
Diagnostic: first-sniff alignment vs. trial occurrence.

Hypothesis: for novel odour first presentations (occ=1) the first detected
inhalation onset is later than for repeated presentations, pushing ep_aligned
too far forward and making neural responses appear before t=0 in sniff-
aligned rasters.

Tests:
  1. ep_nframes distribution per occurrence, novel vs familiar (all sessions)
  2. Statistical test: occ=1 vs occ=2-8 (novel, all sessions pooled)
  3. Alignment variant comparison (current / option-A / option-B)
  4. Side-by-side raster reconstruction for neuron 27 (27_250118_KK152_123)
"""

import sys, os, pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

sys.path.append(r"C:\Users\xaand\Documents\PhD\Analysis\Informatics\Python\Electrophysiology\helpers")

DATA_PATH  = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\data"
SAVE_PATH  = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\Analysis\sniff_alignment_diagnostic\new_pp"
os.makedirs(SAVE_PATH, exist_ok=True)

with open(os.path.join(DATA_PATH, 'sniffs.pkl'), 'rb') as f:
    sniffs = pickle.load(f)
with open(os.path.join(DATA_PATH, 'spks.pkl'), 'rb') as f:
    spks = pickle.load(f)

n_ses        = len(sniffs)
pre_event    = spks['pre_event']       # 4 s
post_event   = spks['post_event']      # 6 s
pre_event_al = spks['pre_event_al']    # 2 s
post_event_al= spks['post_event_al']   # 4 s
ifr_sr       = spks['ifr_sr']          # 60 Hz
nsamp_diff   = int((pre_event - pre_event_al) * ifr_sr)   # 120
nsamples_al  = int((pre_event_al + post_event_al) * ifr_sr)  # 360
T_VEC_AL     = np.linspace(-pre_event_al, post_event_al, nsamples_al)

# ── TOL-vivid palette (same as plotting script) ──────────────────────────────
TOL_VIVID = ['#0077BB','#33BBEE','#009988','#EE7733','#CC3311','#EE3377','#BBBBBB']

# ═══════════════════════════════════════════════════════════════════════════
# Core alignment function — mirrors ep_preprocess.py exactly, but
# accepts parameter overrides so we can test variants.
# ═══════════════════════════════════════════════════════════════════════════

def compute_alignment(sniffs_ses, *, first_after=4.35, max_delay=1.5,
                      est_delay=0.35, cam_lag=-0.2):
    """
    Re-compute ep_aligned, ep_nframes, and tr_delays per trial.
    Returns tr_delays in seconds (= amount added to ep_onsets).
    """
    ntrials = len(sniffs_ses['trial_idx'])
    tr_delays = np.zeros(ntrials)

    for tr in range(ntrials):
        inh_sec = np.asarray(sniffs_ses['ml_inh_onsets'][tr], float) / 60
        inh_sec = inh_sec - first_after
        inh_postodor = inh_sec[inh_sec >= 0]

        if len(inh_postodor) > 0:
            tr_delay = inh_postodor[0]
        else:
            tr_delay = est_delay

        if tr_delay >= max_delay:
            tr_delay = est_delay
        else:
            tr_delay = tr_delay + est_delay

        tr_delay += cam_lag
        tr_delays[tr] = tr_delay

    ep_aligned_out = sniffs_ses['ep_onsets'] + tr_delays
    ep_nframes_out = np.round(tr_delays * ifr_sr).astype(int)
    return ep_aligned_out, ep_nframes_out, tr_delays


VARIANTS = {
    'current\n(max_delay=1.5)': dict(first_after=4.35, max_delay=1.5),
    'opt-A\n(max_delay=0.6)':   dict(first_after=4.35, max_delay=0.6),
    'opt-B\n(first_after=4.0)': dict(first_after=4.0,  max_delay=1.5),
}

# Pre-compute tr_delays for all sessions × variants
all_delays = {label: [] for label in VARIANTS}   # will hold per-trial tr_delay
all_occur  = []                                   # parallel arrays
all_nov    = []
all_ses_id = []

for ses in range(n_ses):
    occ = sniffs[ses]['trial_occur']
    nov = sniffs[ses]['trial_novelty'].astype(bool)
    all_occur.append(occ)
    all_nov.append(nov)
    all_ses_id.append(np.full(len(occ), ses, int))

    for label, params in VARIANTS.items():
        _, _, td = compute_alignment(sniffs[ses], **params)
        all_delays[label].append(td)

all_occur  = np.concatenate(all_occur)
all_nov    = np.concatenate(all_nov)
all_ses_id = np.concatenate(all_ses_id)
for label in VARIANTS:
    all_delays[label] = np.concatenate(all_delays[label])

# Compare against stored ep_nframes to verify current-variant matches original
stored_nframes = np.concatenate([sniffs[s]['ep_nframes'] for s in range(n_ses)])
recomp_nframes = np.round(all_delays['current\n(max_delay=1.5)'] * ifr_sr).astype(int)
n_mismatch = np.sum(stored_nframes != recomp_nframes)
print(f"[Verification] recomputed vs stored ep_nframes mismatches: {n_mismatch} / {len(stored_nframes)}")

# ═══════════════════════════════════════════════════════════════════════════
# Plot 1 — tr_delay by occurrence × novelty (current alignment)
# ═══════════════════════════════════════════════════════════════════════════

max_occ = 8
occ_bins = np.arange(1, max_occ + 1)
td_current = all_delays['current\n(max_delay=1.5)']

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, (cond_label, cond_mask) in zip(axes, [('Novel', all_nov), ('Familiar', ~all_nov)]):
    data_by_occ = [td_current[(all_occur == o) & cond_mask] for o in occ_bins]
    bp = ax.boxplot(data_by_occ, positions=occ_bins, widths=0.6, patch_artist=True,
                    flierprops=dict(marker='.', markersize=2, alpha=0.3),
                    medianprops=dict(color='k', linewidth=1.5))
    for patch in bp['boxes']:
        patch.set_facecolor('#AED6F1')
    ax.set_xlabel('Trial occurrence', fontsize=11)
    ax.set_title(cond_label, fontsize=12, fontweight='bold')
    ax.set_xticks(occ_bins)
    ax.axhline(0, color='gray', lw=0.5, ls='--')
    ax.tick_params(direction='in')

axes[0].set_ylabel('tr_delay (s from poke to sniff ref)', fontsize=10)
fig.suptitle('Sniff alignment delay per trial occurrence — current parameters', fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(SAVE_PATH, 'delay_by_occurrence_current.png'), dpi=200)
plt.show()

# ═══════════════════════════════════════════════════════════════════════════
# Plot 2 — same but per session (novel only)
# ═══════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(3, 4, figsize=(14, 9), sharey=True, sharex=True)
for ses, ax in enumerate(axes.flat):
    ses_mask = (all_ses_id == ses) & all_nov
    data_by_occ = [td_current[(all_occur == o) & ses_mask] for o in occ_bins]
    bp = ax.boxplot(data_by_occ, positions=occ_bins, widths=0.6, patch_artist=True,
                    flierprops=dict(marker='.', markersize=1.5, alpha=0.3),
                    medianprops=dict(color='k', linewidth=1.2))
    for patch in bp['boxes']:
        patch.set_facecolor('#AED6F1')
    ax.axhline(0, color='gray', lw=0.5, ls='--')
    ax.set_title(sniffs[ses].get('folder_identifier', f'ses {ses}'), fontsize=7)
    ax.tick_params(direction='in', labelsize=7)

fig.supxlabel('Trial occurrence', fontsize=10)
fig.supylabel('tr_delay (s)', fontsize=10)
fig.suptitle('Sniff alignment delay — novel trials per session', fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(SAVE_PATH, 'delay_by_occurrence_per_session.png'), dpi=200)
plt.show()

# ═══════════════════════════════════════════════════════════════════════════
# Test 2 — statistical test: occ=1 vs occ=2–8 (novel, all sessions pooled)
# ═══════════════════════════════════════════════════════════════════════════

nov_td = td_current[all_nov]
nov_occ = all_occur[all_nov]

groups = [nov_td[nov_occ == o] for o in occ_bins]
H, p_kw = stats.kruskal(*groups)
print(f"\n[Stats] Novel: Kruskal-Wallis across occurrences 1–{max_occ}: H={H:.2f}, p={p_kw:.4g}")

occ1  = nov_td[nov_occ == 1]
occ2p = nov_td[nov_occ > 1]
U, p_mwu = stats.mannwhitneyu(occ1, occ2p, alternative='two-sided')
med1  = np.median(occ1)
med2p = np.median(occ2p)
print(f"[Stats] Novel occ=1 median delay: {med1*1000:.0f} ms  |  occ≥2: {med2p*1000:.0f} ms")
print(f"[Stats] Mann-Whitney U occ=1 vs occ≥2: U={U:.0f}, p={p_mwu:.4g}")

# Per-session occ=1 vs occ>1 median comparison (for reporting)
print("\n[Per-session] Novel occ=1 vs occ>1 median tr_delay (ms):")
print(f"  {'session':<20}  occ=1   occ>1   diff")
for ses in range(n_ses):
    ses_nov = (all_ses_id == ses) & all_nov
    d1  = td_current[ses_nov & (all_occur == 1)]
    d2p = td_current[ses_nov & (all_occur > 1)]
    if len(d1) and len(d2p):
        m1, m2p = np.median(d1)*1000, np.median(d2p)*1000
        folder = sniffs[ses].get('folder_identifier', f'ses {ses}')
        print(f"  {folder:<20}  {m1:5.0f}   {m2p:5.0f}   {m1-m2p:+.0f}")

# ═══════════════════════════════════════════════════════════════════════════
# Plot 3 — variant comparison (novel occ=1 vs occ>1 per variant)
# ═══════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
for ax, (label, td) in zip(axes, {k: all_delays[k] for k in VARIANTS}.items()):
    nov_td_v = td[all_nov]
    occ1_v  = nov_td_v[nov_occ == 1]
    occ2p_v = nov_td_v[nov_occ > 1]
    ax.boxplot([occ1_v, occ2p_v], labels=['occ=1', 'occ≥2'],
               patch_artist=True,
               flierprops=dict(marker='.', markersize=2, alpha=0.3),
               medianprops=dict(color='k', linewidth=1.5),
               boxprops=dict(facecolor='#AED6F1'))
    U_v, p_v = stats.mannwhitneyu(occ1_v, occ2p_v, alternative='two-sided')
    ax.set_title(f'{label}\np={p_v:.3g}', fontsize=9)
    ax.axhline(0, color='gray', lw=0.5, ls='--')
    ax.tick_params(direction='in')

axes[0].set_ylabel('tr_delay (s)', fontsize=10)
fig.suptitle('Alignment delay occ=1 vs occ≥2 — novel trials, parameter variants', fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(SAVE_PATH, 'variant_comparison_occ1_vs_occ2p.png'), dpi=200)
plt.show()

# ═══════════════════════════════════════════════════════════════════════════
# Plot 4 — Neuron 27 raster: current vs each variant
# ═══════════════════════════════════════════════════════════════════════════

NEURONS = [27, 65, 66, 71]   # 27_250118_KK152_123  |  65_250118_KK152_274

# ═══════════════════════════════════════════════════════════════════════════
# Helpers for raster reconstruction
# ═══════════════════════════════════════════════════════════════════════════

def recompute_raster(nrn, orig_ep_aligned, ep_aligned_new):
    """Shift cntrd_ts_al by the change in alignment reference per trial."""
    delta = ep_aligned_new - orig_ep_aligned
    raster_out = []
    for i in range(len(delta)):
        spk = np.asarray(spks['cntrd_ts_al'][nrn][i], float) - delta[i]
        spk = spk[(spk >= -pre_event_al) & (spk <= post_event_al)]
        raster_out.append(spk)
    return raster_out

def recompute_fr_al(nrn, ep_nframes_new):
    """Re-slice all_fr using updated ep_nframes."""
    fr_poke = spks['all_fr'][nrn]
    n_trials = fr_poke.shape[0]
    out = np.full((n_trials, nsamples_al), np.nan)
    for i in range(n_trials):
        tf = nsamp_diff + ep_nframes_new[i]
        if tf >= 0 and tf + nsamples_al <= fr_poke.shape[1]:
            out[i, :] = fr_poke[i, tf : tf + nsamples_al]
    return out

def draw_novel_col(ax_r, ax_f, trial_chem, trial_idx, nov_mask, raster, fr_al):
    ymax = 0
    for a, odor in enumerate(np.unique(trial_chem[nov_mask])):
        odor_mask = (trial_chem == odor) & nov_mask
        idxs      = trial_idx[odor_mask]
        spk_list  = [raster[i] for i in idxs]
        ypos      = np.arange(ymax, ymax + len(spk_list))
        color     = TOL_VIVID[a % len(TOL_VIVID)]
        ax_r.eventplot(spk_list, color=color, lineoffsets=ypos, linewidths=0.5)
        ax_f.plot(T_VEC_AL, fr_al[odor_mask].mean(0), color=color, linewidth=1)
        ymax += len(spk_list)
    for ax in (ax_r, ax_f):
        ax.axvline(0, color='gray', lw=0.5)
        ax.tick_params(direction='in', labelsize=7)
    ax_r.set_xlim(-pre_event_al, post_event_al)

# Pre-compute variant alignments per session (avoids redundant loops when
# two neurons share the same session)
variant_ep_al = {}   # (ses, label) -> ep_aligned_new
variant_ep_nf = {}   # (ses, label) -> ep_nframes_new
for ses_v in set(int(spks['ses_idx'][n]) for n in NEURONS):
    for label, params in VARIANTS.items():
        ep_al, ep_nf, _ = compute_alignment(sniffs[ses_v], **params)
        variant_ep_al[(ses_v, label)] = ep_al
        variant_ep_nf[(ses_v, label)] = ep_nf

n_cols = len(VARIANTS)   # current is included as the first variant

for NRN in NEURONS:
    ses     = int(spks['ses_idx'][NRN])
    unit_id = f"{NRN}_{spks['ses_id'][NRN]}_{spks['nrn_id'][NRN]}"
    print(f"\n[Raster] Reconstructing raster for {unit_id} (ses_idx={ses})")

    orig_ep_aligned = sniffs[ses]['ep_aligned']
    trial_idx_s     = sniffs[ses]['trial_idx'] - 1
    trial_nov_s     = sniffs[ses]['trial_novelty'].astype(bool)
    trial_chem_s    = sniffs[ses]['trial_chem_id']

    fig, axes = plt.subplots(2, n_cols, figsize=(4 * n_cols, 6),
                             sharey='row', sharex='all')
    fig.suptitle(f'{unit_id}  —  Novel condition\nTop: raster  |  Bottom: mean FR per odour',
                 fontsize=9)

    for col, (label, params) in enumerate(VARIANTS.items()):
        ep_al_new = variant_ep_al[(ses, label)]
        ep_nf_new = variant_ep_nf[(ses, label)]

        if 'current' in label:
            raster = [spks['cntrd_ts_al'][NRN][i] for i in range(len(trial_idx_s))]
            fr_al  = np.array(spks['all_fr_al'][NRN])
        else:
            raster = recompute_raster(NRN, orig_ep_aligned, ep_al_new)
            fr_al  = recompute_fr_al(NRN, ep_nf_new)

        draw_novel_col(axes[0, col], axes[1, col],
                       trial_chem_s, trial_idx_s, trial_nov_s, raster, fr_al)
        axes[0, col].set_title(label.replace('\n', ' '), fontsize=8, fontweight='bold')
        if col == 0:
            axes[0, col].set_ylabel('Trial #', fontsize=8)
            axes[1, col].set_ylabel('Spikes / sec', fontsize=8)
        axes[1, col].set_xlabel('Time [sec]', fontsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(SAVE_PATH, f'raster_comparison_{unit_id}.png'), dpi=200)
    plt.show()

# ═══════════════════════════════════════════════════════════════════════════
# Plot 5 — Per-occurrence PSTH (chronological)
#
# For each neuron: all novel trials in chronological order, showing raster
# and IFR.  Also overlays mean PSTH per occurrence (occ 1–8) so you can see
# whether response latency shifts across presentations.
# ═══════════════════════════════════════════════════════════════════════════

MAX_OCC_PSTH = 8

def plot_chron_occurrence_psth(NRN, raster, fr_al, unit_id):
    """
    Three-panel figure for one neuron:
      Row 0: raster — ALL trials in chronological session order, coloured by
             odour; novel trials full opacity, familiar trials faded
      Row 1: IFR heatmap, sniff-aligned (same chronological ordering)
      Row 2: individual trial IFR poke-aligned (all_fr) for novel trials,
             shaded dark→light by occurrence (occ=1 darkest)
    """
    ses          = int(spks['ses_idx'][NRN])
    trial_idx_s  = sniffs[ses]['trial_idx'] - 1   # already in session order
    trial_nov_s  = sniffs[ses]['trial_novelty'].astype(bool)
    trial_chem_s = sniffs[ses]['trial_chem_id']
    trial_occ_s  = sniffs[ses]['trial_occur']

    n_rows       = len(trial_idx_s)
    unique_chems = np.unique(trial_chem_s)
    chem_to_color = {c: TOL_VIVID[i % len(TOL_VIVID)]
                     for i, c in enumerate(unique_chems)}

    row_colors = [chem_to_color[trial_chem_s[i]] for i in range(n_rows)]
    row_alpha  = [1.0 if trial_nov_s[i] else 0.25 for i in range(n_rows)]

    # Sniff-aligned IFR in chronological order
    fr_chron = fr_al[trial_idx_s]   # (n_rows, nsamples_al)

    # Poke-aligned firing rate and time vector
    fr_poke      = spks['all_fr'][NRN]
    nsamples_poke = fr_poke.shape[1]
    T_VEC_POKE   = np.linspace(-pre_event, post_event, nsamples_poke)

    # Rows 0–1 share sniff-aligned x; Row 2 has independent poke-aligned x
    fig = plt.figure(figsize=(9, 10))
    gs  = fig.add_gridspec(3, 1, height_ratios=[2, 1.2, 1.2], hspace=0.35)
    ax_r = fig.add_subplot(gs[0])
    ax_h = fig.add_subplot(gs[1], sharex=ax_r)
    ax_p = fig.add_subplot(gs[2])   # independent x-axis

    fig.suptitle(f'{unit_id}  —  All trials, chronological order\n'
                 f'Rows 0–1: sniff-aligned  |  Row 2: poke-aligned, '
                 f'novel only, dark = occ 1',
                 fontsize=9)

    # -- Row 0: raster --
    for row in range(n_rows):
        ax_r.eventplot([raster[trial_idx_s[row]]],
                       colors=[row_colors[row]],
                       lineoffsets=[row],
                       linewidths=0.5,
                       alpha=row_alpha[row])
    ax_r.axvline(0, color='gray', lw=0.6)
    ax_r.set_ylabel('Trial # (chronological)', fontsize=8)
    ax_r.tick_params(direction='in', labelsize=7)
    ax_r.set_ylim(-0.5, n_rows - 0.5)

    # -- Row 1: IFR heatmap --
    vmax = np.nanpercentile(fr_chron, 98)
    im = ax_h.imshow(fr_chron, aspect='auto', origin='lower',
                     extent=[-pre_event_al, post_event_al, 0, n_rows],
                     cmap='viridis', vmin=0, vmax=vmax)
    ax_h.axvline(0, color='white', lw=0.6)
    ax_h.set_ylabel('Trial # (chronological)', fontsize=8)
    ax_h.set_xlabel('Time re. sniff [sec]', fontsize=8)
    ax_h.tick_params(direction='in', labelsize=7)
    plt.colorbar(im, ax=ax_h, label='FR (Hz)', pad=0.01, fraction=0.03)

    # -- Row 2: individual trial IFR, poke-aligned, novel only --
    # Colour = odour identity (same palette as raster); alpha dims with occurrence
    def _occ_alpha(occ, max_occ=MAX_OCC_PSTH):
        return 0.9 - 0.6 * (occ - 1) / max(max_occ - 1, 1)

    legend_handles = {}   # chem -> first Line2D for legend
    for row in range(n_rows):
        if not trial_nov_s[row]:
            continue
        ti   = trial_idx_s[row]
        occ  = int(trial_occ_s[row])
        chem = trial_chem_s[row]
        if occ > MAX_OCC_PSTH:
            continue
        color = chem_to_color[chem]
        alpha = _occ_alpha(occ)
        line, = ax_p.plot(T_VEC_POKE, fr_poke[ti], color=color,
                          linewidth=0.7, alpha=alpha)
        if chem not in legend_handles:
            legend_handles[chem] = line
    ax_p.axvline(0, color='gray', lw=0.6)
    ax_p.set_xlabel('Time re. poke [sec]', fontsize=8)
    ax_p.set_ylabel('IFR — novel (Hz)', fontsize=8)
    ax_p.tick_params(direction='in', labelsize=7)
    if legend_handles:
        ax_p.legend(list(legend_handles.values()),
                    [f'odour {c}' for c in legend_handles],
                    fontsize=6, ncol=4, loc='upper right', framealpha=0.6)

    fname = os.path.join(SAVE_PATH, f'chron_occ_psth_{unit_id}.png')
    fig.savefig(fname, dpi=200)
    plt.show()
    print(f"  Saved: {fname}")


for NRN in NEURONS:
    ses     = int(spks['ses_idx'][NRN])
    unit_id = f"{NRN}_{spks['ses_id'][NRN]}_{spks['nrn_id'][NRN]}"
    print(f"\n[Chron PSTH] {unit_id}")

    # Use current (stored) alignment for now; easy to swap in a variant later
    raster_cur = [spks['cntrd_ts_al'][NRN][i]
                  for i in range(len(sniffs[ses]['trial_idx']))]
    fr_al_cur  = np.array(spks['all_fr_al'][NRN])

    plot_chron_occurrence_psth(NRN, raster_cur, fr_al_cur, unit_id)

# ═══════════════════════════════════════════════════════════════════════════
# Plot 6 — Per-odorant × presentation grid
#
# Columns = novel odorants  |  Rows = presentation number (occ 1…MAX_OCC_PSTH)
# Bottom summary row = mean ± SEM across all presentations.
# All panels share the same x- and y-axes for direct comparison.
# Saved to per_neuron_psth/ subfolder.
# ═══════════════════════════════════════════════════════════════════════════

PSTH_SAVE_PATH = os.path.join(SAVE_PATH, 'per_neuron_psth')
os.makedirs(PSTH_SAVE_PATH, exist_ok=True)

def plot_per_odorant_by_occ(NRN, unit_id):
    ses          = int(spks['ses_idx'][NRN])
    trial_idx_s  = sniffs[ses]['trial_idx'] - 1
    trial_nov_s  = sniffs[ses]['trial_novelty'].astype(bool)
    trial_chem_s = sniffs[ses]['trial_chem_id']
    trial_occ_s  = sniffs[ses]['trial_occur']

    # Poke-aligned (odour presentation = t=0)
    fr_al         = np.array(spks['all_fr'][NRN])
    n_samp        = fr_al.shape[1]
    t_vec         = np.linspace(-pre_event, post_event, n_samp)

    nov_chems     = np.unique(trial_chem_s[trial_nov_s])
    n_chems       = len(nov_chems)
    chem_to_color = {c: TOL_VIVID[i % len(TOL_VIVID)]
                     for i, c in enumerate(np.unique(trial_chem_s))}

    # MAX_OCC_PSTH data rows + 1 mean summary row
    n_rows = MAX_OCC_PSTH + 1
    fig, axes = plt.subplots(n_rows, n_chems,
                             figsize=(2.8 * n_chems, 1.6 * n_rows),
                             sharex=True, sharey=True,
                             gridspec_kw={'hspace': 0.05, 'wspace': 0.05})
    # Ensure axes is always 2-D
    if n_chems == 1:
        axes = axes[:, np.newaxis]

    fig.suptitle(unit_id, fontsize=10, fontweight='bold')

    for col, chem in enumerate(nov_chems):
        color = chem_to_color[chem]
        axes[0, col].set_title(f'odour {chem}', fontsize=8,
                               color=color, fontweight='bold')

        all_traces = []
        for occ_idx, occ in enumerate(range(1, MAX_OCC_PSTH + 1)):
            ax = axes[occ_idx, col]
            mask = trial_nov_s & (trial_chem_s == chem) & (trial_occ_s == occ)
            idxs = trial_idx_s[mask]

            ax.axvline(0, color='gray', lw=0.5, ls='--', zorder=0)
            ax.spines[['top', 'right']].set_visible(False)
            show_yticks = (col == 0)
            ax.tick_params(direction='in', labelsize=6,
                           left=show_yticks, labelleft=show_yticks,
                           bottom=False, labelbottom=False)

            if len(idxs):
                trace = fr_al[idxs[0]]
                ax.plot(t_vec, trace, color=color, linewidth=0.9)
                all_traces.append(trace)

            if col == 0:
                ax.set_ylabel(f'{occ}', fontsize=7, rotation=0,
                              ha='right', va='center', labelpad=10)

        # Mean ± SEM summary row
        ax_m = axes[MAX_OCC_PSTH, col]
        if all_traces:
            stacked = np.array(all_traces)
            mfr     = np.nanmean(stacked, axis=0)
            sfr     = np.nanstd(stacked, axis=0) / np.sqrt(stacked.shape[0])
            ax_m.plot(t_vec, mfr, color=color, linewidth=1.4)
            ax_m.fill_between(t_vec, mfr - sfr, mfr + sfr,
                              color=color, alpha=0.2)
        ax_m.axvline(0, color='gray', lw=0.5, ls='--', zorder=0)
        ax_m.spines[['top', 'right']].set_visible(False)
        ax_m.tick_params(direction='in', labelsize=6,
                         left=(col == 0), labelleft=(col == 0))
        ax_m.set_xlabel('t re. odour [sec]', fontsize=7)
        if col == 0:
            ax_m.set_ylabel('mean', fontsize=7, rotation=0,
                            ha='right', va='center', labelpad=10)

    # Shared y-label for the occurrence rows
    fig.text(0.005, 0.55, 'Presentation #', va='center',
             rotation='vertical', fontsize=9)

    fname = os.path.join(PSTH_SAVE_PATH, f'per_odor_occ_{unit_id}.png')
    fig.savefig(fname, dpi=200, bbox_inches='tight')
    plt.show()
    print(f"  Saved: {fname}")


for NRN in NEURONS:
    unit_id = f"{NRN}_{spks['ses_id'][NRN]}_{spks['nrn_id'][NRN]}"
    print(f"\n[Per-odour × occ grid] {unit_id}")
    plot_per_odorant_by_occ(NRN, unit_id)

print(f"\nAll diagnostic plots saved to {SAVE_PATH}")
