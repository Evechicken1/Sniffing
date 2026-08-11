#%%
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\xaand\Documents\PhD\Analysis\Informatics\Python\Sniffing')
import sniff_tools as st
from scipy import stats
import datetime
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

#### Import sniffs dictionary ####

# OFC ephys dataset: KK152 / KK153 / KK154, 250118-250121.
# All 12 sessions are non-opto novel/familiar only, so we use the clean
# (non-opto) novel-vs-familiar pipeline ported from ldt_ep_sniffing.py.
#
# Loaded straight from the per-session postprocessed .mat files (same as
# ldt_ep_sniffing), NOT the pre-bundled sniffs.pkl.

animals = ['KK152', 'KK153', 'KK154']
dates = ['250118', '250119', '250120', '250121']
sess_ids = [f"{d}_{a}" for a in animals for d in dates]

paths_list = [f"C:\\Users\\xaand\\Documents\\PhD\\Experiments\\postprocessed_files\\{sess_id}" for sess_id in sess_ids]

sniffs = st.import_sniff_mat_select(paths_list)

if len(sniffs) > 0:
    print('Succesfully imported', len(sniffs), 'mice', [sniffs[m]['folder_identifier'] for m in range(len(sniffs))])


#%% Structure data into dataframes
##Variables to adjust

trialcorr = 0 #applies baseline correction per odor exposure trial; 1 applies baseline correction per mouse(day)
savefig = 1 #save figures to folder
savefolder = r"C:\Users\xaand\Documents\PhD\Experiments\Ephys OFC\Analysis\Sniffing\old_preprocesssing"
if savefig and not os.path.exists(savefolder):
    os.makedirs(savefolder)

# ---- Optional plot toggles (set to True to generate the plot) ----
PLOT_ABS_SNIFF_PER_ODOR        = False  # Plotting absolute sniffing responses per odor
PLOT_RESP_DISTRIBUTION         = False  # Plotting distribution of odor-responses per condition
PLOT_BREATHING_CHANGE_PER_ODOR = False  # Breathing change per odor per condition
PLOT_INDIV_TRIALS              = False  # Individual trials overlayed on average trace (use only if habituation curves are flat)
PLOT_BASELINE_OVER_PRESENT     = False  # Baseline over presentations

# short label for filenames (full sess_ids list is not a valid filename)
sess_label = "_".join(animals) + f"_n{len(sess_ids)}"

nframes = 713
fps = 713/12

#%% Structuring data for usage
#initialize some variables to not error
isnov=0
isfam=0
nmice = len(sniffs)
baseline_bins = range(0,4)
n_presentations = np.flatnonzero((c := np.bincount(sniffs[0]['trial_occur'])) == c.max())[-1] # takes into account that blank odor can be presented extra times

#non-opto (novel/familiar) variables
nconds = 2
conditions = [isnov, isfam]
colors = ['green', 'purple']
ecolor = ['darkgreen', 'purple']
graph = ['novel', 'familiar']
markercolors = ['purple', 'green']
BoI = range(5, 8)
seq = [1, 0]
comp_grps = ["nov vs fam"]
line_type = ["solid", "solid"]
ngraph = len(graph)

# Defining of the selector (ntrials x conditions) !! important to accurately predefine conditions before selecting!!
ntrials = sniffs[0]["trial_idx"].size
selector = [np.zeros((ntrials, len(conditions)), int) for i in range(nmice)]

for m in range(nmice):
    isnov = (sniffs[m]["trial_novelty"] == True) & (sniffs[m]["trial_opto"] == False) & (sniffs[m]["trial_occur"] <= 3)
    isfam = (sniffs[m]["trial_familiarity"] == True) & (sniffs[m]["trial_opto"] == False) & (sniffs[m]["trial_occur"] <= 3)

    conditions = [isnov, isfam]

    trial_idx = sniffs[m]["trial_idx"] - 1

    for cond in range(selector[m].shape[1]):
        idx_cnds = trial_idx[conditions[cond]]
        selector[m][idx_cnds, cond] = True

nsel = selector[0][:, 0].sum()
ngraph = len(conditions)

#%% Create a list of arrays for each condition, for plotting the different conditions
sniff_lst = [[[] for m in range(nmice)] for g in range(len(graph))]
hist_edges = np.linspace(0, 723, 13)

baseline_lst = [[[] for n in range(nmice)] for g in range(len(graph))]

for g in range(len(graph)):
    for m in range(len(selector)):
        tmp = sniffs[m]["ml_inh_onsets"][selector[m][:, g] == 1]
        for i in range(len(tmp)):
            hist = np.histogram(tmp[i], bins=hist_edges)[0]
            sniff_lst[g][m].append(hist)
            baseline_lst[g][m].append(np.mean(sniff_lst[g][m][i][baseline_bins]))

sniff_arr = np.array(sniff_lst)

#%% Creating a full data array (nmice x nconditions x npresentations x nodors x nbins)
m_data = []
n_frames = 720
n_bins = 12
bin_size = n_frames/n_bins
bins = np.linspace(0,720,13)

for m in range(nmice):
    trial_idx = sniffs[m]['trial_idx'] - 1
    odors = np.unique(sniffs[m]['trial_chem_id'])
    nov_data = []
    fam_data = []
    tmp = []

    nov_trials = [(sniffs[m]['trial_chem_id'] == odors[o]) & (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0) for o in range(len(odors))]
    nov_trials = [i for i in nov_trials if sum(i)>0]
    fam_trials = [(sniffs[m]['trial_chem_id'] == odors[o]) & (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0) for o in range(len(odors))]
    fam_trials = [i for i in fam_trials if sum(i)>0]

    isnov = len(nov_trials) >0
    isfam = len(fam_trials) >0

    if isnov:
        tmp = [[] for n in range(n_presentations)]
        for n in range(n_presentations):
            for o in range(len(nov_trials)):
                tmp[n].append(np.histogram(sniffs[m]['ml_inh_onsets'][trial_idx[nov_trials[o]]][n], bins = bins)[0])
        nov_data = np.array(tmp)

    if isfam:
        tmp = [[] for n in range(n_presentations)]
        for n in range(n_presentations):
            for o in range(len(fam_trials)):
                tmp[n].append(np.histogram(sniffs[m]['ml_inh_onsets'][trial_idx[fam_trials[o]]][n], bins = bins)[0])
        fam_data = np.array(tmp)

    conditions = [nov_data, fam_data]
    m_data.append(conditions)
m_data = np.array(m_data, dtype=object)

n_odors = m_data[0, 0].shape[1]

#%% Plotting average sniffs for all mice, during the trial
bin_edges = np.linspace(-3.5, 7.5, num = 12)

plt.figure(figsize=(6,3.5))
for g in range(ngraph):
    mean = sniff_arr[g].mean(axis = (0))
    baseline = sniff_arr[g,:,:,baseline_bins].mean(axis = (0,1))
    mean = mean - baseline[:,None]
    mean = mean.mean(axis = 0)
    sem = np.std(sniff_arr[g].mean(axis = 1), axis = 0)

    plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
    plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
plt.title("Mean sniffing (first 3 presentations)", pad = 10, fontsize=15, weight = 'bold')

plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
plt.xlabel("time from odor presentation (s)", fontsize=15)
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.legend()

if savefig:
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{savefolder}/mean sniffing to first 3 presentations_sniffs_{sess_label}_{now}.png"
    plt.savefig(filename, dpi=400, bbox_inches="tight")

plt.show()

#%% Plotting average inhalations per mouse, throughout the trial
bin_edges = np.linspace(-3.5, 7.5, num = 12)
for m in range(nmice):
    plt.figure()
    plt.title(f"Mean sniffing first 3 presentations ({sess_ids[m]})")
    plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
    plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
    plt.ylabel("Δ avg inhalations (inh/s)")
    plt.xlabel("time from odor presentation (s)")
    ax = plt.gca()
    ax.spines[['right', 'top']].set_visible(False)

    for g in range(len(graph)):
        mean = sniff_arr[g,m]
        baseline = sniff_arr[g,m,:,baseline_bins].mean(axis = 0)
        mean = mean - baseline[:,None]
        mean = mean.mean(axis = 0)
        sem = np.std(sniff_arr[g].mean(axis = 1), axis = 0)

        plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
        plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
    plt.legend()

    if savefig:
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        plt.savefig(f"{savefolder}/Sniffing to first 3 presentations ({sess_ids[m]})_{now}.png", dpi=400, bbox_inches="tight")
    plt.show()


#%% Comparison of sniffing during specific timebin
############### Calculating data

#t_range = np.array((4.3,6.6))*fps
t_range = np.array((4.3,7.6))*fps
t_time = (t_range[1]-t_range[0])/fps
bl_range = np.array((0,4))*fps
bl_time = (bl_range[1]-bl_range[0])/fps

nov = []
fam = []

for m in range(nmice):

    trial_idx = sniffs[m]["trial_idx"] - 1
    nov_odors = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    fam_odors = (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)

    nov_idxs = trial_idx[nov_odors]
    fam_idxs = trial_idx[fam_odors]

    nov_data_tmp = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][nov_idxs]])/t_time
    fam_data_tmp = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][fam_idxs]])/t_time

    nov_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][nov_idxs]])/bl_time
    fam_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][fam_idxs]])/bl_time

    nov_data = (nov_data_tmp - nov_bl).mean()
    fam_data = (fam_data_tmp - fam_bl).mean()

    nov.append(nov_data)
    fam.append(fam_data)

############### Plotting
# ---- helper functions ----
def p_to_stars(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 5e-2 else 'ns'

def add_sig_bar(ax, x1, x2, y, h, stars):
    xi1, xi2 = xpos[x1], xpos[x2]
    ax.plot([xi1, xi1, xi2, xi2],
            [y, y+h, y+h, y],
            lw=1, c='k',
            clip_on=False, zorder=5)
    ax.text((xi1+xi2)/2, y + h*1.05, stars,
            ha='center', va='bottom',
            clip_on=False, zorder=6)

# ================= PLOT =================
############### Order adjustment
# put familiar first, then novel
graph = ["Familiar", "Novel"]

# reorder rows to match graph order
data = np.array([fam, nov])

############### Plotting
xpos = {label: i for i, label in enumerate(graph)}
nconds = len(graph)
fig_width = max(4, nconds * 1.6)
plt.figure(figsize=(fig_width, 5), dpi=300)
seq_arr = np.zeros(shape=(nmice, nconds))

for m in range(nmice):
    avg = data[:, m]
    seq_arr[m, :] = avg

    # light per-mouse trace
    plt.plot(graph, seq_arr[m], color='gray', linestyle='-', alpha=0.9)

    # markers per group
    plt.plot(graph[0], seq_arr[m, 0], color='purple', marker='o', linestyle='-', alpha = 0.8)  # familiar
    plt.plot(graph[1], seq_arr[m, 1], color='green',  marker='o', linestyle='-', alpha = 0.8)  # novel

plt.ylabel("Δ avg inhalations/sec",  fontsize=17)
plt.title("Mean sniffing change", pad = 40,  fontsize=17, weight = 'bold')
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)

# colors
bar_colors = ['purple', 'green']

# group means and SEM across mice
means = np.nanmean(seq_arr, axis=0)
ns    = np.sum(np.isfinite(seq_arr), axis=0).astype(float)
sems  = np.nanstd(seq_arr, axis=0, ddof=1) / np.sqrt(np.maximum(ns, 1))

# bar plot
bar_width = 0.6
for i, (x, mean, sem, c) in enumerate(zip(graph, means, sems, bar_colors)):
   base_rgba = mcolors.to_rgba(c)
   base_rgba_black = mcolors.to_rgba('black')

   face_rgba = (base_rgba[0], base_rgba[1], base_rgba[2], 0.25)
   edge_rgba = (base_rgba_black[0], base_rgba_black[1], base_rgba_black[2], 0.6)

   plt.bar(i, mean, width=bar_width,
           color=face_rgba, edgecolor=edge_rgba,
           linewidth=1, zorder=0)

   plt.errorbar(i, mean, yerr=sem, capsize=2,
                color=base_rgba_black, lw=1, zorder=2, alpha=0.7)

# ---- significance bars (pairwise, paired tests across mice) ----
n_groups = len(graph)
pairs = []
if n_groups >= 2:
    pairs.append((0,1))

ymin = np.nanmin(seq_arr)
ymax = np.nanmax(seq_arr)
yrng = (ymax - ymin) if np.isfinite(ymax - ymin) and (ymax - ymin) > 0 else 1.0
base_y = ymax + 0.02 * yrng
step_h = 0.08 * yrng
line_h = 0.01 * yrng
margin = 0.1 * (ymax - ymin if ymax > ymin else 1)

raw_ps = []
tests = []
for (i, j) in pairs:
    a = seq_arr[:, i]
    b = seq_arr[:, j]
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() >= 2:
        t, p = stats.wilcoxon(a[mask], b[mask])
    else:
        p = np.nan
    raw_ps.append(p)
    tests.append((i, j, p))

# Bonferroni correction
m = max(1, len([p for p in raw_ps if np.isfinite(p)]))
for k, (i, j, p) in enumerate(tests):
    if not np.isfinite(p):
        stars = 'n/a'
    else:
        p_adj = min(p * m, 1.0)
        stars = p_to_stars(p_adj)

    y = base_y + k * step_h
    add_sig_bar(ax, graph[i], graph[j], y, line_h, stars)

ax.set_ylim(bottom=ymin - 0.05*yrng, top=base_y + margin)

plt.xticks(np.array(graph),  fontsize=17)

if savefig:
    plt.savefig(f"{savefolder}/Mean sniffing change_{sess_label}.png", dpi=300, bbox_inches="tight")

plt.show()


#%% Plotting breathing change per mouse throughout trial per condition for first 3 presentations
bin_edges = np.linspace(-3.5, 7.5, num = 12)

which_presentation = [(range(0,3)),(range(0,3))] #first 3 presentations
plt.figure(figsize=(6,3.5))
for m in range(nmice):
    for g in range(ngraph):
        data = m_data[m,g,which_presentation[g]]
        baseline = m_data[m,g,which_presentation[g],:][:,:,baseline_bins].mean(axis = (-1))
        mean = data -baseline[:,:,None]
        mean = mean.mean(axis = (0,1))
        sem = np.std(data[:,g].mean(axis = 1).astype(float), axis = 0)/np.sqrt(nmice)

        plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
        plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
    plt.title(f"Mean sniffing first 3 presentations {sniffs[m]['folder_identifier']}")

    plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
    plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
    plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
    plt.xlabel("time from odor presentation (s)", fontsize=15)
    ax = plt.gca()
    ax.spines[['right', 'top']].set_visible(False)
    plt.legend()

    if savefig:
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{savefolder}/Sniffing to first 3 presentations_sniffs_{sess_ids[m]}.png"
        plt.savefig(filename, dpi=400, bbox_inches="tight")
    plt.show()

#%% Plotting breathing change all mice first 3 presentations
bin_edges = np.linspace(-3.5, 7.5, num = 12)

which_presentation = [(range(0,3)),(range(0,3))] #first 3 presentations
plt.figure(figsize=(6,3.5))
for g in range(ngraph):
    data = m_data[:,g,which_presentation[g]]
    baseline = m_data[:,g,which_presentation[g],:][:,:,:,baseline_bins].mean(axis = (-1))
    mean_tmp = data -baseline[:,:,:,None]
    mean = mean_tmp.mean(axis = (0,1,2))
    sem = np.std(data[:,g].mean(axis = 1).astype(float), axis=0)/np.sqrt(nmice)

    plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
    plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
plt.title("Mean sniffing (first 3 presentations)", pad = 10, fontsize=15, weight = 'bold')

plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
plt.xlabel("time from odor presentation (s)", fontsize=15)
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.legend()

if savefig:
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = rf"{savefolder}\Sniffing to first 3 presentations_sniffs_{sess_label}.png"
    plt.savefig(filename, dpi=400, bbox_inches="tight")

plt.show()

#%% Plotting average habituation across mice per condition
plt.figure(figsize=(5,5))
m_data.shape # (n_mice, n_conditions, n_presentations, n_odors, n_bins)
for i in range(ngraph):
    x_ticks = np.arange(n_presentations)
    data = m_data[:,i].mean(axis=0) # select condition i across mice
    baseline = data[:,:,baseline_bins].mean(axis=(2))
    mean = data[:,:,BoI].mean(axis=(2))
    mean = np.array(mean, dtype = float)
    sem = np.std(mean, axis=1) / np.sqrt(mean.shape[0])

    corr_mean = mean - baseline
    corr_mean = corr_mean.mean(axis=1)
    corr_mean = np.array(corr_mean, dtype = float)

    plt.plot(x_ticks, corr_mean, color = colors[i], label=graph[i])
    plt.fill_between(x_ticks, corr_mean +sem , corr_mean - sem, color = ecolor[i], alpha = 0.5, linewidth = 0)
plt.axhline(0, color='black', linestyle='--')
plt.title('Habituation')
plt.xlabel("Presentation #")
plt.ylabel("Inhalations/second (baseline-subtracted)")
plt.legend()
if savefig:
    hab_dir = savefolder + '/habituation_curves'
    if not os.path.exists(hab_dir):
        os.makedirs(hab_dir)
    plt.savefig(rf"{hab_dir}/{sess_label}.png", dpi = 300)
plt.show()

#%% Plotting average habituation per mouse per condition
m_data.shape # (n_mice, n_conditions, n_presentations, n_odors, n_bins)
for m in range(nmice):

    plt.figure(figsize=(5,5))
    for i in range(ngraph):

        x_ticks = np.arange(n_presentations)
        data = m_data[m,i] # select condition i for mouse m
        baseline = data[:, :, baseline_bins].mean(axis=(2))
        mean = data[:, :, BoI].mean(axis=(2))
        mean = np.array(mean, dtype = float)
        sem = np.std(mean, axis=1) / np.sqrt(mean.shape[0])

        corr_mean = mean - baseline
        corr_mean = corr_mean.mean(axis=1)
        corr_mean = np.array(corr_mean, dtype = float)

        plt.plot(x_ticks, corr_mean, color = colors[i], label=graph[i], ls = line_type[i])
        plt.fill_between(x_ticks, corr_mean +sem , corr_mean - sem, color = ecolor[i], alpha = 0.3, linewidth = 0, ls=line_type[i])
    plt.axhline(0, color='black', linestyle='--', alpha = 0.6, lw = .9)
    plt.title(rf'Habituation {sess_ids[m]}', pad = 10, fontsize=13)
    plt.xlabel("Presentation #")
    plt.ylabel("Inhalations/second (baseline-subtracted)")
    plt.legend()
    if savefig:
        hab_dir = savefolder + '/habituation_curves'
        if not os.path.exists(hab_dir):
            os.makedirs(hab_dir)
        plt.savefig(rf"{hab_dir}/{sess_ids[m]}.png", dpi = 300)
    plt.show()


#%% Plotting breathing change throughout trial per condition for first 3 novel presentations and all familiar presentations
#  - only use if warranted by habituation curves (flat)
bin_edges = np.linspace(-3.5, 7.5, num = 12)

which_presentation = [(range(0,3)),(range(0,n_presentations))] #first 3 presentations for novel, all presentations for familiar
plt.figure(figsize=(6,3.5))
for g in range(ngraph):
    data = m_data[:,g,which_presentation[g]]
    baseline = m_data[:,g,which_presentation[g],:][:,:,:,baseline_bins].mean(axis = (-1))
    mean = data -baseline[:,:,:,None]
    mean = mean.mean(axis = (0,1,2))
    sem = np.std(sniff_arr[g].mean(axis = 1), axis = 0)/np.sqrt(int(nmice))

    plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
    plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
plt.title("Mean sniffing (first 3 / all presentations)", pad = 10, fontsize=15, weight = 'bold')

plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
plt.xlabel("time from odor presentation (s)", fontsize=15)
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.legend()

if savefig:
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = rf"{savefolder}\Sniffing to first 3_all presentations_sniffs_{sess_label}.png"
    plt.savefig(filename, dpi=400, bbox_inches="tight")

plt.show()


# %% Plotting absolute sniffing responses per odor
if PLOT_ABS_SNIFF_PER_ODOR:
    t_range = np.array((4.6,6.6))*fps
    t_time = (t_range[1]-t_range[0])/fps
    bl_range = np.array((.5,3.5))*fps
    bl_time = (bl_range[1]-bl_range[0])/fps

    for m in range(nmice):

        trial_idx = sniffs[m]["trial_idx"] - 1

        is_nov = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
        is_fam = (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
        conditions = [is_nov, is_fam]
        cond_names = ["novel", "familiar"]
        cond_colors = ['green', 'purple']
        fig, ax = plt.subplots(1, len(conditions), figsize=(10, 5))
        for c in range(len(conditions)):
            odors = np.unique(sniffs[m]['trial_chem_id'][conditions[c]])
            odor_ls = []
            odor_sem_ls = []
            for o in odors:
                cond_odors = (sniffs[m]['trial_chem_id'] == o) & (conditions[c])
                odor_idxs = trial_idx[cond_odors]
                odor_data = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][odor_idxs]])/t_time
                odor_mean = np.nanmean(odor_data)
                odor_sem = np.nanstd(odor_data)/np.sqrt(len(odor_data))

                odor_ls.append(odor_mean)
                odor_sem_ls.append(odor_sem)

            data = odor_ls
            xticks = np.arange(0, len(data))
            ax[c].bar(xticks, data, color=cond_colors[c])
            ax[c].errorbar(xticks, data, yerr=odor_sem_ls, fmt='o', color='black', ecolor='black', elinewidth=1, capsize=3)
            ax[c].set_title(f"{cond_names[c]}", fontweight='bold')
            ax[c].set_xticks(xticks, labels=odors)
            ax[c].set_xlabel("Odor")
            ax[c].set_ylabel("inhalations/sec")
            ax[c].set_ylim(0,10)

            fig.suptitle(f"Absolute breathing per odor {sniffs[m]['folder_identifier']}", fontsize=16, fontweight='bold')
        plt.show()

#%% Plotting distribution of odor-responses per condition
if PLOT_RESP_DISTRIBUTION:
    for m in range(nmice):
        plt.figure()

        trial_idx = sniffs[m]["trial_idx"] - 1

        is_nov = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
        is_fam = (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
        conditions = [is_nov, is_fam]
        cond_names = ["novel", "familiar"]
        cond_colors = ['green', 'purple']
        fig, ax = plt.subplots(1, len(conditions), figsize=(10, 2))
        for c in range(len(conditions)):
            odors = np.unique(sniffs[m]['trial_chem_id'][conditions[c]])
            print(cond_names[c], odors)
            odor_ls = []
            for o in odors:
                cond_odors = (sniffs[m]['trial_chem_id'] == o) & (conditions[c])
                odor_idxs = trial_idx[cond_odors]
                odor_data = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][odor_idxs]])/t_time
                odor_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][odor_idxs]])/bl_time
                odor_data = np.nanmean(odor_data - odor_bl)
                odor_ls.append(odor_data)

            data = odor_ls
            ax[c].hist(data, bins=6)
            ax[c].set_title(f"{cond_names[c]}", fontweight='bold')
            ax[c].set_xlabel("Odor")
            ax[c].set_xlim(-2,6)
            ax[c].set_ylabel("Δ avg inhalations/sec")
            ax[c].axvline(x=0, color="black",linewidth=1, alpha = 0.7)

            fig.suptitle(f"Breathing change per odor {sniffs[m]['animal_name']}", fontsize=16, fontweight='bold')

#%% Breathing change per odor per condition
if PLOT_BREATHING_CHANGE_PER_ODOR:
    t_range = np.array((4.6,6.6))*fps
    t_time = (t_range[1]-t_range[0])/fps
    bl_range = np.array((.5,3.5))*fps
    bl_time = (bl_range[1]-bl_range[0])/fps

    for m in range(nmice):

        trial_idx = sniffs[m]["trial_idx"] - 1

        is_nov = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
        is_fam = (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
        conditions = [is_nov, is_fam]
        cond_names = ["novel", "familiar"]
        cond_colors = ['green', 'purple']
        fig, ax = plt.subplots(1, len(conditions), figsize=(10, 5))
        for c in range(len(conditions)):
            odors = np.unique(sniffs[m]['trial_chem_id'][conditions[c]])
            odor_ls = []
            odor_sem_ls = []
            for o in odors:
                cond_odors = (sniffs[m]['trial_chem_id'] == o) & (conditions[c])
                odor_idxs = trial_idx[cond_odors]
                odor_data = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][odor_idxs]])/t_time
                odor_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][odor_idxs]])/bl_time
                odor_mean = np.nanmean(odor_data - odor_bl)
                odor_sem = np.nanstd(odor_data - odor_bl)/np.sqrt(len(odor_data))

                odor_ls.append(odor_mean)
                odor_sem_ls.append(odor_sem)

            data = odor_ls
            xticks = np.arange(0, len(data))
            ax[c].bar(xticks, data, color=cond_colors[c])
            ax[c].errorbar(xticks, data, yerr=odor_sem_ls, fmt='o', color='black', ecolor='black', elinewidth=1, capsize=3)
            ax[c].set_title(f"{cond_names[c]}", fontweight='bold')
            ax[c].set_xticks(xticks, labels=odors)
            ax[c].set_xlabel("Odor")
            ax[c].set_ylabel("Δ avg inhalations/sec")
            ax[c].set_ylim(-2,6)

            fig.suptitle(f"Breathing change per odor {sniffs[m]['animal_name']}", fontsize=16, fontweight='bold')
            if savefig:
                plt.savefig(f"{savefolder}/Breathing_change_per_odor_{sess_ids[m]}.png", dpi=300, bbox_inches="tight")


#%% Plot individual trials overlayed on average trace per mouse - only use if warranted by habituation curves (flat)
if PLOT_INDIV_TRIALS:
    bin_edges = np.linspace(-3.5, 7.5, num=12)

    for m in range(nmice):

        fig, axes = plt.subplots(1, ngraph, figsize=(8,4), sharey=True)

        for g in range(len(graph)):

            ax = axes[g]

            # ---- pull trials ----
            trials = sniff_arr[g, m, 6:]                        # (n_trials, n_bins)
            baseline = trials[:, baseline_bins].mean(axis=1, keepdims=True)
            trials_bl = trials - baseline

            mean = trials_bl.mean(axis=0)
            sem  = trials_bl.std(axis=0, ddof=1) / np.sqrt(trials_bl.shape[0])

            ax.set_title(graph[g])
            ax.axhline(y=0, color="black", linewidth=1, alpha=0.1, ls="dotted")
            ax.axvline(x=0, color="black", linewidth=1, alpha=0.7)

            ax.set_xlabel("time from odor presentation (s)")
            if g == 0:
                ax.set_ylabel("Δ avg inhalations (inh/s)")

            ax.spines[['right','top']].set_visible(False)

            # ---- plot individual trials ----
            for tr in range(trials_bl.shape[0]):
                ax.plot(
                    bin_edges,
                    trials_bl[tr],
                    color=colors[g],
                    alpha=0.2,
                    linewidth=1,
                    ls=line_type[g]
                )

            # ---- overlay mean ± SEM ----
            ax.plot(bin_edges, mean, color=colors[g], linewidth=2, ls=line_type[g])
            ax.errorbar(
                bin_edges, mean, yerr=sem, fmt='o',
                color=colors[g], ecolor=ecolor[g], elinewidth=1, capsize=3
            )

        fig.suptitle(f"Sniffing individual trials ({sess_ids[m]})")
        plt.tight_layout()

        if savefig:
            now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            plt.savefig(
                f"{savefolder}/Sniffing_trials_subplot_individual_trial_traces_{sess_ids[m]}_{now}.png",
                dpi=400, bbox_inches="tight"
            )

        plt.show()


#%% Baseline over presentations
if PLOT_BASELINE_OVER_PRESENT:
    save_dir = savefolder + '/Baseline habituation'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    baseline_BoI = baseline_bins

    # ---- make subplots ----
    n_cols = 1
    n_rows = nmice

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows), sharex=True, sharey=True)
    axes = np.array(axes).reshape(-1)

    x_ticks = np.arange(n_presentations)
    mean_ls = []

    for m in range(nmice):
        ax = axes[m]

        # mean during baseline across conditions, odors and baseline bins
        baseline_data = m_data[m][:,:,:,baseline_BoI].mean(axis = (0,2,3))

        mean = np.array(baseline_data, dtype=float)
        mean_ls.append(mean)

        ax.plot(x_ticks, mean, color="black", markersize=5, label=sess_ids[m])

        ax.set_title(sess_ids[m])
        ax.set_xlabel("Presentation")
        ax.set_ylabel('Inhalations/second')
    plt.suptitle('Baseline change', fontweight="bold", y = .91)
    plt.show()

    plt.figure(figsize=(5,5))
    plt.plot(np.array(mean_ls).mean(axis=0), color = 'black')
    plt.title('Baseline change across mice', fontweight = 'bold')
    plt.xlabel('Presentation #')
    plt.ylabel('Inhalations/second')
    plt.show()

# %%
