#%%
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

#!!! Session ids should  be of one type (stim, inhibition or non-opto)

##odor sessions - sessions from which to take odors paired with stim

#%% Get odors from previous sessions

#odor_sessions = ['260324_XT013','260324_XT014','260324_XT015', '260324_XT016'] 
odor_sessions = ['260328_XT013','260328_XT014','260328_XT015', '260328_XT016']
#odor_sessions = ['260328_XT013','260328_XT014', '260328_XT016']



odor_paths_list = [f"Y:\\postprocessed_data\\{sess_id}" for sess_id in odor_sessions]
odor_sniffs = st.import_sniff_mat_select(odor_paths_list)
odors_opto = [np.unique(odor_sniffs[m]['trial_chem_id'][odor_sniffs[m]['trial_opto'] ==1]) for m in range(len(odor_sniffs))]
odors_nov = [np.unique(odor_sniffs[m]['trial_chem_id'][odor_sniffs[m]['trial_novelty'] ==1]) for m in range(len(odor_sniffs))]


#%%
##LED 

#all 
#sess_ids = ['260303_XT009','260303_XT010','260303_XT011', '260303_XT012'] 
#sess_ids = ['260325_XT013','260325_XT014','260325_XT015', '260325_XT016']


#on target
#sess_ids = ['260302_XT009','260302_XT010','260302_XT011']

##laser 

#sess_ids = ['260307_XT009','260307_XT010','260307_XT011', '260307_XT012']
sess_ids = ['260329_XT013','260329_XT014','260329_XT015', '260329_XT016']

# check that sessions are from consecutive days
for sess_id in range(len(sess_ids)):
    sess_date = sess_ids[sess_id].split("_")[0]
    odor_date = odor_sessions[sess_id].split("_")[0]
    if int(sess_date) != int(odor_date)+1:
        print(f"Warning: session {sess_ids[sess_id]} is not from the day after its corresponding odor session {odor_sessions[sess_id]}")


paths_list = [f"Y:\\postprocessed_data\\{sess_id}" for sess_id in sess_ids]

sniffs = st.import_sniff_mat_select(paths_list)

if len(sniffs)>0:
    print('Succesfully imported', len(sniffs), 'mice',[sniffs[m]['folder_identifier'] for m in range(len(sniffs))])


#%% Structure data into dataframes
##Variables to adjust

savefig = 1 #save figures to folder
savefolder = r"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\Sniffing to nov fam"

nframes = 713
fps = 713/12

#%% Structuring data for usage
#initialize some variables to not error
isnov=0
isfam=0
isopto=0
isstim = any(sniffs[0]["trial_opto"]==1)
isinh = any(sniffs[0]["trial_opto"]==-1)
nmice = len(sniffs)
baseline_bins = range(0,4)
n_presentations = np.flatnonzero((c := np.bincount(sniffs[0]['trial_occur'])) == c.max())[-1] # takes into account that blank odor can be presented 4 times extra
ngraph = 3

nconds = 3
conditions = [isnov,isfam, isopto]    
colors = ['green' , 'purple', 'violet'] 
ecolor = ["darkgreen","purple", 'violet']
graph = ["novel","familiar", "fam+stim"]
markercolors = ['purple','green','violet'] 
BoI = range(5, 8)
seq = [1,0,2]
comp_grps= ["nov vs fam", "nov vs opto", "fam vs opto"]
line_type = ["solid", "solid","dashed"]

#%% Creating a full data array
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
    opto_data = []

    tmp = []

    nov_trials = [(sniffs[m]['trial_chem_id'] == odors[o]) & (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)for o in range(len(odors))]
    nov_trials = [i for i in nov_trials if sum(i)>0]
    fam_trials = [(sniffs[m]['trial_chem_id'] == odors[o]) & (~np.isin(sniffs[m]['trial_chem_id'], odors_opto[m])) & (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0) for o in range(len(odors))]
    fam_trials = [i for i in fam_trials if sum(i)>0]
    opto_trials = [(sniffs[m]['trial_chem_id'] == odors_opto[m][o]) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_blank'] == 0) for o in range(len(odors_opto[m]))]
    opto_trials = [i for i in opto_trials if sum(i)>0]
    print(opto_trials)

    isnov = len(nov_trials) >0
    isfam = len(fam_trials) >0
    isopto = len(opto_trials) >0
        
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
        
    if isopto:
        tmp = [[] for n in range(n_presentations)]
        for n in range(n_presentations):
            for o in range(len(opto_trials)):
                tmp[n].append(np.histogram(sniffs[m]['ml_inh_onsets'][trial_idx[opto_trials[o]]][n], bins = bins)[0])
        opto_data = np.array(tmp)

    conditions = [nov_data, fam_data, opto_data]
    m_data.append(conditions)
m_data = np.array(m_data, dtype=object)


#%% Comparison of sniffing during specific timebin
############### Calculating data

#t_range = np.array((4.6,6.6))*fps
t_range = np.array((4.3,7.6))*fps
t_time = (t_range[1]-t_range[0])/fps
bl_range = np.array((0,4))*fps
bl_time = (bl_range[1]-bl_range[0])/fps
graph = ["familiar", "fam+stim","novel"]

nov = []
fam = []
opto = []


for m in range(nmice):
    
    trial_idx = sniffs[m]["trial_idx"] - 1 
    #nov_odors = np.logical_and(sniffs[m]['trial_novelty'] == 1,  sniffs[m]['trial_occur'] <= 3, sniffs[m]['trial_opto'] == 0, sniffs[m]['trial_blank'] == 0)
    nov_odors = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    fam_odors = (sniffs[m]['trial_familiarity'] == 1) & (~np.isin(sniffs[m]['trial_chem_id'], odors_opto[m])) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    opto_odors = np.isin(sniffs[m]['trial_chem_id'], odors_opto[m]) # & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_blank'] == 0)
    blank_opto_odors = (sniffs[m]['trial_blank'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0)

    nov_idxs = trial_idx[nov_odors]
    fam_idxs = trial_idx[fam_odors]
    opto_idxs = trial_idx[opto_odors]
    blank_idxs = trial_idx[blank_opto_odors]
    
    nov_data_tmp = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][nov_idxs]])/t_time
    fam_data_tmp = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][fam_idxs]])/t_time
    opto_data_tmp = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][opto_idxs]])/t_time


    nov_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][nov_idxs]])/bl_time
    fam_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][fam_idxs]])/bl_time
    opto_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][opto_idxs]])/bl_time

    nov_data = (nov_data_tmp - nov_bl).mean()
    fam_data = (fam_data_tmp - fam_bl).mean()
    opto_data = (opto_data_tmp - opto_bl).mean()

    nov.append(nov_data)
    fam.append(fam_data)
    opto.append(opto_data)

data = np.array([nov,fam,opto])

############### Plotting
# ---- helper functions ----
# make a mapping: string -> integer position

def p_to_stars(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 5e-2 else 'ns'

def add_sig_bar(ax, x1, x2, y, h, stars):
    xi1, xi2 = xpos[x1], xpos[x2]
    ax.plot([xi1, xi1, xi2, xi2],
            [y, y+h, y+h, y],
            lw=1, c='k',
            clip_on=False, zorder=5)        # allow drawing outside
    ax.text((xi1+xi2)/2, y + h*1.05, stars,
            ha='center', va='bottom',
            clip_on=False, zorder=6)
    
# ================= PLOT =================
############### Order adjustment
# put familiar first, then fam+stim (opto), then novel
graph = ["Familiar", "Familiar\n+ stim", "Novel"]

# data is still [nov, fam, opto] → reorder rows
data = np.array([fam, opto, nov])   # now matches graph order

############### Plotting
xpos = {label: i for i, label in enumerate(graph)}
nconds = len(graph)
fig_width = max(4, nconds * 1.6)
plt.figure(figsize=(fig_width, 5), dpi=300)
seq_arr = np.zeros(shape=(nmice, nconds))

for m in range(nmice):
    avg = data[:, m]
    seq_arr[m, :] = avg  # already ordered

    # light per-mouse trace
    plt.plot(graph, seq_arr[m], color='gray', linestyle='-', alpha=0.9)

    # markers per group
    plt.plot(graph[0], seq_arr[m, 0], color='purple', marker='o', linestyle='-', alpha = 0.8)  # familiar
    plt.plot(graph[1], seq_arr[m, 1], color='purple', marker='o', linestyle='-', alpha = 0.8)  # fam+stim
    plt.plot(graph[2], seq_arr[m, 2], color='green',  marker='o', linestyle='-', alpha = 0.8)  # novel

plt.ylabel("Δ avg inhalations/sec",  fontsize=17)
plt.title("Mean sniffing change", pad = 40,  fontsize=17, weight = 'bold')
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)

# colors
bar_colors = ['purple', 'purple', 'green']

# group means and SEM across mice
means = np.nanmean(seq_arr, axis=0)
ns    = np.sum(np.isfinite(seq_arr), axis=0).astype(float)
sems  = np.nanstd(seq_arr, axis=0, ddof=1) / np.sqrt(np.maximum(ns, 1))



# bar plot
bar_width = 0.6
bar_width = 0.6
for i, (x, mean, sem, c) in enumerate(zip(graph, means, sems, bar_colors)):
    # convert to RGBA
   base_rgba = mcolors.to_rgba(c)
   base_rgba_black = mcolors.to_rgba('black')

   # make two versions: translucent fill, stronger edge
   face_rgba = (base_rgba[0], base_rgba[1], base_rgba[2], 0.25)  # alpha for fill
   edge_rgba = (base_rgba_black[0], base_rgba_black[1], base_rgba_black[2], 0.6)  # alpha for edge

   # draw normal bar
   plt.bar(i, mean, width=bar_width,
           color=face_rgba, edgecolor=edge_rgba,
           linewidth=1, zorder=0)

   # add error bar
   plt.errorbar(i, mean, yerr=sem, capsize=2,
                color=base_rgba_black, lw=1, zorder=2, alpha=0.7)
# ---- significance bars (pairwise, paired t-tests across mice) ----
# decide which pairs to test, based on how many groups you actually plotted
n_groups = len(graph)
pairs = []
if n_groups >= 2:
    pairs.append((0,1))
if n_groups >= 3:
    pairs.extend([(0,2), (1,2)])

# compute y positions
ymin = np.nanmin(seq_arr)
ymax = np.nanmax(seq_arr)
yrng = (ymax - ymin) if np.isfinite(ymax - ymin) and (ymax - ymin) > 0 else 1.0
0.8
# start just above the tallest point
base_y = ymax + 0.02 * yrng
step_h = 0.08 * yrng  # vertical spacing between bars
line_h = 0.01 * yrng  # height of the little bracket
margin = 0.1 * (ymax - ymin if ymax > ymin else 1)

# run tests and draw
m_valid_tests = 0
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

    # stagger bars upward
    y = base_y + k * step_h
    add_sig_bar(ax, graph[i], graph[j], y, line_h, stars)

# make sure everything fits
ax.set_ylim(bottom=ymin - 0.05*yrng, top=base_y + margin)

plt.xticks(np.array(graph),  fontsize=17)

if savefig:
    plt.savefig(f"{savefolder}/Mean sniffing change_{sess_ids}.png", dpi=300, bbox_inches="tight")

plt.show()


#%% Plotting breathing change per mouse throughout trial per condition for first 3 presentations
bin_edges = np.linspace(-3.5, 7.5, num = 12)

which_presentation = [(range(0,3)),(range(0,3)),(range(0,3))] #first 3 presentations
plt.figure(figsize=(6,3.5))
for m in range(nmice):
    for g in range(ngraph):
        data = m_data[m,g,which_presentation[g]]
        baseline = m_data[m,g,which_presentation[g],:][:,:,baseline_bins].mean(axis = (-1))
        mean = data -baseline[:,:,None]
        mean = mean.mean(axis = (0,1))
        sem = np.std(data[:,g].mean(axis = 1), axis = 0)/np.sqrt(nmice)
            
        plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
        plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
    plt.title(f"Mean sniffing first 3 presentations_post-stim day {sniffs[m]['folder_identifier']}")

    plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
    plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
    plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
    plt.xlabel("time from odor presentation (s)", fontsize=15)
    ax = plt.gca() #to remove the top and right spines
    ax.spines[['right', 'top']].set_visible(False) #to remove the top and right spines
    plt.legend()

    if any(sniffs[0]["trial_opto"]==1): #shading for stim
        plt.axvspan(xmin=0.6, xmax=2.6, color="red", alpha=0.3, linewidth = 0)
    elif any(sniffs[0]["trial_opto"]==-1): #shading for inh
        x_values = [0] + [3 + i * 0.1 for i in range(10)]
        alpha_values = [0.33] + [0.3 - i * 0.03 for i in range(10)]
        [plt.axvspan(x_values[i], x_values[i + 1], color="blue", alpha=alpha_values[i], linewidth=0) for i in range(len(x_values) - 1)]
        
    if savefig:
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  # e.g. 20250819_184512
        filename = f"{savefolder}/Sniffing to first 3 presentations_sniffs_post-stim day{sess_ids}.png"
        plt.savefig(filename, dpi=400, bbox_inches="tight")
    plt.show()

#%% Plotting breathing change all mice first 3 presentations

bin_edges = np.linspace(-3.5, 7.5, num = 12)

sem = [[[] for m in range(nmice)] for l in range(ngraph)]
mean = [[[]for m in range(nmice)] for l in range(ngraph)]
baseline = [[[]for m in range(nmice)] for l in range(ngraph)]

which_presentation = [(range(0,3)),(range(0,3)),(range(0,3))] #first 3 presentations for novel, all presentations for familiar and opto
plt.figure(figsize=(6,3.5))
for g in range(ngraph):
    data = m_data[:,g,which_presentation[g]]
    baseline = m_data[:,g,which_presentation[g],:][:,:,:,baseline_bins].mean(axis = (-1))
    mean_tmp = data -baseline[:,:,:,None]
    mean = mean_tmp.mean(axis = (0,1,2))
    sem = np.std(data[:,g].mean(axis = 1).astype(float), axis=0)/np.sqrt(nmice)
        
    plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
    plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
#plt.title(f"Mean sniffing first 3 presentations {sniffs[0]['session_date']}")
plt.title("Mean sniffing (first 3 presentations)", pad = 10, fontsize=15, weight = 'bold')

plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
plt.xlabel("time from odor presentation (s)", fontsize=15)
ax = plt.gca() #to remove the top and right spines
ax.spines[['right', 'top']].set_visible(False) #to remove the top and right spines
plt.legend()

if any(sniffs[0]["trial_opto"]==1): #shading for stim
    plt.axvspan(xmin=0.6, xmax=2.6, color="red", alpha=0.3, linewidth = 0)
elif any(sniffs[0]["trial_opto"]==-1): #shading for inh
    x_values = [0] + [3 + i * 0.1 for i in range(10)]
    alpha_values = [0.33] + [0.3 - i * 0.03 for i in range(10)]
    [plt.axvspan(x_values[i], x_values[i + 1], color="blue", alpha=alpha_values[i], linewidth=0) for i in range(len(x_values) - 1)]
    
if savefig:
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  # e.g. 20250819_184512
    filename = f"{savefolder}\Sniffing to first 3 presentations_sniffs_{sess_ids}.png"
    plt.savefig(filename, dpi=400, bbox_inches="tight")

plt.show()

#%% Plotting average habituation per mouse per condition

plt.figure(figsize=(5,5))
m_data.shape # (n_mice, n_conditions, n_presentations, n_odors, n_bins)
for i in range(len(conditions)):
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
plt.savefig(rf"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\habituation_curves\{sess_ids}.png", dpi = 300)

plt.show()

#%% Plotting average habituation across mice per condition
m_data.shape # (n_mice, n_conditions, n_presentations, n_odors, n_bins)
for m in range(nmice):
        
    plt.figure(figsize=(5,5))
    for i in range(len(conditions)):

        x_ticks = np.arange(n_presentations)
        data = m_data[m,i] # select condition i across mice
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
    plt.savefig(rf"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\habituation_curves\{sess_ids[m]}.png", dpi = 300)

    plt.show()


#%% Plotting breathing change throughout trial per condition for first 3 novel presentation and all familiar/opto presentations
#  - only use if warranted by habituation curves (flat)

bin_edges = np.linspace(-3.5, 7.5, num = 12)

sem = [[[] for m in range(nmice)] for l in range(ngraph)]
mean = [[[]for m in range(nmice)] for l in range(ngraph)]
baseline = [[[]for m in range(nmice)] for l in range(ngraph)]

which_presentation = [(range(0,3)),(range(0,3)),(range(0,3))] #first 4 presentations for novel, all presentations for familiar and opto
plt.figure(figsize=(6,3.5))
for g in range(ngraph):
    data = m_data[:,g,which_presentation[g]]
    baseline = m_data[:,g,which_presentation[g],:][:,:,:,baseline_bins].mean(axis = (-1))
    mean = data -baseline[:,:,:,None]
    mean = mean.mean(axis = (0,1,2))
    sem = np.std(m_data[:,g].mean(axis = 1), axis = 0)/np.sqrt(int(nmice))
        
    plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
    plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
#plt.title(f"Mean sniffing first 3 presentations {sniffs[0]['session_date']}")
plt.title("Mean sniffing (first 3/8/8 presentations)", pad = 10, fontsize=15, weight = 'bold')

plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
plt.xlabel("time from odor presentation (s)", fontsize=15)
ax = plt.gca() #to remove the top and right spines
ax.spines[['right', 'top']].set_visible(False) #to remove the top and right spines
plt.legend()

if any(sniffs[0]["trial_opto"]==1): #shading for stim
    plt.axvspan(xmin=0.6, xmax=2.6, color="red", alpha=0.3, linewidth = 0)
elif any(sniffs[0]["trial_opto"]==-1): #shading for inh
    x_values = [0] + [3 + i * 0.1 for i in range(10)]
    alpha_values = [0.33] + [0.3 - i * 0.03 for i in range(10)]
    [plt.axvspan(x_values[i], x_values[i + 1], color="blue", alpha=alpha_values[i], linewidth=0) for i in range(len(x_values) - 1)]
    
if savefig:
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  # e.g. 20250819_184512
    filename = f"{savefolder}\Sniffing to first 3_3_3 presentations_sniffs_{sess_ids}.png"
    plt.savefig(filename, dpi=400, bbox_inches="tight")

plt.show()


#%% Plotting breathing change per mouse throughout trial per condition for first 3 presentations for novel and 8/8 for familiar and opto
# - only use if warranted by habituation curves (flat)

bin_edges = np.linspace(-3.5, 7.5, num = 12)

which_presentation = [(range(0,4)),(range(0,n_presentations)),(range(0,n_presentations))] #first 4 presentations for novel, all presentations for familiar and opto
plt.figure(figsize=(6,3.5))
for m in range(nmice):
    for g in range(ngraph):
        data = m_data[m,g,which_presentation[g]]
        baseline = m_data[m,g,which_presentation[g],:][:,:,baseline_bins].mean(axis = (-1))
        mean = data -baseline[:,:,None]
        mean = mean.mean(axis = (0,1))
        sem = np.std(sniff_arr[g].mean(axis = 1), axis = 0)/np.sqrt(nmice)
            
        plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
        plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
    plt.title(f"Mean sniffing first 3/8/8 presentations {sniffs[m]['folder_identifier']}")

    plt.axhline(y=0, color="black",linewidth=1, alpha = 0.1, ls = "dotted")
    plt.axvline(x=0, color="black",linewidth=1, alpha = 0.7)
    plt.ylabel("Δ avg inhalations (inh/s)", fontsize=15)
    plt.xlabel("time from odor presentation (s)", fontsize=15)
    ax = plt.gca() #to remove the top and right spines
    ax.spines[['right', 'top']].set_visible(False) #to remove the top and right spines
    plt.legend()

    if any(sniffs[0]["trial_opto"]==1): #shading for stim
        plt.axvspan(xmin=0.6, xmax=2.6, color="red", alpha=0.3, linewidth = 0)
    elif any(sniffs[0]["trial_opto"]==-1): #shading for inh
        x_values = [0] + [3 + i * 0.1 for i in range(10)]
        alpha_values = [0.33] + [0.3 - i * 0.03 for i in range(10)]
        [plt.axvspan(x_values[i], x_values[i + 1], color="blue", alpha=alpha_values[i], linewidth=0) for i in range(len(x_values) - 1)]
        
    if savefig:
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  # e.g. 20250819_184512
        filename = f"{savefolder}/Sniffing to first 3_8_8 presentations_sniffs{sess_ids}.png"
        plt.savefig(filename, dpi=400, bbox_inches="tight")
    plt.show()



#%% Sniffing to opto blanks vs non-opto blanks
############### Calculating data

#t_range = np.array((4.6, 6.6)) * fps
t_range = np.array((4.4, 7.4)) * fps
bl_range = np.array((0.5, 3.5)) * fps

t_time = (t_range[1] - t_range[0]) / fps
bl_time = (bl_range[1] - bl_range[0]) / fps

graph = ["Blank", "Blank\n+ stim"]

opto_blank_vals = []
non_opto_blank_vals = []

for m in range(nmice):

    trial_idx = sniffs[m]["trial_idx"] - 1

    # blank trials excluding chem 117
    opto_blanks = (
        (sniffs[m]['trial_blank'] == 1) &
        (sniffs[m]['trial_opto'] == 1)
    )

    non_opto_blanks = (
        (sniffs[m]['trial_blank'] == 1) &
        (sniffs[m]['trial_opto'] == 0)
    )

    opto_blank_idxs = trial_idx[opto_blanks]
    non_opto_blank_idxs = trial_idx[non_opto_blanks]

    # sniffing during test window
    opto_blank_data = np.array([
        len(np.where((i > t_range[0]) & (i < t_range[1]))[0])
        for i in sniffs[m]['ml_inh_onsets'][opto_blank_idxs]
    ]) / t_time

    non_opto_blank_data = np.array([
        len(np.where((i > t_range[0]) & (i < t_range[1]))[0])
        for i in sniffs[m]['ml_inh_onsets'][non_opto_blank_idxs]
    ]) / t_time

    # sniffing during baseline window
    opto_blank_bl = np.array([
        len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0])
        for i in sniffs[m]['ml_inh_onsets'][opto_blank_idxs]
    ]) / bl_time

    non_opto_blank_bl = np.array([
        len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0])
        for i in sniffs[m]['ml_inh_onsets'][non_opto_blank_idxs]
    ]) / bl_time

    # baseline-subtracted mean per mouse
    if len(opto_blank_data) > 0:
        opto_blank_vals.append((opto_blank_data - opto_blank_bl).mean())
    else:
        opto_blank_vals.append(np.nan)

    if len(non_opto_blank_data) > 0:
        non_opto_blank_vals.append((non_opto_blank_data - non_opto_blank_bl).mean())
    else:
        non_opto_blank_vals.append(np.nan)

data = np.array([non_opto_blank_vals, opto_blank_vals])  # matches graph order

############### Plotting
def p_to_stars(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 5e-2 else 'ns'

def add_sig_bar(ax, x1, x2, y, h, stars):
    xi1, xi2 = xpos[x1], xpos[x2]
    ax.plot([xi1, xi1, xi2, xi2],
            [y, y+h, y+h, y],
            lw=1, c='k', clip_on=False, zorder=5)
    ax.text((xi1+xi2)/2, y + h*1.05, stars,
            ha='center', va='bottom',
            clip_on=False, zorder=6)

xpos = {label: i for i, label in enumerate(graph)}
nconds = len(graph)

fig_width = max(4, nconds * 1.6)
plt.figure(figsize=(fig_width, 5), dpi=300)

seq_arr = np.zeros((nmice, nconds))

for m in range(nmice):
    avg = data[:, m]
    seq_arr[m, :] = avg

    # per-mouse trace
    plt.plot(graph, seq_arr[m], color='gray', linestyle='-', alpha=0.9)

    # markers
    plt.plot(graph[0], seq_arr[m, 0], color='gray',   marker='o', linestyle='-', alpha=0.8)
    plt.plot(graph[1], seq_arr[m, 1], color='purple', marker='o', linestyle='-', alpha=0.8)

plt.ylabel("Δ avg inhalations/sec", fontsize=17)
plt.title("Mean sniffing change - blank trials", pad=40, fontsize=17, weight='bold')

ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)

# colors
bar_colors = ['gray', 'purple']

# group means and SEM across mice
means = np.nanmean(seq_arr, axis=0)
ns = np.sum(np.isfinite(seq_arr), axis=0).astype(float)
sems = np.nanstd(seq_arr, axis=0, ddof=1) / np.sqrt(np.maximum(ns, 1))

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

# ---- significance bar ----
ymin = np.nanmin(seq_arr)
ymax = np.nanmax(seq_arr)
yrng = (ymax - ymin) if np.isfinite(ymax - ymin) and (ymax - ymin) > 0 else 1.0

base_y = ymax + 0.02 * yrng
line_h = 0.01 * yrng
margin = 0.1 * yrng

a = seq_arr[:, 0]
b = seq_arr[:, 1]
mask = np.isfinite(a) & np.isfinite(b)

if mask.sum() >= 2:
    stat, p = stats.wilcoxon(a[mask], b[mask])
    stars = p_to_stars(p)
else:
    stars = 'n/a'

add_sig_bar(ax, graph[0], graph[1], base_y, line_h, stars)

ax.set_ylim(bottom=ymin - 0.05 * yrng, top=base_y + margin)

plt.xticks(np.arange(len(graph)), graph, fontsize=17)

if savefig:
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    plt.savefig(f"{savefolder}/Mean sniffing change_blank_opto_vs_nonopto_{sess_ids}.png",
                dpi=300, bbox_inches="tight")

plt.show()
#%% Plot individual trials overlayed on average trace per mouse, for first 3 novel presentations and all familiar/opto presentations - only use if warranted by habituation curves (flat)
import numpy as np
import matplotlib.pyplot as plt
import datetime

bin_edges = np.linspace(-3.5, 7.5, num=12)

for m in range(nmice):

    fig, axes = plt.subplots(1, 3, figsize=(12,4), sharey=True)

    for g in range(len(graph)):

        ax = axes[g]

        # ---- pull trials ----
        trials = sniff_arr[g, m, 6:]                        # (n_trials, n_bins)
        baseline = trials[:, baseline_bins].mean(axis=1, keepdims=True)
        trials_bl = trials - baseline

        mean = trials_bl.mean(axis=0)
        sem  = trials_bl.std(axis=0, ddof=1) / np.sqrt(trials_bl.shape[0])

        # ---- axes styling ----
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
            bin_edges,
            mean,
            yerr=sem,
            fmt='o',
            color=colors[g],
            ecolor=ecolor[g],
            elinewidth=1,
            capsize=3
        )

        # ---- stim shading ----
        if isstim:
            ax.axvspan(0.6, 2.6, color="red", alpha=0.3, linewidth=0)

        # ---- inhalation shading ----
        if isinh:
            x_values = [0] + [3 + i * 0.1 for i in range(10)]
            alpha_values = [0.33] + [0.3 - i * 0.03 for i in range(10)]

            for i in range(len(x_values)-1):
                ax.axvspan(
                    x_values[i],
                    x_values[i+1],
                    color="blue",
                    alpha=alpha_values[i],
                    linewidth=0
                )

    fig.suptitle(f"Sniffing individual trials ({sess_ids[m]})")

    plt.tight_layout()

    if savefig:
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        plt.savefig(
            f"{savefolder}/Sniffing_trials_subplot_individual_trial_traces_{sess_ids[m]}_{now}.png",
            dpi=400,
            bbox_inches="tight"
        )

    plt.show()



# %% Plotting absolute sniffing responses per odor

t_range = np.array((4.6,6.6))*fps
#t_range = np.array((4.6,10))*fps
t_time = (t_range[1]-t_range[0])/fps
bl_range = np.array((.5,3.5))*fps
bl_time = (bl_range[1]-bl_range[0])/fps

for m in range(nmice):
    
    trial_idx = sniffs[m]["trial_idx"] - 1 

    is_nov = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    is_fam = (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    is_opto = (sniffs[m]['trial_opto'] != 0) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_blank'] == 0)
    conditions = [is_nov, is_fam, is_opto]
    cond_names = ["novel", "familiar", "opto"]
    cond_colors = ['green' , 'purple', 'violet']
    fig, ax = plt.subplots(1, len(conditions), figsize=(15, 5))
    for c in range(len(conditions)):
        odors = np.unique(sniffs[m]['trial_chem_id'][conditions[c]])
        #print(cond_names[c], odors)
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
        #plot bar plots
        ax[c].bar(xticks, data, color=cond_colors[c])
        ax[c].errorbar(xticks, data, yerr=odor_sem_ls, fmt='o', color='black', ecolor='black', elinewidth=1, capsize=3)
        ax[c].set_title(f"{cond_names[c]}", fontweight='bold')
        ax[c].set_xticks(xticks, labels=odors)
        ax[c].set_xlabel("Odor")
        ax[c].set_ylabel("inhalations/sec")
        ax[c].set_ylim(0,10)

        fig.suptitle(f"Breathing change per odor {sniffs[m]['animal_name']}", fontsize=16, fontweight='bold')
        #plt.savefig(rf"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\per odor\Absolute_breathing_per_odor_{sniffs[m]['folder_identifier']}.png", dpi=300, bbox_inches="tight")

#%% Plotting distribution of odor-responses per condition

for m in range(nmice):
    plt.figure()
    
    trial_idx = sniffs[m]["trial_idx"] - 1 

    is_nov = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    is_fam = (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    is_opto = (sniffs[m]['trial_opto'] != 0) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_blank'] == 0)
    conditions = [is_nov, is_fam, is_opto]
    cond_names = ["novel", "familiar", "opto"]
    cond_colors = ['green' , 'purple', 'violet']   
    fig, ax = plt.subplots(1, len(conditions), figsize=(15, 2))
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

        xticks = np.arange(0, len(data))
        #plot bar plots
        ax[c].hist(data, bins=6)
        ax[c].set_title(f"{cond_names[c]}", fontweight='bold')
        ax[c].set_xlabel("Odor")
        ax[c].set_xlim(-2,6)
        ax[c].set_ylabel("Δ avg inhalations/sec")
        #ax[c].set_ylim(-2,6)
        ax[c].axvline(x=0, color="black",linewidth=1, alpha = 0.7)

        fig.suptitle(f"Breathing change per odor {sniffs[m]['animal_name']}", fontsize=16, fontweight='bold')

#%% Breathing change per odor per condition

t_range = np.array((4.6,6.6))*fps
#t_range = np.array((4.6,10))*fps
t_time = (t_range[1]-t_range[0])/fps
bl_range = np.array((.5,3.5))*fps
bl_time = (bl_range[1]-bl_range[0])/fps

for m in range(nmice):
    
    trial_idx = sniffs[m]["trial_idx"] - 1 

    is_nov = (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    is_fam = (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)
    is_opto = (sniffs[m]['trial_opto'] != 0) & (sniffs[m]['trial_occur'] <= 3) & (sniffs[m]['trial_blank'] == 0)
    conditions = [is_nov, is_fam, is_opto]
    cond_names = ["novel", "familiar", "opto"]
    cond_colors = ['green' , 'purple', 'violet']
    fig, ax = plt.subplots(1, len(conditions), figsize=(15, 5))
    for c in range(len(conditions)):
        odors = np.unique(sniffs[m]['trial_chem_id'][conditions[c]])
        #print(cond_names[c], odors)
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
        #plot bar plots
        ax[c].bar(xticks, data, color=cond_colors[c])
        ax[c].errorbar(xticks, data, yerr=odor_sem_ls, fmt='o', color='black', ecolor='black', elinewidth=1, capsize=3)
        ax[c].set_title(f"{cond_names[c]}", fontweight='bold')
        ax[c].set_xticks(xticks, labels=odors)
        ax[c].set_xlabel("Odor")
        ax[c].set_ylabel("Δ avg inhalations/sec")
        ax[c].set_ylim(-2,6)

        fig.suptitle(f"Breathing change per odor {sniffs[m]['animal_name']}", fontsize=16, fontweight='bold')
        plt.savefig(f"{savefolder}/Breathing_change_per_odor_{sess_ids}.png", dpi=300, bbox_inches="tight")


#%% Baseline over presentations
save_dir = savefolder + '/Baseline habituation'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)    

# ---- make subplots ----
n_cols = 1
n_rows = nmice

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows), sharex=True, sharey=True)
axes = np.array(axes).reshape(-1)  # flatten in case of grid

x_ticks = np.arange(n_presentations)
mean_ls = []

for m in range(nmice):
    ax = axes[m]

    # collect per-odor responses
    baseline_data = []
    for o in range(n_odors):
        # mean during baseline
        data = m_data[m][:,:,:,baseline_BoI].mean(axis = (0,2,3))

        baseline_data.append(data)

    baseline_data = np.array(baseline_data)

    # mean and SEM across odors
    mean = baseline_data.mean(axis=0)
    sem = baseline_data.std(axis=0, ddof=1) / np.sqrt(baseline_data.shape[0])
    mean_ls.append(mean)

    # shaded SEM
    ax.plot(x_ticks, mean, color="black", markersize=5, label=sess_ids[m])
    ax.fill_between(
        x_ticks,
        mean - sem,
        mean + sem,
        color="gray",
        alpha=0.3,
        linewidth=0
    )
    #print(baseline_data[0:5])
    
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
