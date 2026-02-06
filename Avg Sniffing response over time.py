import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\xaand\Documents\PhD\Analysis\Informatics\PhD-patch-seq\Sniffing')
import sniff_tools as st
from scipy import stats
import datetime
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

#### Import sniffs dictionary ####

#!!! Session ids should  be of one type (stim, inhibition or non-opto)
#!!! Imports list of session id strings in this format ['250118_KK152', '250118_KK152']
#!!! IS the path correct?

#sess_ids = ['231020_KK066', '231215_KK082', '231215_KK085', '240211_KK092', '240211_KK093'] # DA inh (Xander) 
#sess_ids = ['231024_KK066', '231219_KK082', '231219_KK085', '240215_KK092', '240215_KK093'] # DA stim (Xander) 
#sess_ids = ['240201_KK087', '240201_KK088', '240201_KK089', '240301_KK099', '240301_KK100'] # DA CTRL inh (Xander)
#sess_ids = ['240205_KK087', '240205_KK088', '240205_KK089', '240305_KK099', '240305_KK100'] # DA CTRL stim (Xander)
#sess_ids = ['250118_KK152', '250118_KK153', '250118_KK154','250119_KK152','250119_KK153','250119_KK154','250120_KK152','250120_KK153','250120_KK154',"250121_KK152", "250121_KK153", "250121_KK154"]
#sess_ids = ['250818_XT001','250818_XT002','250818_XT004','250818_XT005']
#sess_ids = ['250822_XT001','250822_XT002','250822_XT003','250822_XT004','250822_XT005']
#sess_ids = ['251004_XT006','251004_XT007','251004_XT008']
sess_ids = ['250818_XT001','250818_XT002','250818_XT004','250818_XT005','251004_XT006','251004_XT007','251004_XT008']
#sess_ids = ['250818_XT005']

paths_list = [f"Z:\postprocessed_data\{sess_id}" for sess_id in sess_ids]
sniffs = st.import_sniff_mat_select(paths_list)

if len(sniffs)>0:
    print('Succesfully imported', len(sniffs), 'mice',[sniffs[m]['folder_identifier'] for m in range(len(sniffs))])

#%% Structure data into dataframes
##Variables to adjust

trialcorr = 0 #applies baseline correction per odor exposure trial; 1 applies baseline correction per mouse(day)
savefig = 0 #save figures to folder
savefolder = r"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\Sniffing to nov fam"
n_presentations = 6
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




#Define all variables based on if sniffs[0] is an opto stim or inh or non-opto trial
if isstim: #stim variables
    nconds = 3
    conditions = [isnov,isfam, isopto]    
    colors = ['green' , 'purple', 'purple'] 
    ecolor = ["darkgreen","purple", 'purple']
    graph = ["novel","familiar", "fam+stim"]
    markercolors = ['purple','green','purple'] 
    BoI = range(5, 8)
    seq = [1,0,2]
    comp_grps= ["nov vs fam", "nov vs opto", "fam vs opto"]
    line_type = ["solid", "solid","dotted"]


elif isinh: #inh variables
    nconds = 3
    conditions = [isnov,isfam,isopto]    
    colors = ['green' , 'purple', 'green']
    ecolor = ["darkgreen","purple", 'green']
    graph = ["novel","familiar", 'nov+inh']
    markercolors = [ 'purple', 'green' , 'purple'] 
    BoI = range(5, 8)
    seq = [1,0,2]
    comp_grps= ["nov vs fam", "nov vs opto", "fam vs opto"]
    line_type["solid", "solid","dotted"]

else: #non-opto variables
    nconds = 2
    conditions = [isnov,isfam]    
    graph = ["novel","familiar"]
    conditions = [isnov,isfam]
    colors = ['green' , 'purple']
    ecolor = ['darkgreen','purple']
    graph = ['familiar', 'novel']
    BoI = range(5, 8)
    seq = [1,0]
    comp_grps= ["nov vs fam"]
    line_type["solid", "solid"]

# Defining of the selector (ntrials x conditions) !! important to accurately predefine conditions before selecting!!
ntrials = sniffs[0]["trial_idx"].size
selector = [np.zeros((ntrials, len(conditions)), int) for i in range(nmice)]   

for m in range(nmice):
    #print("Mouse", m)
    isnov = (sniffs[m]["trial_novelty"]== True) & (sniffs[m]["trial_opto"]== False)  & (sniffs[m]["trial_occur"] <= 3) #& (sniffs[m]["trial_blank"] == 0)
    isfam = (sniffs[m]["trial_familiarity"]== True) & (sniffs[m]["trial_opto"]== False) & (sniffs[m]["trial_occur"] <= 3) #& (sniffs[m]["trial_blank"] == 0)
    isopto = (sniffs[m]["trial_opto"] != False) & (sniffs[m]["trial_reward"] == False) & (sniffs[m]["trial_blank"] == 0) & (sniffs[m]["trial_chem_id"] != 117) & (sniffs[m]["trial_occur"] <= 3) 
    
#    isnov = (sniffs[m]["trial_novelty"]== True) & (sniffs[m]["trial_opto"]== False)  & (sniffs[m]["trial_occur"] >= 4)
#    isfam = (sniffs[m]["trial_familiarity"]== True) & (sniffs[m]["trial_opto"]== False) & (sniffs[m]["trial_occur"] >= 4)
#    isopto = (sniffs[m]["trial_opto"] != False) & (sniffs[m]["trial_reward"] == False) & (sniffs[m]["trial_blank"] == 0) & (sniffs[m]["trial_chem_id"] != 117) & (sniffs[m]["trial_occur"] >= 4)
    

    conditions = [isnov, isfam, isopto] if (isstim or isinh) else [isnov, isfam]

    # join condition1 with another condition - minimal trial occurence
    trial_idx = sniffs[m]["trial_idx"]-1
    
    for cond in range(selector[m].shape[1]):
        idx_cnds = trial_idx[conditions[cond]]
#        print(cond,idx_cnds)
        selector[m][idx_cnds,cond] = True
#        print(selector[m][:,cond])

nsel = selector[0][:,0].sum()
ngraph = len(conditions)

#%% Create a list of arrays(36x12) for each condition, which we can use for plotting the different conditions
sniff_lst = [[[]for m in range(nmice)] for g in range(len(graph))] # make a list of arrays in size of selector sum of selected category
hist_edges = np.linspace(0, 723,13)

mean_arr = np.zeros((len(graph), (nmice),  12))
corrmean_arr = np.zeros((len(graph), nmice,  12))

std_arr = np.zeros((len(graph), (nmice),  12))
baseline_lst = [[[] for n in range(nmice)] for g in range(len(graph))]
baseline_arr = np.zeros((len(graph), (nmice),  18))

for g in range(len(graph)):
#    print("graph", graph[g])
    for m in range(len(selector)):
#        print("Mouse", m)
        
        tmp = sniffs[m]["ml_inh_onsets"][selector[m][:,g] == 1]
        for i in range(len(tmp)):
            hist = np.histogram(tmp[i], bins=hist_edges)[0]
            #print(hist, i) 
            sniff_lst[g][m].append(hist)
            baseline_lst[g][m].append(np.mean(sniff_lst[g][m][i][baseline_bins]))

sniff_arr = np.array(sniff_lst)

#%% Plotting average sniffs for all mice, during the trial

bin_edges = np.linspace(-3.5, 7.5, num = 12)

sem = [[[] for m in range(nmice)] for l in range(ngraph)]
mean = [[[]for m in range(nmice)] for l in range(ngraph)]
baseline = [[[]for m in range(nmice)] for l in range(ngraph)]

plt.figure(figsize=(6,3.5))
for g in range(ngraph):
    mean = sniff_arr[g].mean(axis = (0))
    baseline = sniff_arr[g,:,:,baseline_bins].mean(axis = (0,1))
    mean = mean - baseline[:,None]
    mean = mean.mean(axis = 0)
    sem = np.std(sniff_arr[g].mean(axis = 1), axis = 0)
        
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
    filename = f"{savefolder}/Sniffing to first 3 presentations_sniffs[0]['session_date']_{now}.svg"
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
    ax = plt.gca() #to remove the top and right spines
    ax.spines[['right', 'top']].set_visible(False) #to remove the top and right spines

    for g in range(len(graph)):       
        mean = sniff_arr[g,m]
        baseline = sniff_arr[g,m,:,baseline_bins].mean(axis = 0)
        mean = mean - baseline[:,None]
        mean = mean.mean(axis = 0)
        sem = np.std(sniff_arr[g].mean(axis = 1), axis = 0)

        plt.plot(bin_edges, mean, color= colors[g], label = graph[g], ls = line_type[g])
        plt.errorbar(bin_edges, mean, yerr=sem, fmt='o', color=colors[g], ls = line_type[g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)

    if isstim: #shading for stim
       plt.axvspan(xmin=0.6, xmax=2.6, color="red", alpha=0.3, linewidth = 0)
    if isinh: #shading for inh
        x_values = [0] + [3 + i * 0.1 for i in range(10)]
        alpha_values = [0.33] + [0.3 - i * 0.03 for i in range(10)]
        [plt.axvspan(x_values[i], x_values[i + 1], color="blue", alpha=alpha_values[i], linewidth=0) for i in range(len(x_values) - 1)]    
    
    if savefig:
        now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        plt.savefig(f"{savefolder}/Sniffing to first 3 presentations ({sess_ids[m]})_{now}.svg", dpi=400, bbox_inches="tight")
    plt.show()

#%% Comparison of sniffing during specific timebin
############### Calculating data

t_range = np.array((4.6,6.6))*fps
t_time = (t_range[1]-t_range[0])/fps
bl_range = np.array((.5,3.5))*fps
bl_time = (bl_range[1]-bl_range[0])/fps
graph = ["familiar", "novel","fam+stim"]

nov = []
fam = []
opto = []


for m in range(nmice):
    
    trial_idx = sniffs[m]["trial_idx"] - 1 
    nov_odors = np.logical_and(sniffs[m]['trial_novelty'] == 1,  sniffs[m]['trial_occur'] <= 8)
    fam_odors = np.logical_and(sniffs[m]['trial_familiarity'] == 1,  sniffs[m]['trial_occur'] <= 8)
    opto_odors = np.logical_and(sniffs[m]['trial_opto'] == 1,  sniffs[m]['trial_occur'] <= 8)
    blank_odors = np.logical_and(sniffs[m]['trial_blank'] == 1,  sniffs[m]['trial_occur'] <= 8)
    
    nov_idxs = trial_idx[nov_odors]
    fam_idxs = trial_idx[fam_odors]
    opto_idxs = trial_idx[opto_odors]
    blank_idxs = trial_idx[blank_odors]
    
    nov_data = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][nov_idxs]])/t_time
    fam_data = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][fam_idxs]])/t_time
    opto_data = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][opto_idxs]])/t_time

    nov_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][nov_idxs]])/bl_time
    fam_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][fam_idxs]])/bl_time
    opto_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][opto_idxs]])/bl_time

    nov_data = (nov_data - nov_bl).mean()
    fam_data = (fam_data - fam_bl).mean()
    opto_data = (opto_data - opto_bl).mean()

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
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    plt.savefig(f"{savefolder}/Mean sniffing change_{now}.svg", dpi=300, bbox_inches="tight")

plt.show()

#%% Stimulation effect per odor

t_range = np.array((4.6,6.6))*fps
t_time = (t_range[1]-t_range[0])/fps
bl_range = np.array((.5,3.5))*fps
bl_time = (bl_range[1]-bl_range[0])/fps
graph = ["familiar", "novel","fam+stim"]

for m in range(nmice):
    plt.figure()
    trial_idx = sniffs[m]["trial_idx"] - 1 
    
    odor_ls =[]
    odors = np.unique(sniffs[m]['trial_chem_id'][sniffs[m]['trial_opto'] == 1])
    for o in odors:
        opto_odors = (sniffs[m]['trial_opto'] == 1) & (sniffs[m]['trial_occur'] >= 6) & (sniffs[m]['trial_chem_id'] == o) & (sniffs[m]['trial_blank'] == 0)
        odor_idxs = trial_idx[opto_odors]
        odor_data = np.array([len(np.where((i >t_range[0]) & (i < t_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][odor_idxs]])/t_time
        odor_bl = np.array([len(np.where((i > bl_range[0]) & (i < bl_range[1]))[0]) for i in sniffs[m]['ml_inh_onsets'][odor_idxs]])/bl_time
        odor_data = np.nanmean(odor_data - odor_bl)
        odor_ls.append(odor_data)
        
        data = odor_ls
    
    xticks = np.arange(0, len(data))
    plt.bar(xticks, data)
    plt.title(sniffs[m]['animal_name'])
    plt.ylim(-2,4)
    plt.show()


























