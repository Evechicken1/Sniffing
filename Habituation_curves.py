"""
Habituation curves for all odors (green = novel; purple = familiar)

@author: Xander
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\xaand\Documents\PhD\Analysis\Informatics\PhD-patch-seq\Sniffing')
import sniff_tools as st
from scipy import stats
import pandas as pd
import os

#### Import sniffs dictionary ####

#!!! Session ids should  be of one type (stim, inhibition or non-opto, familiarisation)
#!!! Imports list of session id strings in this format ['250118_KK152', '250118_KK152']
#!!! The days need to be consecutive and for the same mouse
#!!! IS the path correct?

#odor_LUT_dict = 

#sess_ids = ['250118_KK152', '250118_KK153', '250118_KK154','250119_KK152','250119_KK153','250119_KK154','250120_KK152','250120_KK153','250120_KK154',"250121_KK152", "250121_KK153", "250121_KK154"]
#sess_ids = ['250814_XT001']
#sess_ids = ['250814_XT001','250814_XT002','250814_XT003','250814_XT004','250814_XT005',]
#sess_ids = ['250818_XT001','250818_XT002','250818_XT004','250818_XT005']
#sess_ids = ['250822_XT001','250822_XT002','250822_XT003','250822_XT004','250822_XT005']
#sess_ids = ['250930_XT007','250930_XT008']
#sess_ids = ['260302_XT009','260302_XT010','260302_XT011', '260302_XT012']
#sess_ids = ['250818_XT001','250818_XT002','250818_XT004','250818_XT005','251004_XT006','251004_XT007','251004_XT008','260302_XT009','260302_XT010','260302_XT011', '260302_XT012']
sess_ids = ['260306_XT009','260306_XT010','260306_XT011', '260306_XT012']

paths_list = [f"Z:\postprocessed_data\{sess_id}" for sess_id in sess_ids]
sniffs = st.import_sniff_mat_select(paths_list)
if len(sniffs)>0:
    print('Succesfully imported', len(sniffs), 'mice',[sniffs[m]["animal_name"] for m in range(len(sniffs))])

#%% Structure data into dataframes
##Variables to adjust

trialcorr =0 # 0 applies baseline correction per odor exposure trial; 1 applies baseline correction per mouse(day)
savefig = 0 #save figures to folder
savefolder = "C:/Users/xaand/Documents/PhD/Experiments/Ephys OFC/Analysis"
save_dir = savefolder + '/Habituation'
if not os.path.exists(save_dir):
    os.makedirs(save_dir) 

#%% Structuring data for usage
#initialize some variables to not error
isstim = any(sniffs[0]["trial_opto"]==1)
isinh = any(sniffs[0]["trial_opto"]==-1)
n_mice = len(sniffs)
n_presentations = np.flatnonzero((c := np.bincount(sniffs[0]['trial_occur'])) == c.max())[-1] # takes into account that blank odor can be presented 4 times extra
n_trials = sniffs[0]["trial_idx"].size
n_frames = 720
n_bins = 12
bin_size = n_frames/n_bins
bins = np.linspace(0,720,13)
baseline_BoI = range(0,4)


#Define all variables based on if sniffs[0] is an opto stim or inh or non-opto trial
if isstim: #stim variables
    nconds = 3
#    conditions = [isnov,isfam, isopto]    
    colors = ['green' , 'purple', 'violet'] 
    ecolor = ["darkgreen","purple", 'violet']
    graph = ["novel","familiar", "fam+stim"]
    markercolors = ['purple','green','purple'] 
    BoI = range(5, 8)
    seq = [1,0,2]
    comp_grps= ["nov vs fam", "nov vs opto", "fam vs opto"]


elif isinh: #inh variables
    nconds = 3
#    conditions = [isnov,isfam,isopto]    
    colors = ['green' , 'purple', 'green']
    ecolor = ["darkgreen","purple", 'green']
    graph = ["novel","familiar", 'nov+inh']
    markercolors = [ 'purple', 'green' , 'purple'] 
    BoI = range(5, 8)
    seq = [1,0,2]
    comp_grps= ["nov vs fam", "nov vs opto", "fam vs opto"]


else: #non-opto variables
    nconds = 2
#    conditions = [isnov,isfam]    
    graph = ["familiar","novel"]
#    conditions = [isnov,isfam]
    colors = ['green' , 'purple']
    ecolor = ['darkgreen','purple']
    graph = ['familiar', 'novel']
    BoI = range(5, 8)
    seq = [1,0]
    comp_grps= ["nov vs fam"]


# Defining of the selector (ntrials x conditions) !! important to accurately predefine conditions before selecting!!
selector_crits = {}

sniff_crits = {'is_nov': None, 'is_fam': None, 'is_blank': None, 'is_opto': None, 'is_123': None}
sniffs_selector = np.zeros(shape= (n_mice,n_trials, len(sniff_crits)))

## nov-fam-opto-blank selector
for m in range(n_mice):
    trial_idx = sniffs[m]['trial_idx'] - 1
    sniff_crits['is_nov'] = sniffs[m]['trial_novelty'] == 1
    sniff_crits['is_fam'] = sniffs[m]['trial_familiarity'] == 1
    sniff_crits['is_blank'] = sniffs[m]['trial_chem_id'] == 117
    if isstim or isinh:
        sniff_crits['is_opto'] = sniffs[m]['trial_opto'] != 0
    sniff_crits['is_123'] = sniffs[m]['trial_occur']<=3

    sniff_crits = pd.DataFrame(sniff_crits)
    selector_tmp = pd.DataFrame(data = np.zeros(shape= (n_trials,sniff_crits.shape[1])), columns= sniff_crits.keys())
    for c in sniff_crits:  
        idxs = trial_idx[np.where(sniff_crits[c] == True)]
        selector_tmp.loc[idxs,c] = 1
    selector_tmp = selector_tmp.to_numpy()
    sniffs_selector[m] = selector_tmp
    
### presentation TRIAL selector array
presentation_list = [[i for i in range(1,n_presentations+1)] for m in range(n_mice)]
pres_selector = np.zeros(shape= (n_mice,n_trials, n_presentations ))

for m in range(n_mice):
    trial_idx = sniffs[m]['trial_idx'] - 1
    pres_crits = {o: sniffs[m]['trial_occur'] == o for o in presentation_list[m]}
    pres_crits = pd.DataFrame(pres_crits)
    pres_tmp = pd.DataFrame(data = np.zeros(shape= (n_trials,pres_crits.shape[1])), columns= presentation_list)
    for c in presentation_list[m]:  
        idxs = trial_idx[np.where(pres_crits[c] == True)]
        pres_tmp.loc[idxs,c] = 1
    
    pres_selector[m] = pres_tmp
    
### odor TRIAL selector array
odor_list = [list(set(sniffs[m]['trial_chem_id'])) for m in range(n_mice)]
for i in range(len(odor_list)): odor_list[i].remove(117)
odor_selector = np.zeros(shape= (n_mice,n_trials, len(odor_list[0])))
n_odors = len(odor_list[0])

for m in range(n_mice):
    trial_idx = sniffs[m]['trial_idx'] - 1
    odor_crits = {o: sniffs[m]['trial_chem_id'] == o for o in odor_list[m]}
    odor_crits = pd.DataFrame(odor_crits)
    odor_tmp = pd.DataFrame(data = np.zeros(shape= (n_trials,odor_crits.shape[1])), columns= odor_crits.keys())
    for c in odor_crits:  
        idxs = trial_idx[np.where(odor_crits[c] == True)]
        odor_tmp.loc[idxs,c] = 1
    selector_tmp = odor_tmp.to_numpy()
    odor_selector[m] = odor_tmp      


#%% Plotting habituation per odor per mouse
m_data = []

for m in range(n_mice):
    plt.figure(figsize = (n_mice*1.5,5))
    trial_idx = sniffs[m]['trial_idx'] - 1
    odors = np.unique(sniffs[m]['trial_chem_id'])
    nov_data = []
    fam_data = []
    opto_data = []
    tmp = []

    nov_trials = [(sniffs[m]['trial_chem_id'] == odors[o]) & (sniffs[m]['trial_novelty'] == 1) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0)for o in range(len(odors))]
    nov_trials = [i for i in nov_trials if sum(i)>0]
    fam_trials = [(sniffs[m]['trial_chem_id'] == odors[o]) & (sniffs[m]['trial_familiarity'] == 1) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_opto'] == 0) & (sniffs[m]['trial_blank'] == 0) for o in range(len(odors))]
    fam_trials = [i for i in fam_trials if sum(i)>0]
    opto_trials = [(sniffs[m]['trial_chem_id'] == odors[o]) & (sniffs[m]['trial_chem_id'] != 117) & (sniffs[m]['trial_opto'] != 0) & (sniffs[m]['trial_blank'] == 0) for o in range(len(odors))]
    opto_trials = [i for i in opto_trials if sum(i)>0]
    
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
    


    fig, ax = plt.subplots(1, len(conditions), figsize=(20,5), sharex=True, sharey=True)
    
    for i in range(len(conditions)):
        x_ticks = np.arange(n_presentations)
        data = conditions[i]
        baseline = data[:,:,baseline_BoI].mean(axis=(2))
        mean = data[:,:,BoI].mean(axis=(2))
        sem = np.std(mean, axis=1) / np.sqrt(mean.shape[0])

        corr_mean = mean - baseline
        corr_mean = corr_mean.mean(axis=1)        

        
        ax[i].plot(x_ticks, corr_mean, color = colors[i])
        ax[i].fill_between(x_ticks, corr_mean +sem , corr_mean - sem, color = ecolor[i], alpha = 0.5, linewidth = 0)
        ax[i].axhline(0, color='black', linestyle='--')
        ax[i].set_title(graph[i])
        
    fig.suptitle(f'Habituation {sess_ids[m]}', fontweight = 'bold')
    fig.supxlabel("Presentation #")
    fig.supylabel("")
    plt.savefig(r"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\habituation_curves" + "\ " + sess_ids[m] + '.png', dpi = 300)
    plt.show()    
    
    m_data.append(conditions)
m_data = np.array(m_data) # convert to array for easier indexing

#%% Plotting average habituation across mice per condition
plt.figure(figsize=(5,5))
m_data.shape # (n_mice, n_conditions, n_presentations, n_odors, n_bins)
for i in range(len(conditions)):
    x_ticks = np.arange(n_presentations)
    data = m_data[(0,1),i].mean(axis=0) # select condition i across mice
    baseline = data[:,:,baseline_BoI].mean(axis=(2))
    mean = data[:,:,BoI].mean(axis=(2))
    sem = np.std(mean, axis=1) / np.sqrt(mean.shape[0])

    corr_mean = mean - baseline
    corr_mean = corr_mean.mean(axis=1)        

    plt.plot(x_ticks, corr_mean, color = colors[i], label=graph[i])
    plt.fill_between(x_ticks, corr_mean +sem , corr_mean - sem, color = ecolor[i], alpha = 0.5, linewidth = 0)
plt.axhline(0, color='black', linestyle='--')  
plt.title('Habituation')
plt.xlabel("Presentation #")
plt.ylabel("Inhalations/second (baseline-subtracted)")
plt.legend()
plt.savefig(r"C:\Users\xaand\Documents\PhD\Experiments\Opto OFC-LDTg\Analysis\habituation_curves\all_mice.png", dpi = 300)

plt.show()


#%% Plotting average habituation per mouse
# ---- make subplots ----
n_cols = 1
n_rows = n_mice

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows), sharex=True, sharey=True)
axes = np.array(axes).reshape(-1)  # flatten in case of grid

x_ticks = np.arange(n_presentations)


for m in range(n_mice):
    ax = axes[m]

    # collect responses per odor per presentation
    per_odor_data = []
    n_conds = len(m_data[m])
    for c in range(n_conds):
        n_presentations = m_data[m][c].shape[0]
    
        # mean during BoI
        data_tmp = m_data[m][c][:, :,BoI].mean(axis=2)

        # mean during baseline
        data_tmp_bl = m_data[m][c][:, :, baseline_BoI].mean(axis=2)
    
        # subtract baseline
        data_sub = data_tmp - data_tmp_bl

        per_odor_data.append(data_sub)
    per_odor_data = np.concatenate(per_odor_data, axis=1)  # concatenate odors along columns

    # mean and SEM across odors
    mean = per_odor_data.mean(axis=1)
    sem = per_odor_data.std(axis=1) / np.sqrt(per_odor_data.shape[1])  # SEM across odors

    # shaded SEM
    ax.plot(x_ticks, mean, color="black", markersize=5, label=f"Mouse {m+1}")
    ax.fill_between(
        x_ticks,
        mean - sem,
        mean + sem,
        color="gray",
        alpha=0.3,
        linewidth=0
    )

    ax.set_title(sess_ids[m])
    ax.set_xlabel("Presentation")
    ax.set_ylabel("Inhalations/second (baseline-subtracted)")

plt.suptitle('Habituation', fontweight="bold", y = .91)
plt.savefig(save_dir + "\ " + str(fr"odor responses {sniffs[m]['folder_identifier']}") + '.png', dpi = 300)
plt.show()
   
#%% baseline evolution over presentations     
    
save_dir = savefolder + '/Baseline habituation'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)    

# ---- make subplots ----
n_cols = 1
n_rows = n_mice

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows), sharex=True, sharey=True)
axes = np.array(axes).reshape(-1)  # flatten in case of grid

x_ticks = np.arange(n_presentations)
mean_ls = []

for m in range(n_mice):
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

   


    
#%% baseline evolution over presentations (blanks) 

save_dir = savefolder + '/habituation to blanks'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)    
    
# ---- make subplots ----
n_cols = 1
n_rows = n_mice

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows), sharex=True, sharey=True)
axes = np.array(axes).reshape(-1)  # flatten in case of grid

x_ticks = np.arange(n_presentations)

for m in range(n_mice):
    ax = axes[m]

    # collect per-odor responses
    baseline_data = []
    for o in range(n_odors):
        # mean during baseline
        data = m_data[m][:,:,:,baseline_BoI].mean(axis = (0,2,3))

        baseline_data.append(data)

    baseline_data = np.array(baseline_data)

    # mean and SEM across odors
    mean = baseline_data
    sem = baseline_data.std(axis=0, ddof=1) / np.sqrt(baseline_data.shape[0])

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
    ax.set_ylabel("Breathing Change")
plt.suptitle('Baseline change over exposures', fontweight="bold", y = .91)


plt.show()


    
    


















