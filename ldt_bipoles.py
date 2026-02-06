# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 10:07:17 2021

@author: kkondrakiewicz
"""
import numpy as np
import matplotlib.pyplot as plt
import copy
import sys
sys.path.append(r'C:\Users\xaand\OneDrive\PhD\Analysis\Sniffing\Code')
import sniff_tools as st
import seaborn as sns
from scipy import stats
import matplotlib as mpl

#%% Specify paths and some global analysis parameteres
#sniff_path = r'C:\Users\Kacper\Desktop\LDTg_Chat_BIPOLES\data_all'
#sniff_path = r'C:\Users\Kacper\Desktop\LDTg_VGlut_BIPOLES\data_all'
#sniff_path = r'C:\Users\Kacper\Desktop\pm_AON_DA\sniff_data'
#sniff_path = r'C:\Users\Kacper\Desktop\preBotz\data'
sniff_path = r'C:\Users\xaand\OneDrive\unief\2023-2024\Masters_Thesis\Data_Analysis\KK62_63\data'

#fig_path = r'C:\Users\Kacper\Desktop\LDTg_Chat_BIPOLES\inh_all'
#fig_path = r'C:\Users\Kacper\Desktop\LDTg_Chat_BIPOLES\stim2_all'
#fig_path = r'C:\Users\Kacper\Desktop\LDTg_VGlut_BIPOLES\stim_all_all'
#fig_path = r'C:\Users\Kacper\Desktop\pm_AON_DA\figs_behavior\stim2'
#fig_path = r'C:\Users\Kacper\Desktop\preBotz\inh'
fig_path = r'C:\Users\xaand\OneDrive\unief\2023-2024\Masters_Thesis\Data_Analysis\KK62_63\Figures'

#sess_ids = ['221114_KK036', '221114_KK037', '230221_KK042', '230221_KK043', '230221_KK044']  
#sess_ids = ['221118_KK036', '221118_KK037', '230225_KK042', '230225_KK043', '230225_KK044']

#sess_ids = ['230221_KK042', '230221_KK043', '230221_KK044']
#sess_ids = ['230227_KK042', '230227_KK043', '230227_KK044']  

#sess_ids = ['230321_KK047', '230321_KK048', '230621_KK054', \
#            '230818_KK059', '230818_KK060', '230818_KK061']  # '230621_KK055', '230621_KK056',

#sess_ids = ['230317_KK047', '230317_KK048', '230617_KK054', \
#    '230814_KK059', '230814_KK060', '230814_KK061'] # '230617_KK055', '230617_KK056', 
#sess_ids = ['230525_KK051', '230525_KK052', '230525_KK053']  
#sess_ids = ['230324_KK045', '230324_KK046']

sess_ids = ['231006_KK062', '231006_KK063']  


nframes = 723 # how many camera frames per trial you expect
pup_nframes = 373 # the same for pupil camera
sr = 60
pup_sr = pup_nframes/12
sigma = 0.25
binsize = 2 # for binned analysis, bin size in seconds
est_lat = 0 # estimated olfactometer latency
odor_start = 4 + est_lat
odor_end = 6 + est_lat
bsln_start = 1
ndays = 1
sniff_the_bin = [4.5, 7.5] # concentrate on this part - from 1 sec to 3 sec after odor presentation
pup_bin = [5.5, 8.5] # this can be different for pupil, which has slower dynamics
mpl.rcParams['svg.fonttype'] = 'none' # this should make text editable in .svg figures

nov_color = '#006400'
fam_color = '#580F41'

sniff_dirs = []
for ses in sess_ids:
    tmp = sniff_path + '\\' + ses
    sniff_dirs.append(tmp)

#%% Import odor trial data
sniffs = st.import_sniff_mat_select(sniff_dirs)


#%% Exctract some basic info from the imported data
nses = len(sniffs)
nmice = int(nses/ndays)
ntrials = sniffs[0]['trial_idx'].size
npres = max(sniffs[0]['trial_occur'])
sr = sniffs[0]['samp_freq']
figure_size = (4, nmice*3)


#%% Which odors are opto?
odors_stim = []
odors_inh = []
for m in range(nses):
    tmp = np.logical_and(sniffs[m]['trial_opto']==1, sniffs[m]['trial_reward']==0)
    odors_stim.append(np.unique(sniffs[m]['trial_type'][tmp]))
    tmp = np.logical_and(sniffs[m]['trial_opto']==-1, sniffs[m]['trial_reward']==0)
    odors_inh.append(np.unique(sniffs[m]['trial_type'][tmp]))
    
#odors_stim = [np.array([ 2,  8, 16, 19, 24, 26, 30]), np.array([ 2,  5,  9, 17, 22, 27, 28]), \
#              np.array([ 1, 10, 12, 16, 26, 27, 30])] # for KK051-53
#%%  
if len(odors_inh[0]) > 0:
    odors_opto = odors_inh
    var_color = nov_color
    var_label = 'NOV + INH'
    var_label2 = 'BLANK + INH'
elif len(odors_stim[0]) > 0:
    odors_opto = odors_stim
    var_color = fam_color
    var_label = 'FAM + STIM'
    var_label2 = 'BLANK + STIM'
    

#%% Restructure sniffing data into 3-dim array: trials x time point x mice and calculate breathing rate (multiple methods)
sniff_ons, sniff_list, sniff_bins, sniff_delbins, sniff_mybin = st.bin_sniff(sniffs, \
    nframes, bsln_start, odor_start, sniff_the_bin, binsize, classic_tracker=False)
sniff_gauss, sniff_delta = st.ins_sniff(sniff_ons, bsln_start, odor_start, sigma, sr)

#%% Create odor category matrix, indicating for each trial which odor type it is     
incl_descr = 'First 3 odor presentations'
tr_cat, tr_incl = st.select_trials_nov(sniffs, fam_min=1, fam_max=3, nov_min=1, nov_max=3)

for m in range(nmice):
    is_opto = np.isin(sniffs[m]['trial_type'], odors_opto[m])
    is_ctrl = (is_opto == 0)*1
    is_opto = np.reshape(is_opto, [np.shape(is_opto)[0], 1]) *1
    tr_cat[m] = np.hstack([tr_cat[m], is_opto]) # now column 3 is opto stimuli
    tr_cat[m][:,1] = np.logical_and(tr_cat[m][:,1], is_ctrl) # and columns 1  is only novel without inhibition
    tr_cat[m][:,0] = np.logical_and(tr_cat[m][:,0], is_ctrl) # the same applies to familiar - without stimulation
    # Just repeat the novelty / familiarity column for trial include matrix
    if odors_opto == odors_inh:
        tr_incl[m] = np.hstack([tr_incl[m], np.reshape(tr_incl[m][:,1], [tr_incl[m].shape[0], 1])]) 
    elif odors_opto == odors_stim:
        tr_incl[m] = np.hstack([tr_incl[m], np.reshape(tr_incl[m][:,0], [tr_incl[m].shape[0], 1])])
    
    # AND correct this matrix so for example novel (col 1) are only novels without inhibition
    tr_incl[m] = np.logical_and(tr_incl[m], tr_cat[m])*1

ncat = tr_cat[0].shape[1]

#%% Calculate mean and SEM for each occurence of a given trial type
sniff_1bin_av, sniff_1bin_n, sniff_1bin_sem = st.av_by_occur(sniffs, sniff_mybin, tr_cat)

#%% Calculate mean sniffing across time for selected presentation
sniff_av = np.zeros([nframes, nses, ncat])
sniff_n = np.sum(tr_incl[0], 0)
sniff_sem = sniff_av.copy()

for m in range(nses):
    for cat in range(ncat):
        which_rows = tr_incl[m][:,cat] == 1
        which_incl = sniffs[m]['trial_idx'][which_rows] - 1 # IN MATLAB TRIAL INDEXES START FROM 1!!
        
        #print(np.where(which_rows==True))
        # To manually verify that we selected right trials, print event indexes in each cat
        to_check = which_incl*4+1 # go from trial index to event number
        to_check.sort()
        #print(to_check)
        
        tmp_data = sniff_delta[which_incl, :, m].T
        #tmp_data = sniff_gauss[which_incl, :, m].T
        
        # Look at raw data from the selected trials
        #plt.figure(); plt.plot(tmp_data)
        #plt.ylim([-2, 8])
        
        sniff_av[:,m,cat] = np.mean(tmp_data, 1)
        sniff_sem[:,m,cat] = np.std(tmp_data, 1) / np.sqrt(sniff_n[cat])
             

#%% Plot breathing across time for some selected trials
fig, axes = plt.subplots(int(nmice), 1, sharex = 'all', sharey='all', figsize = figure_size, dpi = 250)
axes = axes.flatten()
tvec = np.linspace(-4, 8, nframes) - est_lat
    
for m in range(nmice):
    
    axes[m].plot(tvec, sniff_av[:,m,0], label = 'FAM', color = fam_color, linestyle = '-')
    axes[m].fill_between(tvec, sniff_av[:,m,0] + sniff_sem[:,m,0], sniff_av[:,m,0] - sniff_sem[:,m,0], alpha = 0.2, color = fam_color)
    
    axes[m].plot(tvec, sniff_av[:,m,1], label = 'NOV', color = nov_color, linestyle = '-')
    axes[m].fill_between(tvec, sniff_av[:,m,1] + sniff_sem[:,m,1], sniff_av[:,m,1] - sniff_sem[:,m,1], alpha = 0.2, color = nov_color)
    
    axes[m].plot(tvec, sniff_av[:,m,3], label = var_label, color = var_color, linestyle = '--')
    axes[m].fill_between(tvec, sniff_av[:,m,3] + sniff_sem[:,m,3], sniff_av[:,m,3] - sniff_sem[:,m,3], alpha = 0.2, color = var_color)
  

    axes[m].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
    axes[m].axvline(x = 2, linestyle = '--', color = 'gray', linewidth = 1)

    axes[m].set_ylabel(u"\u0394" + ' sniffing [inh/sec]')
    axes[m].set_xlim([-3,7])
    #axes[m].set_ylim([-1,3.5])
    
    ax2 = axes[m].twinx()
    ax2.set_yticks([])
    mouse_id = sniffs[m]['unique_id'][7:12]
    ax2.set_ylabel('Mouse ' + mouse_id)
    
    
axes[m].legend(loc=2)
axes[m].set_xlabel('Time from odor presentation [sec]')
fig.suptitle(incl_descr)
plt.savefig(fig_path + '\\Sniff_' + incl_descr + '.png', bbox_inches = 'tight', dpi=250)
plt.savefig(fig_path + '\\Sniff_' + incl_descr + '.svg', bbox_inches = 'tight', transparent = True)


#%% Plot habituation curve for each mouse
fig, axes = plt.subplots(int(nmice), 1, sharex = 'all', sharey='all', figsize = figure_size)
axes = axes.flatten()

pvec_f = np.arange(1,9)
pvec_n = np.arange(1,9)

for m in range(nmice):
    axes[m].errorbar(pvec_f, sniff_1bin_av[:8,m,0], sniff_1bin_sem[:8,m,0], label = 'FAM', color = fam_color, linestyle = '-', marker = 'o')
    axes[m].errorbar(pvec_n, sniff_1bin_av[:8,m,1], sniff_1bin_sem[:8,m,1], label = 'NOV', color = nov_color, linestyle = '-', marker = 'o')        
    axes[m].errorbar(pvec_n, sniff_1bin_av[:8,m,3], sniff_1bin_sem[:8,m,3], label = var_label, color = var_color, linestyle = '--', marker = 'o')  
    
    axes[m].set_ylabel(u"\u0394" + ' sniffing [inh/sec]')
    
    ax2 = axes[m].twinx()
    ax2.set_yticks([])
    mouse_id = sniffs[m]['unique_id'][7:12]
    ax2.set_ylabel('Mouse ' + mouse_id)
    
    
axes[m].legend()
axes[m].set_xlabel('Presentation number')
fig.suptitle('Habituation curve')

plt.savefig(fig_path + '\\Sniff_hab.png', bbox_inches = 'tight', dpi=250)
plt.savefig(fig_path + '\\Sniff_hab.svg', bbox_inches = 'tight', transparent = True)


#%% Plot average habituation curve
plt.figure(dpi = 250)

hab_av = np.mean(sniff_1bin_av, 1)
hab_sem = np.mean(sniff_1bin_av, 1) / np.sqrt(nmice)

plt.errorbar(pvec_f-0.1, hab_av[:8,0], hab_sem[:8,0], label = 'FAM', color = fam_color, linestyle = '-', fmt = 'o')
plt.errorbar(pvec_f, hab_av[:8,1], hab_sem[:8,1], label = 'NOV', color = nov_color, linestyle = '-', fmt = 'o')
plt.errorbar(pvec_f+0.1, hab_av[:8,3], hab_sem[:8,3], label = var_label, color = var_color, linestyle = '--', fmt = 'o', mfc='white')

#plt.ylim([-0.3,4])
plt.ylabel(u"\u0394" + ' sniffing [inh/sec]')    
plt.legend()
plt.xlabel('Presentation number')

plt.savefig(fig_path + '\\Sniff_hab_average.png')
plt.savefig(fig_path + '\\Sniff_hab_average.svg', bbox_inches = 'tight', transparent = True)


#%% Import pupil dilation data and parse it into trials

if 'data_path' in locals():
    pup_raw, pup_ts = st.import_pupil(data_path)
    
    if np.isnan(pup_ts):
        print('Assuming there are no lost frames in pupil cam')
        pup_ts = []
        pup_dt = 1/pup_sr
        for m in range(nmice):
            pup_trial_ts = np.linspace(pup_dt, pup_dt*pup_nframes, pup_nframes)
            pup_ts.append(np.zeros(pup_raw[m].shape[0] + pup_nframes))
            
            ii = 0
            for fr in range(0, pup_raw[m].shape[0], pup_nframes):
                pup_ts[m][fr:fr+pup_nframes] = (ii * 1e11) + pup_trial_ts # pretend trials are evenly spaced
                ii += 1
            
            pup_ts[m] = pup_ts[m][0:pup_raw[m].shape[0]]
            
    pup_m = st.parse_pupil(pup_raw, pup_ts, ntrials, pup_nframes, nses, smoothen=1)
    
    #%% Normalize and average pupil data
    pup_delta = pup_m.copy()
    pup_mybin = np.zeros([ntrials, nses])
    
    for m in range(nses):
        for tr in range(ntrials):
            tmp_data = pup_m[tr, :, m]
            bsl = np.nanmean(tmp_data[int(bsln_start*pup_sr) : int(odor_start*pup_sr)])
            pup_delta[tr, 0:tmp_data.size, m] = (tmp_data - bsl)
            
            pup_mybin[tr, m] = np.nanmean(pup_delta[tr, int(pup_bin[0]*pup_sr):int(pup_bin[1]*pup_sr), m])
    
      
    pup_1bin_av, pup_1bin_n, pup_1bin_sem = st.av_by_occur(sniffs, pup_mybin, tr_cat)
      
    #%% The same, but for mean pupil
    pup_av = np.zeros([pup_nframes, nses, ncat])
    pup_n = np.sum(tr_incl[0], 0)
    pup_sem = pup_av.copy()
    
    for m in range(nses):
        for cat in range(ncat):
            which_rows = np.logical_and(tr_incl[m][:,cat] == 1, tr_cat[m][:, cat] == 1) 
            which_incl = sniffs[m]['trial_idx'][which_rows] - 1 # IN MATLAB TRIAL INDEXES START FROM 1!!
            
            tmp_data = pup_delta[which_incl, :, m].T
            #tmp_data = pup_m[which_incl, :, m].T
            
            pup_av[:,m,cat] = np.nanmean(tmp_data, 1)
            pup_sem[:,m,cat] = np.nanstd(tmp_data, 1) / np.sqrt(pup_n[cat])
            


#%% Plot mean breathing trace from all animals
plt.figure()
tvec = np.linspace(-4, 8, nframes) - est_lat

tmp_av = np.mean(sniff_av[:,:,0],1)
tmp_sem = np.std(sniff_av[:,:,0],1) / np.sqrt(sniff_av[:,:,0].shape[1])
plt.plot(tvec, tmp_av, label = 'FAM', color = fam_color, linestyle = '-')
plt.fill_between(tvec, tmp_av + tmp_sem, tmp_av - tmp_sem, alpha = 0.2, color = fam_color)

tmp_av = np.mean(sniff_av[:,:,1],1)
tmp_sem = np.std(sniff_av[:,:,1],1) / np.sqrt(sniff_av[:,:,1].shape[1])
plt.plot(tvec, tmp_av, label = 'NOV', color = nov_color, linestyle = '-')
plt.fill_between(tvec, tmp_av + tmp_sem, tmp_av - tmp_sem, alpha = 0.2, color = nov_color)

tmp_av = np.mean(sniff_av[:,:,3],1)
tmp_sem = np.std(sniff_av[:,:,3],1) / np.sqrt(sniff_av[:,:,3].shape[1])
plt.plot(tvec, tmp_av, label = var_label, color = var_color, linestyle = '--')
plt.fill_between(tvec, tmp_av + tmp_sem, tmp_av - tmp_sem, alpha = 0.2, color = var_color)
  

plt.axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
plt.axvline(x = 2, linestyle = '--', color = 'gray', linewidth = 1)

plt.ylabel(u"\u0394" + ' sniffing [inh/sec]')
plt.xlim([-2,6])

    
plt.legend()
plt.xlabel('Time from odor presentation [sec]')
plt.suptitle(incl_descr)

plt.savefig(fig_path + '\\Sniff_average_' + incl_descr + '.png')
plt.savefig(fig_path + '\\Sniff_average_' + incl_descr + '.svg', bbox_inches = 'tight', transparent = True)


#%% Plot pupil across time for some selected trials

if 'data_path' in locals():

    fig, axes = plt.subplots(int(nmice), 1, sharex = 'all', sharey='all', figsize = figure_size)
    axes = axes.flatten()
    tvec_p = np.linspace(-4, 8, pup_nframes) - est_lat
    
    for m in range(nmice):
        
        axes[m].plot(tvec_p, pup_av[:,m,0], label = 'FAM', color = fam_color, linestyle = '-')
        axes[m].fill_between(tvec_p, pup_av[:,m,0] + pup_sem[:,m,0], pup_av[:,m,0] - pup_sem[:,m,0], alpha = 0.2, color = fam_color)
        
        axes[m].plot(tvec_p, pup_av[:,m,1], label = 'NOV', color = nov_color, linestyle = '-')
        axes[m].fill_between(tvec_p, pup_av[:,m,1] + pup_sem[:,m,1], pup_av[:,m,1] - pup_sem[:,m,1], alpha = 0.2, color = nov_color)
        
        axes[m].plot(tvec_p, pup_av[:,m,3], label = var_label, color = var_color, linestyle = '--')
        axes[m].fill_between(tvec_p, pup_av[:,m,3] + pup_sem[:,m,3], pup_av[:,m,3] - pup_sem[:,m,3], alpha = 0.2, color = var_color)
    
    
        axes[m].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
        axes[m].axvline(x = 2, linestyle = '--', color = 'gray', linewidth = 1)
    
        axes[m].set_ylabel(u"\u0394" + ' pupil size [au]')
        axes[m].set_xlim([-3,7])
        
        axes[m].yaxis.set_major_formatter('{x:,.2f}')
        ax2 = axes[m].twinx()
        ax2.set_yticks([])
        mouse_id = sniffs[m]['unique_id'][7:12]
        ax2.set_ylabel('Mouse ' + mouse_id)
        
        
    axes[m].legend()
    axes[m].set_xlabel('Time from odor presentation [sec]')
    fig.suptitle(incl_descr)
    
    plt.savefig(fig_path + '\\Pupil_' + incl_descr + '.png')

    #%% Plot habituation curve for pupil
    fig, axes = plt.subplots(int(nmice), 1, sharex = 'all', sharey='all', figsize = figure_size)
    axes = axes.flatten()
    
    pvec_f = np.arange(1,9)
    pvec_n = np.arange(1,9)
    
    for m in range(nmice):
        axes[m].errorbar(pvec_f, pup_1bin_av[:8,m,0], pup_1bin_sem[:8,m,0], label = 'FAM', color = fam_color, linestyle = '-', marker = 'o')
        axes[m].errorbar(pvec_n, pup_1bin_av[:8,m,1], pup_1bin_sem[:8,m,1], label = 'NOV', color = nov_color, linestyle = '-', marker = 'o')
        axes[m].errorbar(pvec_n, pup_1bin_av[:8,m,3], pup_1bin_sem[:8,m,3], label = var_label, color = var_color, linestyle = '--', marker = 'o')    
        
        axes[m].set_ylabel(u"\u0394" + ' pupil size [au]')
        
        ax2 = axes[m].twinx()
        ax2.set_yticks([])
        mouse_id = sniffs[m]['unique_id'][7:12]
        ax2.set_ylabel('Mouse ' + mouse_id)
        
        
    axes[m].legend()
    #axes[m].set_xticklabels(['PB1', 'PB2', 'PB4', '2', '4', '6', '8', '10'])
    axes[m].set_xlabel('Presentation number')
    fig.suptitle('Habituation curve')
    
    plt.savefig(fig_path + '\\Pupil_hab.png')

#%% Plot mean breathing for control (blank) trials with and wo. stimulation
ctrl_blanks = np.zeros(nmice)
opto_blanks = np.zeros(nmice)
av_ctrl = np.zeros([nmice, 2])

for m in range(nmice):    
    ctrl_blanks[m] = np.unique(sniffs[m]['trial_type'][np.logical_and(sniffs[m]['trial_opto']==0, sniffs[m]['trial_blank']==1)])[0]
    if odors_opto == odors_inh:
        opto_blanks[m] = np.unique(sniffs[m]['trial_type'][np.logical_and(sniffs[m]['trial_opto']==-1, sniffs[m]['trial_blank']==1)])[0]
    elif odors_opto == odors_stim:
        opto_blanks[m] = np.unique(sniffs[m]['trial_type'][np.logical_and(sniffs[m]['trial_opto']==1, sniffs[m]['trial_blank']==1)])[0]

exclude_pre = 4

fig, axes = plt.subplots(int(nmice), 1, sharex = 'all', sharey='all', figsize = (6, 15))
axes = axes.flatten()
tvec = np.linspace(-4, 8, nframes) - est_lat

for m in range(nmice):
    
    # Control blank
    which_rows = sniffs[m]['trial_type'] == ctrl_blanks[m]
    which_incl = sniffs[m]['trial_idx'][which_rows] - 1 # IN MATLAB TRIAL INDEXES START FROM 1!!
    tmp_data = sniff_delta[which_incl[exclude_pre:], :, m].T
    
    av_ctrl[m,0] = np.mean(np.mean(tmp_data[int(sniff_the_bin[0]*sr) : int(sniff_the_bin[1]*sr), :], 0))
    
    sniff_blank_ctrl_av = np.mean(tmp_data, 1)
    sniff_blank_ctrl_se = np.std(tmp_data, 1) / np.sqrt(tmp_data.shape[1])
    
    axes[m].plot(tvec, sniff_blank_ctrl_av, label = 'BLANK', color = 'C7')
    axes[m].fill_between(tvec, sniff_blank_ctrl_av + sniff_blank_ctrl_se, sniff_blank_ctrl_av - sniff_blank_ctrl_se, alpha = 0.2, color = 'C7')
    
    # Opto blank
    which_rows = sniffs[m]['trial_type'] == opto_blanks[m]
    which_incl = sniffs[m]['trial_idx'][which_rows] - 1 # IN MATLAB TRIAL INDEXES START FROM 1!!
    tmp_data = sniff_delta[which_incl, :, m].T
    
    av_ctrl[m,1] = np.mean(np.mean(tmp_data[int(sniff_the_bin[0]*sr) : int(sniff_the_bin[1]*sr), :], 0))
    
    sniff_blank_opto_av = np.mean(tmp_data, 1)
    sniff_blank_opto_se = np.std(tmp_data, 1) / np.sqrt(tmp_data.shape[1])
    
    axes[m].plot(tvec, sniff_blank_opto_av, label = var_label2, color = 'C5', linestyle = '--')
    axes[m].fill_between(tvec, sniff_blank_opto_av + sniff_blank_opto_se, sniff_blank_opto_av - sniff_blank_opto_se, alpha = 0.2, color = 'C5')
    
    # Visuals
    axes[m].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
    axes[m].axvline(x = 2, linestyle = '--', color = 'gray', linewidth = 1)

    axes[m].set_ylabel(u"\u0394" + ' sniffing [inh/sec]')
    #axes[m].set_ylabel(u"\u0394" + ' pupil [au]')
    
    ax2 = axes[m].twinx()
    axes[m].set_xlim([-3,7])
    ax2.set_yticks([])
    mouse_id = sniffs[m]['unique_id'][7:12]
    ax2.set_ylabel('Mouse ' + mouse_id)
    
    
axes[m].legend()
axes[m].set_xlabel('Time from odor presentation [sec]')
incl_descr = 'No odor trials'
fig.suptitle(incl_descr)

plt.savefig(fig_path + '\\Sniff_ctrl_' + incl_descr + '.png')
plt.savefig(fig_path + '\\Sniff_ctrl_' + incl_descr + '.svg', bbox_inches = 'tight', transparent = True)

    
#%% The same for pupil

if 'data_path' in locals():
    
    fig, axes = plt.subplots(int(nmice), 1, sharex = 'all', sharey='all', figsize = figure_size)
    axes = axes.flatten()
    tvec = np.linspace(-4, 8, pup_nframes) - est_lat
    
    for m in range(nmice):
        
        # Control blank
        which_rows = sniffs[m]['trial_type'] == ctrl_blanks[m]
        which_incl = sniffs[m]['trial_idx'][which_rows] - 1 # IN MATLAB TRIAL INDEXES START FROM 1!!
        tmp_data = pup_delta[which_incl[exclude_pre:], :, m].T
        
        pup_blank_ctrl_av = np.mean(tmp_data, 1)
        pup_blank_ctrl_se = np.std(tmp_data, 1) / np.sqrt(tmp_data.shape[1])
        
        axes[m].plot(tvec, pup_blank_ctrl_av, label = 'BLANK', color = 'C7')
        axes[m].fill_between(tvec, pup_blank_ctrl_av + pup_blank_ctrl_se, pup_blank_ctrl_av - pup_blank_ctrl_se, alpha = 0.2, color = 'C7')
        
        # Opto blank
        which_rows = sniffs[m]['trial_type'] == opto_blanks[m]
        which_incl = sniffs[m]['trial_idx'][which_rows] - 1 # IN MATLAB TRIAL INDEXES START FROM 1!!
        tmp_data = pup_delta[which_incl, :, m].T
        
        pup_blank_opto_av = np.mean(tmp_data, 1)
        pup_blank_opto_se = np.std(tmp_data, 1) / np.sqrt(tmp_data.shape[1])
        
        axes[m].plot(tvec, pup_blank_opto_av, label = var_label2, color = 'C5', linestyle = '--')
        axes[m].fill_between(tvec, pup_blank_opto_av + pup_blank_opto_se, pup_blank_opto_av - pup_blank_opto_se, alpha = 0.2, color = 'C5')
        
        # Visuals
        axes[m].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
        axes[m].axvline(x = 2, linestyle = '--', color = 'gray', linewidth = 1)
    
        axes[m].yaxis.set_major_formatter('{x:,.2f}')
        axes[m].set_ylabel(u"\u0394" + ' pupil [au]')
        axes[m].set_xlim([-3,7])
        
        ax2 = axes[m].twinx()
        ax2.set_yticks([])
        mouse_id = sniffs[m]['unique_id'][7:12]
        ax2.set_ylabel('Mouse ' + mouse_id)
        
        
    axes[m].legend()
    axes[m].set_xlabel('Time from odor presentation [sec]')
    incl_descr = 'No odor trials'
    fig.suptitle(incl_descr)
    
    plt.savefig(fig_path + '\\Pupil_ctrl_' + incl_descr + '.png')

#%% Preprocess licks data
lick_sr = 100 # the value to which you want to DOWNSAMPLE licking data
lick_ons = np.zeros([ntrials, 12*lick_sr, nses])
binsize = 0.5
bin_edges = np.arange(0, 12+binsize, binsize)
lick_hist = np.zeros([ntrials, len(bin_edges)-1, nses])

for m in range(nses):
    sniffs[m]['licks'] = [np.atleast_1d(x) for x in sniffs[m]['licks']]

    for tr in range(ntrials):
        tmp_idx = np.round(sniffs[m]['licks'][tr] * lick_sr).astype('int')
        tmp_idx = tmp_idx[tmp_idx < lick_ons.shape[1]]
        lick_ons[tr, tmp_idx, m] = 1 # code sniff onsets as 1
        lick_hist[tr,:,m] = np.histogram(sniffs[m]['licks'][tr], bin_edges)[0] / binsize

lick_gauss, lick_delta = st.ins_sniff(lick_ons, bsln_start = 1, odor_start = 4, sigma = 0.25, sr = lick_sr)

#%% Create reward vs. opto caegory matrix to loop through trials
rew_cat = []
rew_cat_names = ['SMALL REWARD', 'BIG REWARD', 'SMALL REWARD + INH', 'BIG REWARD + INH', \
                 'SMALL REWARD + STIM', 'BIG REWARD + STIM']
opto_labels = ['No laser', 'Inhibition', 'Stimulation']
    
for m in range(nmice):
    rew_cat.append(np.ones([ntrials, 6], np.int8))
    rew_cat[m][:,0] = np.logical_and(sniffs[m]['trial_reward'] == 1, sniffs[m]['trial_opto'] == 0)
    rew_cat[m][:,1] = np.logical_and(sniffs[m]['trial_reward'] == 4, sniffs[m]['trial_opto'] == 0)
    rew_cat[m][:,2] = np.logical_and(sniffs[m]['trial_reward'] == 1, sniffs[m]['trial_opto'] == -1)
    rew_cat[m][:,3] = np.logical_and(sniffs[m]['trial_reward'] == 4, sniffs[m]['trial_opto'] == -1)
    rew_cat[m][:,4] = np.logical_and(sniffs[m]['trial_reward'] == 1, sniffs[m]['trial_opto'] == 1)
    rew_cat[m][:,5] = np.logical_and(sniffs[m]['trial_reward'] == 4, sniffs[m]['trial_opto'] == 1)
    

#%% Plot licking & calculate proportion of anticipatory
anlicks = np.zeros([nmice, rew_cat[0].shape[1]])
tvec_licks = np.arange(-4, 8, 1/lick_sr)  

for m in range(nmice):
    fig, axes = plt.subplots(2, 2, sharex = 'all')
    axes[1,0].get_shared_y_axes().join(axes[1,0], axes[1,1])

    for c in range(np.shape(rew_cat[m])[1]):
        which_tr = sniffs[m]['trial_idx'][rew_cat[m][:,c]==True] - 1
        data_tmp = [sniffs[m]['licks'][i] -4 for i in which_tr]
        data_tmp2 = lick_gauss[which_tr,:,m]
        
        # Count proportion of trials with anticipatory licks
        data_tmp3 = [np.sum((x>0) & (x<3.8)) for x in data_tmp]
        data_tmp3 = np.sum(np.array(data_tmp3) > 0) / len(data_tmp3)
        anlicks[m,c] = data_tmp3
        
        if len(data_tmp) > 0:
            tmp2_mean = np.mean(data_tmp2, 0)
            tmp2_se = np.std(data_tmp2, 0) / np.sqrt(np.shape(data_tmp2)[0])
            
            col = c%2
            row = 0
            lbl_idx = int(np.round(c/2))
            if c > 1: # so, if opto
                row = 1
                ymin = 15
                use_color = 'C3'
            else: # if control
                ymin = 0
                use_color = 'C0'
                
            ymax = ymin + len(data_tmp)  
            ypos = np.arange(ymin, ymax) + 1
                
            axes[0,col].eventplot(data_tmp, color = use_color, lineoffsets = ypos, linewidth = 0.6)
            line, = axes[1, col].plot(tvec_licks, tmp2_mean, color = use_color)
            line.set_label(opto_labels[lbl_idx])
            axes[1, col].fill_between(tvec_licks, tmp2_mean + tmp2_se, tmp2_mean - tmp2_se, alpha = 0.2, color = use_color)

            
        axes[1, col].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
        axes[1, col].axvline(x = 4, linestyle = '--', color = 'gray', linewidth = 1)
        axes[0, col].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
        axes[0, col].axvline(x = 4, linestyle = '--', color = 'gray', linewidth = 1)
            
    axes[1,0].legend()
    axes[0,0].set_xlim([-3,7])
    axes[0,0].set_title('Small reward')
    axes[0,1].set_title('Big reward')
    axes[1,0].set_ylabel('Licks/sec')
    axes[0,0].set_ylabel('Trial')
    #axes[1,0].annotate('Odor', [0.1, 5])
    #axes[1,0].annotate('Reward', [4.1, 5])

    incl_descr = 'reward_licking'
    mouse_id = sniffs[m]['unique_id'][7:12]
    plt.savefig(fig_path + '\\' + incl_descr + '_' + mouse_id + '.png', dpi = 250)
    plt.savefig(fig_path + '\\' + incl_descr + '_' + mouse_id +  '.svg', bbox_inches = 'tight', transparent = True)


#%% Repeat, but for sniffing

tvec_licks = np.arange(-4, 8, 1/lick_sr)  
tvec_licks = np.linspace(-4, 8, nframes) 
for m in range(nmice):
    fig, axes = plt.subplots(2, 2, sharex = 'all')
    axes[1,0].get_shared_y_axes().join(axes[1,0], axes[1,1])

    for c in range(np.shape(rew_cat[m])[1]):
        which_tr = sniffs[m]['trial_idx'][rew_cat[m][:,c]==True] - 1
        data_tmp = [sniffs[m]['ml_inh_onsets'][i]/sr -4 for i in which_tr]
        data_tmp2 = sniff_delta[which_tr,:,m]
        
        if len(data_tmp) > 0:
            tmp2_mean = np.mean(data_tmp2, 0)
            tmp2_se = np.std(data_tmp2, 0) / np.sqrt(np.shape(data_tmp2)[0])
            
            col = c%2
            row = 0
            lbl_idx = int(np.round(c/2))
            if c > 1: # so, if opto
                row = 1
                ymin = 15
                use_color = 'C3'
            else: # if control
                ymin = 0
                use_color = 'C0'
                
            ymax = ymin + len(data_tmp)  
            ypos = np.arange(ymin, ymax) + 1
                
            axes[0,col].eventplot(data_tmp, color = use_color, lineoffsets = ypos, linewidth = 0.6)
            line, = axes[1, col].plot(tvec_licks, tmp2_mean, color = use_color)
            line.set_label(opto_labels[lbl_idx])
            axes[1, col].fill_between(tvec_licks, tmp2_mean + tmp2_se, tmp2_mean - tmp2_se, alpha = 0.2, color = use_color)

            
        axes[1, col].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
        axes[1, col].axvline(x = 4, linestyle = '--', color = 'gray', linewidth = 1)
        axes[0, col].axvline(x = 0, linestyle = '--', color = 'gray', linewidth = 1)
        axes[0, col].axvline(x = 4, linestyle = '--', color = 'gray', linewidth = 1)
            
    axes[1,0].legend()
    axes[0,0].set_xlim([-3,7])
    axes[0,0].set_title('Small reward')
    axes[0,1].set_title('Big reward')
    axes[1,0].set_ylabel('\u0394 sniffing [inh/sec]')
    axes[0,0].set_ylabel('Trial')
    axes[1,0].annotate('Odor', [0.1, 5])
    axes[1,0].annotate('Reward', [4.1, 5])

    incl_descr = 'reward_sniffing'
    mouse_id = sniffs[m]['unique_id'][7:12]
    plt.savefig(fig_path + '\\' + incl_descr + '_' + mouse_id + '.png', dpi = 250)
    plt.savefig(fig_path + '\\' + incl_descr + '_' + mouse_id +  '.svg', bbox_inches = 'tight', transparent = True)


#%% Plot effect summary (sniffing)

tmp_data = np.zeros([nmice, 3])
tmp_data[:,0] = np.mean(sniff_1bin_av[0:3,:,0],0)
tmp_data[:,1] = np.mean(sniff_1bin_av[0:3,:,1],0)
tmp_data[:,2] = np.mean(sniff_1bin_av[0:3,:,3],0)

grav = np.reshape(tmp_data, [np.size(tmp_data),],'F')
group = np.array([*[0]*nmice, *[1]*nmice, *[2]*nmice])
ind_color = [*[fam_color]*nmice, *[nov_color]*nmice, *[var_color]*nmice]
mfc_color = [*[fam_color]*nmice, *[nov_color]*nmice, *['white']*nmice]

#stim_color = [*[fam_color]*2, *['gray']*3]

plt.figure(figsize = (3.5,6), dpi = 250)
#plt.plot(tmp_data.T, color = 'gray', zorder=0)
for m in range(nmice):
    plt.plot(tmp_data[m,:], color = 'gray', zorder=0)
plt.scatter(group, grav, marker = 'o', color = ind_color, facecolors = mfc_color, s = 70, zorder=1)

plt.xticks([0,1,2], ['FAM', 'NOV', var_label], fontsize = 12, fontweight = 'bold')
plt.ylabel('\u0394 sniffing [inh/sec]')

#plt.ylim([-0.5,5.5])
plt.xlim([-0.3, 2.3])


# Perform statistical testing and annotate the results
line_h = np.max(tmp_data)+0.15 # how high to plot annotation lines
line_h = [line_h, line_h, line_h+0.25]

for ii, jj in enumerate([1, 2, 0]): # do pair-wise comparisons
    plt.plot([ii+0.05, jj-0.05], [line_h[ii], line_h[ii]], linewidth = 2.5, color = 'k')
    
    stat, p = stats.ttest_rel(grav[group==ii], grav[group==jj])
    p = p*3 # bonferroni correction
    if (p < 0.05) & (p > 0.01):
        sig_label = '*'; ss = 20; hh = line_h[ii]
    elif (p < 0.01) & (p > 0.001):
        sig_label = '**'; ss = 20; hh = line_h[ii]
    elif (p < 0.01) & (p > 0.001):
        sig_label = '***'; ss = 20; hh = line_h[ii]
    elif p > 0.05:
        sig_label = 'ns'; ss = 15; hh = line_h[ii] + 0.1

    plt.annotate(sig_label, [np.mean([ii, jj]), hh], size = ss, ha = 'center')
    print(p)

plt.savefig(fig_path + '\\sniff_effect_summary_' + '.png', bbox_inches = 'tight')
plt.savefig(fig_path + '\\sniff_effect_summary_' + '.svg', bbox_inches = 'tight', transparent = True)

#%% The same for control (blank) trials
plt.figure(figsize = (2,6), dpi = 250)
av_ctrl_vec = np.reshape(av_ctrl, [np.size(av_ctrl),],'F')
groupc = np.array([*[0]*nmice, *[1]*nmice])
ind_color = [*['k']*nmice*2]
mfc_color = [*['k']*nmice, *['white']*nmice]

stat, p = stats.ttest_rel(av_ctrl_vec[groupc==0], av_ctrl_vec[groupc==1])

if (p < 0.05) & (p > 0.01):
    sig_label = '*'; ss = 20
elif (p < 0.01) & (p > 0.001):
    sig_label = '**'; ss = 20
elif (p < 0.01) & (p > 0.001):
    sig_label = '***'; ss = 20
elif p > 0.05:
    sig_label = 'ns'; ss = 15

for m in range(nmice):
    plt.plot(av_ctrl[m,:], color = 'gray', zorder=0)
plt.scatter(groupc, av_ctrl_vec, marker = 'o', color = ind_color, facecolors = mfc_color, s = 70, zorder=1)

plt.xticks([0,1], ['CTRL', 'OPTO'], fontsize = 12, fontweight = 'bold')
plt.ylabel('\u0394 sniffing [inh/sec]')

line_h = np.max(av_ctrl_vec)+0.15 # how high to plot annotation lines
plt.plot([0.05, 0.95], [line_h, line_h], linewidth = 2.5, color = 'k')
plt.annotate(sig_label, [0.5, line_h+0.05], size = ss, ha = 'center')

#plt.ylim([-2, 2])
plt.xlim([-0.3, 1.3])

plt.savefig(fig_path + '\\sniff_ctrl_summary_' + '.png', bbox_inches = 'tight')
plt.savefig(fig_path + '\\sniff_ctrl_summary_' + '.svg', bbox_inches = 'tight', transparent = True)
