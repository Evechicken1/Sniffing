import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\xaand\OneDrive\unief\2023-2024\Masters_Thesis\Data_Analysis')
import sniff_tools as st
from scipy import stats

#%% Familiarization per mouse

#sess_ids = ['231020_KK066', '231215_KK082', '231215_KK085', '240211_KK092', '240211_KK093','240201_KK087', '240201_KK088', '240201_KK089', '240301_KK099', '240301_KK100'] # all mice inhibition
sess_ids = ['231024_KK066', '231219_KK082', '231219_KK085', '240215_KK092', '240215_KK093','240205_KK087', '240205_KK088', '240205_KK089', '240305_KK099', '240305_KK100'] # all mice stimulation
paths_list = [f"C:/Users/xaand/OneDrive/unief/2023-2024/Masters_Thesis/Data_Analysis/data/{sess_id}" for sess_id in sess_ids]
print(paths_list )
sniffs = st.import_sniff_mat_select(paths_list)
n_DATcre = 5

#%%

ntrials = sniffs[0]["trial_idx"].size
selector = [np.zeros((ntrials, 1), int) for i in range(len(sniffs))]   


for m in range(len(sniffs)):
    print("Mouse", m)
    isopto = (sniffs[m]["trial_opto"]== True) 
    isfam = (sniffs[m]["trial_familiarity"]== True) & (sniffs[m]["trial_opto"]== False) & (sniffs[m]["trial_occur"] <= 3)
    
    conditions = [isopto]
    
    trial_idx = sniffs[m]["trial_idx"]-1
    
    for cond in range(selector[m].shape[1]):
        idx_cnds = trial_idx[conditions[cond]]
        print(cond,idx_cnds)
        selector[m][idx_cnds,cond] = True
        print(selector[m][:,cond])
        
#%%
graph = ["DAT-Cre","Control"] #idx of condition must equal columnidx in selector[m] column from selector which you want to select , e.g. 0(int) for novel smells
sniff_lst = [[[]for m in range(n_DATcre)] for g in range(len(graph))] # make a list of arrays in size of selector sum of selected category
hist_edges = np.linspace(0, 723,13)




mean_arr = np.zeros((len(graph), n_DATcre,  12))
std_arr = np.zeros((len(graph), n_DATcre,  12))

for m in range(len(sniffs)):
    if m <=4:
        
        print("Mouse", m)
        print("graph", graph[0])
        print(paths_list[m])
        for idx in range(len(selector[m])):
            tmp = sniffs[m]["ml_inh_onsets"][idx]
            hist = np.histogram(tmp, bins=hist_edges)[0]
            cond1 = selector[m][idx,0] == 0

            if cond1 == True:
                tmp = sniffs[m]["ml_inh_onsets"][idx][(cond1)]
                hist = np.histogram(tmp, bins=hist_edges)[0]
#                    print(idx, hist)
                sniff_lst[0][m].append(hist)
#            print(idx, hist)
        mean_arr[0][m] = np.mean(sniff_lst[0][m], axis = 0)
        std_arr[0][m] = np.std(sniff_lst[0][m], axis = 0)


    else: 
        print("Mouse", m)
        print("graph", graph[1])
        print(paths_list[m])
        for idx in range(len(selector[m])):
            tmp = sniffs[m]["ml_inh_onsets"][idx]
            hist = np.histogram(tmp, bins=hist_edges)[0]
            cond1 = selector[m][idx,0] == 0
            if cond1 == True:
                tmp = sniffs[m]["ml_inh_onsets"][idx][(cond1)]
                hist = np.histogram(tmp, bins=hist_edges)[0]
#                    print(idx, hist)
                sniff_lst[1][m-n_DATcre].append(hist)
#            print(idx, hist)
        mean_arr[1][m-n_DATcre] = np.mean(sniff_lst[1][m-n_DATcre], axis = 0)
        std_arr[1][m-n_DATcre] = np.std(sniff_lst[1][m-n_DATcre], axis = 0)               
     
sniff_lst=np.array(sniff_lst)

#%%

bin_edges = np.linspace(-3.5, 7.5, num = 12)
colors = ['red' , 'black'] 
ecolor = ["red","black"]

sem = [[] for l in range(len(sniff_lst))]
mean = [[] for l in range(len(sniff_lst))]


plt.figure()
for g in range(len(graph)):
    print("graph", g)
    baseline = np.mean(mean_arr[g][:,0:3])
    sem[g] = np.std(mean_arr[g], axis = 0) / np.sqrt(n_DATcre) #calc SEM     
    mean[g] = np.mean(mean_arr[g], axis=0) # - baseline

    plt.plot(bin_edges, mean[g], color= colors[g], label = graph[g], ls = ["solid", "solid", "dashed"][g])
    plt.title("All trials - opto")
    plt.legend(loc = "lower right")
    plt.errorbar(bin_edges, mean[g], yerr=sem[g], fmt='o', color=colors[g], ls = ["solid", "solid", "dashed"][g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
    plt.axvspan(xmin=0, xmax=2, color="grey", alpha=0.1)
    plt.ylabel(" sniffs")
    plt.xlabel("time(s)")
    ax=plt.gca() #to remve the top and right spines
    ax.spines[['right', 'top']].set_visible(False) #to remve the top and right spines
    
plt.show()

#%% wilcoxon signed rank test

for bin in np.arange(0,12):
    if (stats.shapiro(mean_arr[0][:,bin])[1] >= 0.05) and (stats.shapiro(mean_arr[1][:,bin])[1] >= 0.05) and (stats.ttest_ind(mean_arr[0][:,bin], mean_arr[1][:,bin])[1] <=0.05):
           
        print(bin, stats.ttest_ind(mean_arr[0][:,bin], mean_arr[1][:,bin], ))
        print (stats.shapiro(mean_arr[0][:,bin])[1], stats.shapiro(mean_arr[1][:,bin])[1])
        
#%%

