import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'C:\Users\xaand\OneDrive\unief\2023-2024\Masters_Thesis\Data_Analysis')
import sniff_tools as st

#%% Familiarization per mouse

sess_ids = ['231020_KK066', '231215_KK082', '231215_KK085', '240211_KK092', '240211_KK093','240201_KK087', '240201_KK088', '240201_KK089', '240301_KK099', '240301_KK100'] # all mice

paths_list = [f"C:/Users/xaand/OneDrive/unief/2023-2024/Masters_Thesis/Data_Analysis/data/{sess_id}" for sess_id in sess_ids]
print(paths_list )
sniffs = st.import_sniff_mat_select(paths_list)
n_DATcre = 5

#%%

ntrials = sniffs[0]["trial_idx"].size
selector = [np.zeros((ntrials, 1), int) for i in range(len(sniffs))]   


for m in range(len(sniffs)):
    print("Mouse", m)
    isnov = (sniffs[m]["trial_novelty"]== True) & (sniffs[m]["trial_opto"]== False) & (sniffs[m]["trial_occur"] <= 3)
    isfam = (sniffs[m]["trial_familiarity"]== True) & (sniffs[m]["trial_opto"]== False) & (sniffs[m]["trial_occur"] <= 3)
    
    conditions = [isnov]
    
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
            cond1 = selector[m][idx,0] == 1

            if cond1 == True:
                tmp = sniffs[m]["ml_inh_onsets"][idx][(cond1)]
                hist = np.histogram(tmp, bins=hist_edges)[0]
#                    print(idx, hist)
                sniff_lst[0][m].append(hist)
        mean_arr[0][m] = np.mean(sniff_lst[0][m], axis = 0)
        std_arr[0][m] = np.std(sniff_lst[0][m], axis = 0)


    else: 
        print("Mouse", m)
        print("graph", graph[1])
        print(paths_list[m])
        for idx in range(len(selector[m])):
            cond1 = selector[m][idx,0] == 1

            if cond1 == True:
                tmp = sniffs[m]["ml_inh_onsets"][idx][(cond1)]
                hist = np.histogram(tmp, bins=hist_edges)[0]
                print(idx, hist)
                sniff_lst[1][m-n_DATcre].append(hist)
        mean_arr[1][m-n_DATcre] = np.mean(sniff_lst[1][m-n_DATcre], axis = 0)
        std_arr[1][m-n_DATcre] = np.std(sniff_lst[1][m-n_DATcre], axis = 0)               
      
sniff_lst=np.array(sniff_lst)

#%%

bin_edges = np.linspace(-3.5, 7.5, num = 12)
colors = ['red' , 'black'] 
ecolor = ["darkred","black"]

sem = [[] for l in range(len(sniff_lst))]
mean = [[] for l in range(len(sniff_lst))]

#%%
plt.figure()
for g in range(len(graph)):
    print("graph", g)
    baseline = (np.mean(mean_arr[g][:,0:3]))
    sem[g] = np.std(mean_arr[g], axis = 0) / np.sqrt(n_DATcre) #calc SEM     
    mean[g] = np.mean(mean_arr[g], axis=0) - baseline

    plt.plot(bin_edges, mean[g], color= colors[g], label = graph[g], ls = ["solid", "solid", "dashed"][g])
    plt.title("First 3 novel trials")
    plt.legend()
    plt.errorbar(bin_edges, mean[g], yerr=sem[g], fmt='o', color=colors[g], ls = ["solid", "solid", "dashed"][g] , ecolor= ecolor[g] , elinewidth=1, capsize=3)
    plt.axvspan(xmin=0, xmax=2, color="grey", alpha=0.1)
    plt.axhline(y=0,color="black",linewidth=1, alpha = 0.1, ls = "dotted")
    plt.ylabel(" Δ avg sniffs")
    plt.xlabel("time(s)")
    ax=plt.gca() #to remve the top and right spines
    ax.spines[['right', 'top']].set_visible(False) #to remve the top and right spines
    
plt.show()

#%%
graph = ["DAT-Cre","Control"]
#%% plot comparing sniffing to novel between DATcre and control

tmp = np.zeros(((n_DATcre*18)), dtype= int)
data = np.zeros(( n_DATcre*18,2))


colors = ['green' , 'purple'] 
ecolor = ["darkgreen","purple"]
BoI = np.linspace(0,12,13,dtype=int)[5:8]

plt.figure()
for g in range(len(graph)):
    print(g)
    tmp = []
    for m in range(n_DATcre):
        print(m)
        
        tmp[slice(18*m,18*(m+1))] = np.mean(sniff_lst[g][m][:,BoI], axis=1)
    data[:,g] = tmp
    
    plt.boxplot(data, notch= True)
plt.legend(graph)

plt.show()