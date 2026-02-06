# plot to perform and compare clustering performed by MANGEM (https://www.biorxiv.org/content/biorxiv/early/2023/04/04/2023.04.03.535322.full.pdf)

import pickle
import numpy as np
import pandas as pd
import sys
sys.path.append(
    r'C:\Users\xaand\Documents\PhD\Analysis\Informatics\PhD-patch-seq\Keimpe')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import umap


pickle_path = "C:/Users/xaand/Documents/PhD/Analysis/Informatics/PhD-patch-seq/data/ephys_data.pkl"

# loading the pickle into ps
with open(pickle_path, 'rb') as f:
    ps = pickle.load(f)
  
#%% exporting data
df_corr = ps["df"]
df_corr["Resting Vm"] = df_corr["Resting Vm"] *-1

crit = np.all(ps["selector"][["Sequencing","Rheobase","Resting Vm","Input Resistance","singles"]] ==1, axis =1)
#crit = np.all(ps["selector"][["Sequencing","Rheobase","Resting Vm","Input Resistance"]] ==1, axis =1)

df_corr.loc[crit,:].to_csv(r"C:\Users\xaand\Documents\PhD\Experiments\Patch seq OFC-LDTg\data\filtered csv\ephys_df_filtered.csv")
ps["counts"].loc[crit,:].to_csv(r"C:\Users\xaand\Documents\PhD\Experiments\Patch seq OFC-LDTg\data\filtered csv\counts_filtered.csv")
ps["metadata"].loc[crit,:].to_csv(r"C:\Users\xaand\Documents\PhD\Experiments\Patch seq OFC-LDTg\data\filtered csv\counts_metadata_filtered.csv")

print(df_corr.loc[crit,:],
      df_corr["Resting Vm"],
      pd.isna(df_corr.loc[crit,:]).sum()
        )



#%% Looking where values of 1 clustering lie in the other
path = r"C:\Users\xaand\Documents\PhD\Experiments\Patch seq OFC-LDTg\data\Clustering_Mangem"
ids = ["clusters_Counts_1.csv","clusters_Ephys df_1.csv","clusters_Counts_2.csv","clusters_Ephys df_2.csv"]
paths_list = [f"{path}\{sess_id}" for sess_id in ids]
clustering = []
for i in paths_list:
    tmp = pd.read_csv(i,index_col= "cell_id")
    tmp.drop(tmp.columns[[0,1]], axis = 1)
    clustering.append(tmp)

cell_list1 = df_corr.index[np.all(ps["selector"][["Sequencing","Rheobase","Resting Vm","Input Resistance","singles"]] ==1, axis =1)]
cell_list2 = df_corr.index[np.all(ps["selector"][["Sequencing","Rheobase","Resting Vm","Input Resistance"]] ==1, axis =1)]
#%% plotting dimensions with coloring on input or not

colors = pd.DataFrame(index= cell_list2)
df = pd.DataFrame({'Exists': [cell in cell_list1 for cell in cell_list2]})
plt.figure()
plt.scatter(clustering[2]['0'],clustering[2]['1'], color = "green")

plt.scatter(clustering[3]['0'],clustering[3]['1'], color = "purple")
plt.show()
#%%

bool_series = clustering[2].index.isin(cell_list1)
plt.figure()

df = pd.DataFrame({'Exists': [cell in cell_list1 for cell in cell_list2]})

# Scatter plot for clustering[2] with conditional coloring
plt.scatter(
#    x = clustering[2]['0'][bool_series], y =clustering[2]['1'][bool_series], 
            x= clustering[0]['0'], y= clustering[0]['1'],

#            c= clustering[2]['cluster'],
#            c=['green' if is_in else 'red' for is_in in bool_series_2]
            c="green" 
           )

# Scatter plot for clustering[3] with conditional coloring
plt.scatter(
#    x = clustering[3]['0'][bool_series], y = clustering[3]['1'][bool_series], 
    x= clustering[1]['0'], y= clustering[1]['1'],
#            c= clustering[3]['cluster']
#            c=['purple' if is_in else 'orange' for is_in in bool_series_3]
            c='purple'
            )

plt.show()

#%% TEST
crit = np.all(ps["selector"][["Sequencing","Rheobase","Resting Vm","Input Resistance"]] ==1, axis =1)

modality_1 = df_corr.loc[crit,:]
modality_2 = ps["counts"].loc[crit,:]

# Normalize each modality
scaler_1 = StandardScaler()
scaler_2 = StandardScaler()

modality_1_scaled = scaler_1.fit_transform(modality_1)
modality_2_scaled = scaler_2.fit_transform(modality_2)

# Merge into one dataset (for manifold alignment)
data_combined = np.hstack([modality_1_scaled, modality_2_scaled])

# Feature names
feature_names_1 = [f"Gene_{i+1}" for i in range(10)]
feature_names_2 = [f"Electro_{i+1}" for i in range(10)]
feature_names = feature_names_1 + feature_names_2

# UMAP for nonlinear alignment
umap_model = umap.UMAP(n_components=2, metric='euclidean', min_dist=0.3, n_neighbors=15)
embedding = umap_model.fit_transform(data_combined)  # Aligned representation

# Approximate feature importance in aligned space
loadings = umap_model.embedding_[:20]  # First 20 rows → Feature contributions


#Plot the Biplot
fig, ax = plt.subplots(figsize=(8, 6))

# Scatter plot of samples
ax.scatter(embedding[:, 0], embedding[:, 1], alpha=0.7, edgecolors='k', cmap="coolwarm")

# Add arrows for feature loadings
for i in range(10):  # Genes
    ax.arrow(0, 0, loadings[i, 0], loadings[i, 1], color='red', alpha=0.7, head_width=0.05)
    ax.text(loadings[i, 0] * 1.15, loadings[i, 1] * 1.15, feature_names[i], color='red', fontsize=10)

for i in range(10, 20):  # Electrophysiology
    ax.arrow(0, 0, loadings[i, 0], loadings[i, 1], color='blue', alpha=0.7, head_width=0.05)
    ax.text(loadings[i, 0] * 1.15, loadings[i, 1] * 1.15, feature_names[i], color='blue', fontsize=10)

# Formatting
ax.axhline(0, color='grey', linestyle='--', linewidth=0.5)
ax.axvline(0, color='grey', linestyle='--', linewidth=0.5)
ax.set_xlabel("UMAP1 (Aligned)")
ax.set_ylabel("UMAP2 (Aligned)")
ax.set_title("Cross-Modal Manifold Alignment Biplot")

plt.show()





