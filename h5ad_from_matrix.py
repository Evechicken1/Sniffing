#this script will make an annotated file of a column of genes
#the columns are gene_symbol	gene_id	chromosome	gene_entrez_id	gene_name

from Bio import Entrez
import pandas as pd
import sys
import sqlalchemy
import cruzdb
from cruzdb import Genome


#%%

gene_name = "Gm1992"

Entrez.email = 'xander.temmerman@nerf.be'

query = f'{gene_name}[Gene Name]'
handle = Entrez.esearch(db='gene', term=query)
record = Entrez.read(handle)
handle.close()

#%%

path = ["C:/Users/xaand/Documents/PhD/Experiments/Patch seq OFC-LDTg/data/All_counts.csv"]
gene_list = pd.read_csv(path[0]
                     , usecols = [0])

#%%

Entrez.email = 'xander.temmerman@nerf.be'
Entrez_ids = []

#query = [f'{gene_list.iloc[i,0]}[Gene Name]' for i in range(len(gene_list))]
test = [f'{gene_list.iloc[i,0]}[Gene Name] AND Mus musculus[Organism]' for i in range(10)]

for i in range(10):
        
    handle = Entrez.esearch(db = "gene", term = test[i])
    record = Entrez.read(handle)
    Entrez_ids.append(record["IdList"])
    
    handle.close()

