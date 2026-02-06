# -*- coding: utf-8 -*-
"""
Quick script to adjust the timing of a experimental txt file for in labview
@author: Xander
"""

import pandas as pd

# Path to your file
file_path = r"C:\Users\xaand\OneDrive - KU Leuven\PhD\Experiments\Experiment Files\XT001-005\Day5_XT003_1.txt"
last_event = 408 #need to adjust this

# Load file (no header, first column as integers)
df = pd.read_csv(file_path, skiprows=1, header=None, sep="\s+")  # adjust sep if needed (e.g. ',' or '\t')


header = ['Delta T [ms]','16-bit pattern','Board ID,Valve ID']
time_to_remove = df.iloc[last_event-1,0]
# Subtract the constant from the first column
df[0] = df[0][1:].astype(int) - time_to_remove
df.columns = header
df = df.loc[last_event]
# Save to new file
output_path = file_path.replace(".txt", "_subtracted.txt")
df.to_csv(output_path, index=False, header=False, sep="\t")
print(f"Restart at time {time_to_remove}")
print(f"Saved result to: {output_path}")