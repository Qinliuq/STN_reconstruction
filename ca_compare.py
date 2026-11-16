import pandas as pd
import numpy as np

# Load the files
base_file = pd.read_csv('cells/sth-data/cell-gcaL_HVA', sep=' ', header=None)
dl09_file = pd.read_csv('cells/sth-data/cell-gcaL_HVA-dl0.9', sep=' ', header=None) 
dl11_file = pd.read_csv('cells/sth-data/cell-gcaL_HVA-dl1.1', sep=' ', header=None)

# Get values column (last column)
base_vals = base_file.iloc[:,-1]
dl09_vals = dl09_file.iloc[:,-1]
dl11_vals = dl11_file.iloc[:,-1]

# Compare soma values (first row)
print("Soma comparisons:")
print(f"Base soma value: {base_vals[0]:.6f}")
print(f"dl0.9 soma value: {dl09_vals[0]:.6f} (ratio: {dl09_vals[0]/base_vals[0]:.3f}x)")
print(f"dl1.1 soma value: {dl11_vals[0]:.6f} (ratio: {dl11_vals[0]/base_vals[0]:.3f}x)")

# Compare dendrite values (element-wise ratios then average)
print("\nDendrite comparisons (element-wise ratios then average):")
dl09_ratios = dl09_vals[1:] / base_vals[1:]
dl11_ratios = dl11_vals[1:] / base_vals[1:]
print(f"dl0.9 average ratio: {dl09_ratios.mean():.3f}x")
print(f"dl1.1 average ratio: {dl11_ratios.mean():.3f}x")

locs_vals = zip(base_file.iloc[:,2].values, dl09_ratios.values)
# print as list of tuples
print(list(locs_vals))
