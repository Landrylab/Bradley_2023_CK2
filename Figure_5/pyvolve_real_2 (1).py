import os
import sys
import pyvolve
import numpy as np
import pandas as pd

# Read in the phylogeny
my_tree = pyvolve.read_tree(file=sys.argv[1],scale_tree = sys.argv[2])

# Define CustomFrequencies object (we had to make a slight modification to 'Y' so that the numbers sum exactly to 1)
f = pyvolve.CustomFrequencies("amino_acid", freq_dict = {"A":0.0828953, "C":0.00893229, "D":0.0536049, "E":0.0794287, "F":0.0247466, "G":0.0772864, "H":0.0217862, "I":0.0413678,"K":0.0848567, "L":0.0630974, "M":0.0197123, "N":0.0405072, "P":0.0714715, "Q":0.0498945, "R":0.0530623, "S":0.0810773, "T":0.0558998, "V":0.0585352, "W":0.00708994, "Y":0.0247476700000002})
                                                          
frequencies = f.compute_frequencies()

# Construct the custom matrix

disorder_mat = pd.read_csv('Anisimova_corrected_symmetrical_single_letter_order.txt', sep=" ")

# We need to convert the Panda DF into a NumPy array

custom_matrix = disorder_mat.to_numpy()

# Specify the custom model

custom_model = pyvolve.Model("custom", {"state_freqs": frequencies, "matrix":custom_matrix})

# Specify the ancestral sequence for your partition

partition = pyvolve.Partition(models = custom_model, root_sequence = sys.argv[3])


# Simulate evolution

my_evolver = pyvolve.Evolver(partitions = partition, tree = my_tree)

# Loop it

for i in range(100):
	my_evolver(seqfile = 'p'+sys.argv[4] + '_' + sys.argv[5] + '_' + str(i+100) + ".fasta") 
