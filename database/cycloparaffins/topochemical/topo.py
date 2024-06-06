import csv

from padelpy import from_smiles
from rdkit import Chem
from rdkit.Chem import AllChem

def has_duplicates(seq):
    return len(seq) != len(set(seq))
    
def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

# ++++++ C6 cycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCCCC1', 'C[C@@H]1CC[C@H]1C', 'CC1(C)CCC1', 'C[C@@H]1CC[C@@H]1C', 'C[C@H]1C[C@@H](C)C1', 'CC1CC(C)C1', 'CCC1CCC1', 'CC(C)C1CC1', 'CCC1(C)CC1', 'CC[C@H]1C[C@H]1C', 'CC1CCCC1', 'CC[C@@H]1C[C@H]1C', 'CCCC1CC1', 'CC1CC1(C)C', 'CC1C(C1C)C', 'CC1C(C)C1C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cyclohexane', '1,trans-2-dimethylcyclobutane', '1,1-dimethylcyclobutane', '1,cis-2-dimethylcyclobutane', '1,cis-3-dimethylcyclobutane', '1,trans-3-dimethylcyclobutane', 'ethylcyclobutane', 'isopropylcyclopropane', '1-methyl-1-ethylcyclopropane', '1-methyl-cis-2-ethylcyclopropane', 'methylcyclopentane', '1-methyl-trans-2-ethylcyclopropane', 'propylcyclopropane', '1,1,2-trimethylcyclopropane', '1,cis-2,cis-3-trimethylcyclopropane', '1,cis-2,trans-3-trimethylcyclopropane']

eta_star = []

for i in range(len(descriptors)):
	
	eta_B = float(descriptors[i]['ETA_EtaP_B_RC'])
				
	# Add explicit hydrogens
	mol = Chem.MolFromSmiles(smiles_list[i])
	mol = Chem.AddHs(mol)
	# Generate 3D conformers
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	# Calculate molecular volume
	vol = AllChem.ComputeMolVolume(mol,gridSpacing=.05)
				
	eta_star.append(eta_B*vol)

if has_duplicates(eta_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_star):
		if eta_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))

# name of csv file
filename = "C6.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C7 cycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCCCCC1', 'C[C@@H]1CCC[C@@H]1C', 'C[C@@H]1CC[C@@H](C1)C', 'C[C@@H]1CC[C@H](C1)C', 'C[C@@H]1CCC[C@H]1C', 'CC1(C)CCCC1', 'CCC1CCCC1', 'CC1CCCCC1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cycloheptane', 'cis-1,2-dimethylcyclopentane', 'cis-1,3-dimethylcyclopentane', 'trans-1,3-dimethylcyclopentane', 'trans-1,2-dimethylcyclopentane', '1,1-dimethylcyclopentane', 'ethylcyclopentane', 'methylcyclohexane']

eta_star = []

for i in range(len(descriptors)):
	
	eta_B = float(descriptors[i]['ETA_EtaP_B_RC'])
				
	# Add explicit hydrogens
	mol = Chem.MolFromSmiles(smiles_list[i])
	mol = Chem.AddHs(mol)
	# Generate 3D conformers
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	# Calculate molecular volume
	vol = AllChem.ComputeMolVolume(mol,gridSpacing=.05)
				
	eta_star.append(eta_B*vol)

if has_duplicates(eta_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_star):
		if eta_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))

# name of csv file
filename = "C7.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C8 cycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCCCCCC1', 'C[C@H]1CC[C@@H](C)CC1', 'CC1(C)CCCCC1', 'C[C@@H]1CCCC[C@@H]1C', 'C[C@@H]1CCCC[C@H]1C', 'C[C@@H]1CCC[C@H](C1)C', 'C[C@H]1CC[C@@H](CC1)C', 'C[C@@H]1CCC[C@@H](C1)C', 'CCC1CCCCC1', 'CC(C)C1CCCC1', 'CCC1(C)CCCC1', 'CC[C@@H]1CCC[C@@H]1C', 'CC[C@H]1CC[C@@H](C)C1', 'CC[C@H]1CCC[C@@H]1C', 'CC[C@H]1CC[C@H](C)C1', 'CCCC1CCCC1', 'C[C@@H]1CC(C[C@@H]1C)C', 'C[C@H]1CC[C@@H](C)[C@@H]1C', 'CC1CCC(C)(C)C1', 'C[C@@H]1CC[C@@H](C)C1C', 'CC1CCCC1(C)C', 'CC1C[C@@H](C)[C@H](C)C1', 'C[C@H]1CC[C@@H](C)[C@H]1C', 'C[C@H]1C[C@@H](C)[C@@H](C)C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cyclooctane', 'cis-1,4-dimethylcyclohexane', '1,1-dimethylcyclohexane', 'cis-1,2-dimethylcyclohexane', 'trans-1,2-dimethylcyclohexane', 'trans-1,3-dimethylcyclohexane', 'trans-1,4-dimethylcyclohexane', 'cis-1,3-dimethylcyclohexane', 'ethylcyclohexane', 'isopropylcyclopentane', '1-methyl-1-ethylcyclopentane', '1-methyl-cis-2-ethylcyclopentane', '1-methyl-cis-3-ethylcyclopentane', '1-methyl-trans-2-ethylcyclopentane', '1-methyl-trans-3-ethylcyclopentane', 'propylcyclopentane', '1,cis-2,trans-4-trimethylcyclopentane', '1,trans-2,cis-3-trimethylcyclopentane', '1,1,3-trimethylcyclopentane', '1,cis-2,trans-3-trimethylcyclopentane', '1,1,2-trimethylcyclopentane', '1,trans-2,cis-4-trimethylcyclopentane', '1,cis-2,cis-3-trimethylcyclopentane', '1,cis-2,cis-4-trimethylcyclopentane']

eta_star = []

for i in range(len(descriptors)):
	
	eta_B = float(descriptors[i]['ETA_EtaP_B_RC'])
				
	# Add explicit hydrogens
	mol = Chem.MolFromSmiles(smiles_list[i])
	mol = Chem.AddHs(mol)
	# Generate 3D conformers
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	# Calculate molecular volume
	vol = AllChem.ComputeMolVolume(mol,gridSpacing=.05)
				
	eta_star.append(eta_B*vol)

if has_duplicates(eta_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_star):
		if eta_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))

# name of csv file
filename = "C8.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C9 cycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCCCCCCC1', 'CC(C)(C)C1CCCC1', 'CCCCC1CCCC1', 'CCC(C)C1CCCC1', 'CC[C@H]1CCC[C@H]1CC', 'CCC1(CC)CCCC1', 'CC[C@H]1CC[C@H](CC)C1', 'CC[C@H]1CCC[C@@H]1CC', 'CC[C@H]1CC[C@@H](CC)C1', 'CC[C@]1(C)CC[C@H](C)C1', 'CC[C@]1(C)CC[C@@H](C)C1', 'CC[C@]1(C)CCC[C@@H]1C', 'CC[C@]1(C)CCC[C@H]1C', 'CCC1CCCC1(C)C', 'CCC1CCC(C)(C)C1', 'CCC1[C@H](C)CC[C@H]1C', 'CC[C@H]1[C@@H](C)CC[C@H]1C', 'CC[C@H]1CC[C@@H](C)[C@H]1C', 'CC[C@@H]1CC[C@@H]([C@H]1C)C', 'CC[C@@H]1C[C@H](C)[C@H](C)C1', 'CC[C@@H]1C[C@H](C)C[C@H]1C', 'CC[C@H]1C[C@@H](C)C[C@H]1C', 'CCC1C[C@@H](C)[C@H](C)C1', '[C@@H]1([C@@H]([C@H](CC1)C)CC)C', 'CC[C@H]1CC[C@@H]([C@H]1C)C', 'CC[C@H]1CC[C@H](C)[C@H]1C', 'CC[C@H]1C[C@H](C)[C@H](C)C1', 'CC[C@@H]1C[C@@H](C)C[C@H]1C', 'CC[C@@H]1C[C@@H](C)C[C@@H]1C', 'CC(C)CC1CCCC1', 'CC(C)C1CCCCC1', 'CCC1(C)CCCCC1', 'CC(C)C1(C)CCCC1', 'CCCC1(C)CCCC1', 'CC[C@H]1CCCC[C@H]1C', 'CC(C)[C@H]1CCC[C@H]1C', 'CCC[C@@H]1CCC[C@@H]1C', 'CC[C@@H]1CCC[C@H](C)C1', 'C[C@H]1CC[C@H](C1)C(C)C', 'CCC[C@H]1CC[C@@H](C)C1', 'CC[C@@H]1CC[C@H](C)CC1', 'CC[C@H]1CCCC[C@@H]1C', 'CC(C)[C@@H]1CCC[C@H]1C', 'CCC[C@H]1CCC[C@@H]1C', 'CC[C@H]1CCC[C@H](C)C1', 'CC(C)[C@@H]1CC[C@@H](C)C1', 'CCC[C@H]1CC[C@H](C)C1', 'CC[C@H]1CC[C@H](C)CC1', 'CCCC1CCCCC1', '[C@@H]1([C@H]([C@@H]([C@@H](C1)C)C)C)C', '[C@@H]1([C@H]([C@@H]([C@H](C1)C)C)C)C', 'C[C@@H]1C[C@H](C)[C@@H](C)[C@@H]1C', 'C[C@@H]1C[C@@H](C)[C@H](C)[C@@H]1C', 'C1(C[C@H]([C@@H](C1)C)C)(C)C', '[C@@H]1(C([C@@H](CC1)C)(C)C)C', 'C1([C@@H](C[C@@H](C1)C)C)(C)C', 'CC1(C)CCC(C)(C)C1', 'C1(C[C@@H]([C@@H](C1)C)C)(C)C', '[C@@H]1(C([C@H](CC1)C)(C)C)C', 'C1([C@@H]([C@@H](CC1)C)C)(C)C', 'C[C@@H]1CCC(C)(C)[C@H]1C', 'CC1(C)CCCC1(C)C', 'C[C@H]1C[C@@H](C)[C@H](C)[C@@H]1C', 'C[C@H]1C[C@H](C)[C@@H](C)[C@@H]1C', '[C@@H]1([C@@H](C[C@@H](CC1)C)C)C', '[C@@H]1([C@H](C[C@@H](CC1)C)C)C', '[C@@H]1([C@H](CCC[C@H]1C)C)C', '[C@@H]1(C[C@@H](C[C@H](C1)C)C)C', 'C[C@H]1C[C@H](C[C@H](C1)C)C', 'C[C@@H]1CC[C@H](C)[C@@H](C)C1', 'C[C@H]1C[C@@H](C[C@@H](C1)C)C', 'CC1CCCC(C)(C)C1', 'C[C@H]1CCC[C@H](C)C1C', 'CC1CCC(C)(C)CC1', 'CC1CCCCC1(C)C', 'C[C@@H]1CC[C@H](C)[C@H](C)C1', 'C[C@H]1CCC[C@@H](C)[C@@H]1C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cyclononane', 'tert-butylcyclopentane', 'butylcyclopentane', 'sec-butylcyclopentane', '1,cis-2-diethylcyclopentane', '1,1-diethylcyclopentane', '1,trans-3-diethylcyclopentane', '1,trans-2-diethylcyclopentane', '1,cis-3-diethylcyclopentane', '1,cis-3-dimethyl-1-ethylcyclopentane', '1,trans-3-dimethyl-1-ethylcyclopentane', '1,cis-2-dimethyl-1-ethylcyclopentane', '1,trans-2-dimethyl-1-ethylcyclopentane', '1,1-dimethyl-2-ethylcyclopentane', '1,1-dimethyl-3-ethylcyclopentane', '1,trans-3-dimethyl-cis-2-ethylcyclopentane', '1,cis-3-dimethyl-cis-2-ethylcyclopentane', '1,cis-2-dimethyl-cis-3-ethylcyclopentane', '1,trans-2-dimethyl-cis-3-ethylcyclopentane', '1,cis-2-dimethyl-cis-4-ethylcyclopentane', '1,trans-3-dimethyl-cis-4-ethylcyclopentane', '1,cis-3-dimethyl-cis-4-ethylcyclopentane', '1,trans-2-dimethyl-cis-4-ethylcyclopentane', '1,cis-3-dimethyl-trans-2-ethylcyclopentane', '1,trans-2-dimethyl-trans-3-ethylcyclopentane', '1,cis-2-dimethyl-trans-3-ethylcyclopentane', '1,cis-2-dimethyl-trans-4-ethylcyclopentane', '1,cis-3-dimethyl-trans-4-ethylcyclopentane', '1,trans-3-dimethyl-trans-4-ethylcyclopentane', 'isobutylcyclopentane', 'isopropylcyclohexane', '1-methyl-1-ethylcyclohexane', '1-methyl-1-isopropylcyclopentane', '1-methyl-1-propylcyclopentane', '1-methyl-cis-2-ethylcyclohexane', '1-methyl-cis-2-isopropylcyclopentane', '1-methyl-cis-2-propylcyclopentane', '1-methyl-cis-3-ethylcyclohexane', '1-methyl-cis-3-isopropylcyclopentane', '1-methyl-cis-3-propylcyclopentane', '1-methyl-cis-4-ethylcyclohexane', '1-methyl-trans-2-ethylcyclohexane', '1-methyl-trans-2-isopropylcyclopentane', '1-methyl-trans-2-propylcyclopentane', '1-methyl-trans-3-ethylcyclohexane', '1-methyl-trans-3-propylcyclopentane', '1-methyl-trans-3-isopropylcyclopentane', '1-methyl-trans-4-ethylcyclohexane', 'propylcyclohexane', '1,cis-2,cis-3,cis-4-tetramethylcyclopentane', '1,cis-2,trans-3,trans-4-tetramethylcyclopentane', '1,cis-2,trans-3,cis-4-tetramethylcyclopentane', '1,cis-2,cis-3,trans-4-tetramethylcyclopentane', '1,1,cis-3,cis-4-tetramethylcyclopentane', '1,2,2,trans-3-tetramethylcyclopentane', '1,1,cis-2,cis-4-tetramethylcyclopentane', '1,1,3,3-tetramethylcyclopentane', '1,1,cis-3,trans-4-tetramethylcyclopentane', '1,2,2,cis-3-tetramethylcyclopentane', '1,1,cis-2,trans-3-tetramethylcyclopentane', '1,1,cis-2,cis-3-tetramethylcyclopentane', '1,1,2,2-tetramethylcyclopentane', '1,trans-2,trans-3,cis-4-tetramethylcyclopentane', '1,trans-2,cis-3,trans-4-tetramethylcyclopentane', '1,trans-2,trans-4-trimethylcyclohexane', '1,cis-2,trans-4-trimethylcyclohexane', '1,cis-2,cis-3-trimethylcyclohexane', 'cis,trans-1,3,5-trimethylcyclohexane', '1,cis-3,cis-5-trimethylcyclohexane', '1,trans-2,cis-4-trimethylcyclohexane', '1,cis-3,trans-5-trimethylcyclohexane', '1,1,3-trimethylcyclohexane', '1,cis-2,trans-3-trimethylcyclohexane', '1,1,4-trimethylcyclohexane', '1,1,2-trimethylcyclohexane', '1,cis-2,cis-4-trimethylcyclohexane', '1,trans-2,cis-3-trimethylcyclohexane']

eta_star = []

for i in range(len(descriptors)):
	
	eta_B = float(descriptors[i]['ETA_EtaP_B_RC'])
				
	# Add explicit hydrogens
	mol = Chem.MolFromSmiles(smiles_list[i])
	mol = Chem.AddHs(mol)
	# Generate 3D conformers
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	# Calculate molecular volume
	vol = AllChem.ComputeMolVolume(mol,gridSpacing=.05)
				
	eta_star.append(eta_B*vol)

if has_duplicates(eta_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_star):
		if eta_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))

# name of csv file
filename = "C9.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C10 cycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCCCCCCCC1', 'CCCCC1CCCCC1', 'CCCCCC1CCCC1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cyclodecane', 'butylcyclohexane', 'pentylcyclopentane']

eta_star = []

for i in range(len(descriptors)):
	
	eta_B = float(descriptors[i]['ETA_EtaP_B_RC'])
				
	# Add explicit hydrogens
	mol = Chem.MolFromSmiles(smiles_list[i])
	mol = Chem.AddHs(mol)
	# Generate 3D conformers
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	# Calculate molecular volume
	vol = AllChem.ComputeMolVolume(mol,gridSpacing=.05)
				
	eta_star.append(eta_B*vol)

if has_duplicates(eta_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_star):
		if eta_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))

# name of csv file
filename = "C10.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C11 cycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCCCCCCCCC1', 'CCCCCC1CCCCC1', 'CCCCCCC1CCCC1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cycloundecane', 'pentylcyclohexane', 'hexylcyclopentane']

eta_star = []

for i in range(len(descriptors)):
	
	eta_B = float(descriptors[i]['ETA_EtaP_B_RC'])
				
	# Add explicit hydrogens
	mol = Chem.MolFromSmiles(smiles_list[i])
	mol = Chem.AddHs(mol)
	# Generate 3D conformers
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	# Calculate molecular volume
	vol = AllChem.ComputeMolVolume(mol,gridSpacing=.05)
				
	eta_star.append(eta_B*vol)

if has_duplicates(eta_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_star):
		if eta_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))

# name of csv file
filename = "C11.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C12 cycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCCCCCCCCCC1', 'CCCCCCCC1CCCC1', 'CCCCCCC1CCCCC1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cyclododecane', 'heptylcyclopentane', 'hexylcyclohexane']

eta_star = []

for i in range(len(descriptors)):
	
	eta_B = float(descriptors[i]['ETA_EtaP_B_RC'])
				
	# Add explicit hydrogens
	mol = Chem.MolFromSmiles(smiles_list[i])
	mol = Chem.AddHs(mol)
	# Generate 3D conformers
	AllChem.EmbedMolecule(mol)
	AllChem.UFFOptimizeMolecule(mol)
	# Calculate molecular volume
	vol = AllChem.ComputeMolVolume(mol,gridSpacing=.05)
				
	eta_star.append(eta_B*vol)

if has_duplicates(eta_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_star):
		if eta_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))

# name of csv file
filename = "C12.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)