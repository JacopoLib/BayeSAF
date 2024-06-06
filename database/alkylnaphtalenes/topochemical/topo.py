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

# ++++++ C10 alkylnaphtalenes ++++++

# list of SMILES
smiles_list = ['C1=CC2=C(C=C1)C=CC=C2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['naphthalene']

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
    
# ++++++ C11 alkylnaphtalenes ++++++

# list of SMILES
smiles_list = ['CC1=CC2=C(C=CC=C2)C=C1', 'CC1=C2C=CC=CC2=CC=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2-methylnaphthalene', '1-methylnaphthalene']

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
    
# ++++++ C12 alkylnaphtalenes ++++++

# list of SMILES
smiles_list = ['CC1=CC=CC2=C(C)C=CC=C12', 'CC1=CC2=C(C=C1)C=CC(C)=C2', 'CC1=CC2=C(C=C1)C=C(C)C=C2', 'CC1=C(C)C=C2C=CC=CC2=C1', 'CC1=CC=CC2=CC=CC(C)=C12', 'CC1=CC2=C(C=C1)C(C)=CC=C2', 'CC1=C2C=CC=CC2=C(C)C=C1', 'CC1=CC2=C(C=CC=C2)C(C)=C1', 'CC1=CC2=C(C=CC=C2C)C=C1', 'CC1=C(C)C2=C(C=CC=C2)C=C1', 'CCC1=CC2=C(C=CC=C2)C=C1', 'CCC1=C2C=CC=CC2=CC=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1,5-dimethylnaphthalene', '2,7-dimethylnaphthalene', '2,6-dimethylnaphthalene', '2,3-dimethylnaphthalene', '1,8-dimethylnaphthalene', '1,6-dimethylnaphthalene', '1,4-dimethylnaphthalene', '1,3-dimethylnaphthalene', '1,7-dimethylnaphthalene', '1,2-dimethylnaphthalene', '2-ethylnaphthalene', '1-ethylnaphthalene']

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

# ++++++ C13 alkylnaphtalenes ++++++

# list of SMILES
smiles_list = ['CC(C)C1=C2C=CC=CC2=CC=C1', 'CC(C)C1=CC2=C(C=CC=C2)C=C1', 'CCC1=C2C=CC=CC2=CC=C1C', 'CCC1=C(C)C2=C(C=CC=C2)C=C1', 'CCC1=C(C)C=C2C=CC=CC2=C1', 'CCC1=CC2=C(C=CC=C2)C(C)=C1', 'CCC1=C2C=CC=CC2=C(C)C=C1', 'CCC1=C2C=CC=CC2=CC(C)=C1', 'CCC1=C2C=CC=C(C)C2=CC=C1', 'CCC1=C2C=CC(C)=CC2=CC=C1', 'CCC1=CC2=C(C=C1)C=C(C)C=C2', 'CCC1=CC2=C(C=C1)C(C)=CC=C2', 'CCC1=CC2=C(C=CC=C2C)C=C1', 'CCC1=CC2=C(C=CC(C)=C2)C=C1', 'CCC1=C2C=C(C)C=CC2=CC=C1', 'CCC1=C2C(C)=CC=CC2=CC=C1', 'CCCC1=CC2=C(C=CC=C2)C=C1', 'CCCC1=C2C=CC=CC2=CC=C1', 'CC1=CC2=CC=C(C)C(C)=C2C=C1', 'CC1=CC2=C(C=C1)C(C)=CC=C2C', 'CC1=C(C)C(C)=C2C=CC=CC2=C1', 'CC1=CC2=C(C=C1)C=C(C)C(C)=C2', 'CC1=CC(C)=C(C)C2=C1C=CC=C2', 'CC1=CC=CC2=C1C=CC(C)=C2C', 'CC1=CC=CC2=CC(C)=C(C)C=C12', 'CC1=CC2=C(C=C1)C=CC(C)=C2C', 'CC1=CC=CC2=C1C(C)=C(C)C=C2', 'CC1=CC2=C(C=CC=C2C)C(C)=C1', 'CC1=CC2=C(C=C1)C(C)=CC(C)=C2', 'CC1=CC2=C(C(C)=CC=C2)C(C)=C1', 'CC1=CC=CC2=C(C)C=CC(C)=C12', 'CC1=CC2=C(C=C(C)C=C2)C(C)=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1-isopropylnaphthalene', '2-isopropylnaphthalene', '2-methyl-1-ethylnaphthalene', '1-methyl-2-ethylnaphthalene', '2-methyl-3-ethylnaphthalene', '1-methyl-3-ethylnaphthalene', '1-methyl-4-ethylnaphthalene', '2-methyl-4-ethylnaphthalene', '1-methyl-5-ethylnaphthalene', '2-methyl-5-ethylnaphthalene', '2-methyl-6-ethylnaphthalene', '1-methyl-6-ethylnaphthalene', '1-methyl-7-ethylnaphthalene', '2-methyl-7-ethylnaphthalene', '2-methyl-8-ethylnaphthalene', '1-methyl-8-ethylnaphthalene', '2-propylnaphthalene', '1-propylnaphthalene', '1,2,6-trimethylnaphthalene', '1,4,6-trimethylnaphthalene', '1,2,3-trimethylnaphthalene', '2,3,6-trimethylnaphthalene', '1,2,4-trimethylnaphthalene', '1,2,5-trimethylnaphthalene', '1,6,7-trimethylnaphthalene', '1,2,7-trimethylnaphthalene', '1,2,8-trimethylnaphthalene', '1,3,5-trimethylnaphthalene', '1,3,6-trimethylnaphthalene', '1,3,8-trimethylnaphthalene', '1,4,5-trimethylnaphthalene', '1,3,7-trimethylnaphthalene']

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
filename = "C13.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C14 alkylnaphtalenes ++++++

# list of SMILES
smiles_list = ['CCCCC1=C2C=CC=CC2=CC=C1', 'CCCCC1=CC2=C(C=CC=C2)C=C1', 'CCC1=CC2=C(C=CC=C2CC)C=C1', 'CCC1=C(CC)C=C2C=CC=CC2=C1', 'CCC1=C2C=CC=CC2=C(CC)C=C1', 'CCC1=C(CC)C2=C(C=CC=C2)C=C1', 'CCC1=CC2=C(C=C1)C(CC)=CC=C2', 'CCC1=C2C=CC(C)=CC2=C(C)C=C1', 'CCC1=C2C=CC=CC2=C(C)C=C1C', 'CCC1=C(C)C=C2C=CC=C(C)C2=C1', 'CCC1=C2C=CC=CC2=C(C)C(C)=C1', 'CCC1=C2C=C(C)C=C(C)C2=CC=C1', 'CCC1=C2C(C)=CC=C(C)C2=CC=C1', 'CCC1=CC2=C(C=C1)C(C)=CC(C)=C2', 'CCC1=CC2=C(C=C1)C(C)=CC=C2C', 'CCC1=CC2=C(C=C1)C=CC(C)=C2C', 'CC(C)CC1=CC2=C(C=CC=C2)C=C1', 'CC(C)CC1=C2C=CC=CC2=CC=C1', 'CC1=CC=C2C=CC=CC2=C1C(C)C', 'CCCC1=C2C=CC=CC2=CC=C1C', 'CC1=C(C(C)C)C=C2C=CC=CC2=C1', 'CCCC1=C(C)C=C2C=CC=CC2=C1', 'CCCC1=C(C)C2=C(C=CC=C2)C=C1', 'CC1=CC(C(C)C)=CC2=C1C=CC=C2', 'CC1=C2C=CC=CC2=C(C(C)C)C=C1', 'CC1=CC=CC2=C1C=CC(C(C)C)=C2', 'CC1=CC=CC2=C1C=C(C(C)C)C=C2', 'CC1=CC2=C(C(C)C)C=CC=C2C=C1', 'CCC(C)C1=CC2=C(C=CC=C2)C=C1', 'CCC(C)C1=C2C=CC=CC2=CC=C1', 'CC(C)(C)C1=CC2=C(C=CC=C2)C=C1', 'CC(C)(C)C1=C2C=CC=CC2=CC=C1', 'CC1=CC2=C(C(C)=CC=C2C)C(C)=C1', 'CC1=CC2=C(C=C1)C(C)=C(C)C=C2C', 'CC1=CC2=CC(C)=C(C)C(C)=C2C=C1', 'CC1=CC=C(C)C2=C(C)C=CC(C)=C12', 'CC1=C(C)C2=C(C=CC=C2)C(C)=C1C', 'CC1=CC=C(C)C2=C1C=CC(C)=C2C', 'CC1=CC2=C(C=C1)C(C)=CC(C)=C2C', 'CC1=C(C)C2=C(C=C1)C(C)=C(C)C=C2', 'CC1=CC2=C(C)C=CC(C)=C2C=C1C', 'CC1=CC2=CC=C(C)C(C)=C2C(C)=C1', 'CC1=CC(C)=C(C)C2=C1C=CC=C2C', 'CC1=CC2=C(C=C(C)C(C)=C2)C(C)=C1', 'CC1=CC2=CC(C)=C(C)C=C2C=C1C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1-butylnaphthalene', '2-butylnaphthalene', '1,7-diethylnaphthalene', '2,3-diethylnaphthalene', '1,4-diethylnaphthalene', '1,2-diethylnaphthalene', '1,6-diethylnaphthalene', '4,6-dimethyl-1-ethylnaphthalene', '2,4-dimethyl-1-ethylnaphthalene', '2,5-dimethyl-3-ethylnaphthalene', '1,2-dimethyl-4-ethylnaphthalene', '1,3-dimethyl-5-ethylnaphthalene', '1,4-dimethyl-5-ethylnaphthalene', '1,3-dimethyl-6-ethylnaphthalene', '1,4-dimethyl-6-ethylnaphthalene', '1,2-dimethyl-7-ethylnaphthalene', '2-isobutylnaphthalene', '1-isobutylnaphthalene', '2-methyl-1-isopropylnaphthalene', '2-methyl-1-propylnaphthalene', '3-methyl-2-isopropylnaphthalene', '3-methyl-2-propylnaphthalene', '1-methyl-2-propylnaphthalene', '1-methyl-3-isopropylnaphthalene', '1-methyl-4-isopropylnaphthalene', '1-methyl-6-isopropylnaphthalene', '1-methyl-7-isopropylnaphthalene', '2-methyl-8-isopropylnaphthalene', '2-sec-butylnaphthalene', '1-sec-butylnaphthalene', '2-tert-butylnaphthalene', '1-tert-butylnaphthalene', '1,3,5,8-tetramethylnaphthalene', '1,2,4,6-tetramethylnaphthalene', '1,2,3,6-tetramethylnaphthalene', '1,4,5,8-tetramethylnaphthalene', '1,2,3,4-tetramethylnaphthalene', '1,2,5,8-tetramethylnaphthalene', '1,2,4,7-tetramethylnaphthalene', '1,2,5,6-tetramethylnaphthalene', '1,4,6,7-tetramethylnaphthalene', '1,2,6,8-tetramethylnaphthalene', '1,2,4,8-tetramethylnaphthalene', '1,3,6,7-tetramethylnaphthalene', '2,3,6,7-tetramethylnaphthalene']

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
filename = "C14.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C15 alkylnaphtalenes ++++++

# list of SMILES
smiles_list = ['CCCCCC1=CC2=C(C=CC=C2)C=C1', 'CCCCCC1=C2C=CC=CC2=CC=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2-pentylnaphthalene', '1-pentylnaphthalene']

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
filename = "C15.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)