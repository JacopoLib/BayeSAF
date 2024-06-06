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

# ++++++ C6 alkylbenzenes ++++++

# list of SMILES
smiles_list = ['C1=CC=CC=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['benzene']

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
    
# ++++++ C7 alkylbenzenes ++++++

# list of SMILES
smiles_list = ['CC1=CC=CC=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['toluene']

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
    
# ++++++ C8 alkylbenzenes ++++++

# list of SMILES
smiles_list = ['CCC1=CC=CC=C1', 'CC1=CC=CC(C)=C1', 'CC1=C(C)C=CC=C1', 'CC1=CC=C(C)C=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['ethylbenzene', 'm-xylene', 'o-xylene', 'p-xylene']

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
    
# ++++++ C9 alkylbenzenes ++++++

# list of SMILES
smiles_list = ['CC(C)C1=CC=CC=C1', 'CCC1=CC=C(C)C=C1', 'CCC1=CC(C)=CC=C1', 'CCC1=C(C)C=CC=C1', 'CC1=CC(C)=CC(C)=C1', 'CCCC1=CC=CC=C1', 'CC1=CC=CC(C)=C1C', 'CC1=CC(C)=C(C)C=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cumene', 'p-ethyltoluene', 'm-ethyltoluene', 'o-ethyltoluene', 'mesitylene', 'propylbenzene', '1,2,3-trimethylbenzene', '1,2,4-trimethylbenzene']

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

# ++++++ C10 alkylbenzenes ++++++

# list of SMILES
smiles_list = ['CC(C)(C)C1=CC=CC=C1', 'CCC(C)C1=CC=CC=C1', 'CCCCC1=CC=CC=C1', 'CC1=CC=CC(C(C)C)=C1', 'CC1=CC=C(C(C)C)C=C1', 'CC1=C(C(C)C)C=CC=C1', 'CCC1=CC=CC(CC)=C1', 'CCC1=C(CC)C=CC=C1', 'CCC1=CC=C(CC)C=C1', 'CCC1=CC(C)=CC(C)=C1', 'CCC1=C(C)C=C(C)C=C1', 'CCC1=C(C)C=CC=C1C', 'CCC1=CC(C)=C(C)C=C1', 'CCC1=C(C)C(C)=CC=C1', 'CCC1=C(C)C=CC(C)=C1', 'CC(C)CC1=CC=CC=C1', 'CCCC1=C(C)C=CC=C1', 'CCCC1=CC(C)=CC=C1', 'CCCC1=CC=C(C)C=C1', 'CC1=CC(C)=C(C)C(C)=C1', 'CC1=CC(C)=C(C)C=C1C', 'CC1=C(C)C(C)=C(C)C=C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['tert-butylbenzene', 'sec-butylbenzene', 'butylbenzene', 'm-cymene', 'p-cymene', 'o-cymene', 'm-diethylbenzene', 'o-diethylbenzene', 'p-diethylbenzene', '5-ethyl-m-xylene', '4-ethyl-m-xylene', '2-ethyl-m-xylene', '4-ethyl-o-xylene', '3-ethyl-o-xylene', '2-ethyl-p-xylene', 'isobutylbenzene', '1-methyl-2-propylbenzene', '1-methyl-3-propylbenzene', '1-methyl-4-propylbenzene', '1,2,3,5-tetramethylbenzene', '1,2,4,5-tetramethylbenzene', '1,2,3,4-tetramethylbenzene']

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