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

# ++++++ C10 dicycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CC[C@@H]2CCCC[C@@H]2C1', 'C1CC[C@H]2CCCC[C@@H]2C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['cis-decahydronaphthalene', 'trans-decahydronaphthalene']

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
    
# ++++++ C11 dicycloparaffins ++++++

# list of SMILES
smiles_list = ['CC1CCC[C@@H]2CCCC[C@H]12', 'CC1CC[C@@H]2CCCC[C@@H]2C1', 'C1CCC[C@]2(CCCC[C@@H]12)C', 'CC1CCC[C@H]2CCCC[C@H]12', 'C[C@]12CCCC[C@H]2CCCC1', 'C[C@H]1CC[C@H]2CCCC[C@@H]2C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1-methyl-cis-decahydronaphthalene', '2-methyl-cis-decahydronaphthalene', '4a-methyl-cis-decahydronaphthalene', '1-methyl-trans-decahydronaphthalene', '4a-methyl-trans-decahydronaphthalene', '2-methyl-trans-decahydronaphthalene']

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
    
# ++++++ C12 dicycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCC(CC1)C2CCCCC2', 'CC1CCCC2(C1CCCC2)C', 'CCC1CCC[C@H]2CCCC[C@@H]12', 'CCC12CCCCC1CCCC2', 'CCC1CC[C@@H]2CCCC[C@@H]2C1', 'CCC1CCCC2CCCCC12', 'CCC1CCC2CCCCC2C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['bicyclohexyl', '1,4a-dimethyl-trans-decahydronaphthalene', '1-ethyl,cis-decahydronaphthalene', '4a-ethyl-cis-decahydronaphthalene', '2-ethyl-cis-decahydronaphthalene', '1-ethyl-trans-decahydronaphthalene', '2-ethyl-trans-decahydronaphthalene']

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
    
# ++++++ C13 dicycloparaffins ++++++

# list of SMILES
smiles_list = ['C1CCC(CC1)CC2CCCCC2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['dicyclohexylmethane']

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
    
# ++++++ C14 dicycloparaffins ++++++

# list of SMILES
smiles_list = ['CC(C1CCCCC1)C2CCCCC2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1,1-dicyclohexylethane']

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
    
# ++++++ C15 dicycloparaffins ++++++

# list of SMILES
smiles_list = ['CCC(C1CCCCC1)C2CCCCC2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1,1-dicyclohexylpropane']

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