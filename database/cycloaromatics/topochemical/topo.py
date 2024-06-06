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

# ++++++ C9 cycloaromatic compounds ++++++

# list of SMILES
smiles_list = ['C1=CC2=C(C=C1)CCC2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['indane']

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
    
# ++++++ C10 cycloaromatic compounds ++++++

# list of SMILES
smiles_list = ['C1=CC2=C(C=C1)CCCC2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1,2,3,4-tetrahydronaphthalene']

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
    
# ++++++ C11 cycloaromatic compounds ++++++

# list of SMILES
smiles_list = ['CC1CCCC2=C1C=CC=C2', 'CC1CCC2=C(C=CC=C2)C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1-methyl-[1,2,3,4-tetrahydronaphthalene]', '2-methyl-[1,2,3,4-tetrahydronaphthalene]']

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

# ++++++ C12 cycloaromatic compounds ++++++

# list of SMILES
smiles_list = ['CCC1CCCC2=C1C=CC=C2', 'CCC1CCC2=C(C=CC=C2)C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1-ethyl-[1,2,3,4-tetrahydronaphthalene]', '2-ethyl-[1,2,3,4-tetrahydronaphthalene]']

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
    
# ++++++ C13 cycloaromatic compounds ++++++

# list of SMILES
smiles_list = ['CCCC1CCC2=C(C=CC=C2)C1', 'CCCC1CCCC2=C1C=CC=C2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2-propyl-[1,2,3,4-tetrahydronaphthalene]', '1-propyl-[1,2,3,4-tetrahydronaphthalene]']

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
    
# ++++++ C14 cycloaromatic compounds ++++++

# list of SMILES
smiles_list = ['CCCCC1CCC2=C(C=CC=C2)C1', 'CCCCC1CCCC2=C1C=CC=C2']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2-butyl-[1,2,3,4-tetrahydronaphthalene]', '1-butyl-[1,2,3,4-tetrahydronaphthalene]']

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
    
# ++++++ C15 cycloaromatic compounds ++++++

# list of SMILES
smiles_list = ['CCCCCC1CCCC2=C1C=CC=C2', 'CCCCCC1CCC2=C(C=CC=C2)C1']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['1-pentyl-[1,2,3,4-tetrahydronaphthalene]', '2-pentyl-[1,2,3,4-tetrahydronaphthalene]']

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