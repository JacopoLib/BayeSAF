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

# ++++++ C7 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCC(C)(C)CC', 'CCCC(C)(C)C', 'CC(C)CC(C)C', 'CCC(C)C(C)C', 'CCC(CC)CC', 'CCCC(C)CC', 'CCCCC(C)C', 'CC(C)C(C)(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['3,3-dimethylpentane', '2,2-dimethylpentane', '2,4-dimethylpentane', '2,3-dimethylpentane', '3-ethylpentane', '3-methylhexane', '2-methylhexane', '2,2,3-trimethylbutane']

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
filename = "iC7.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C8 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCC(C)C(C)CC', 'CCCC(C)(C)CC', 'CC(C)CCC(C)C', 'CCCCC(C)(C)C', 'CCCC(C)C(C)C', 'CCC(C)CC(C)C', 'CCC(CC)C(C)C', 'CCCC(CC)CC', 'CCC(C)(CC)CC', 'CCCCCC(C)C', 'CCCC(C)CCC', 'CCCCC(C)CC', 'CC(C)(C)C(C)(C)C', 'CCC(C)C(C)(C)C', 'CCC(C)(C)C(C)C', 'CC(C)C(C)C(C)C', 'CC(C)CC(C)(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['3,4-dimethylhexane', '3,3-dimethylhexane', '2,5-dimethylhexane', '2,2-dimethylhexane', '2,3-dimethylhexane', '2,4-dimethylhexane', '3-ethyl-2-methylpentane', '3-ethylhexane', '3-methyl-3-ethylpentane', '2-methylheptane', '4-methylheptane', '3-methylheptane', '2,2,3,3-tetramethylbutane', '2,2,3-trimethylpentane', '2,3,3-trimethylpentane', '2,3,4-trimethylpentane', '2,2,4-trimethylpentane']

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
filename = "iC8.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C9 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCC(CC)(CC)CC', 'CCC(CC)C(C)(C)C', 'CCC(C(C)C)C(C)C', 'CCC(C)(CC)C(C)C', 'CCC(C)CCC(C)C', 'CC(C)CCCC(C)C', 'CCCC(C)C(C)CC', 'CCC(C)CC(C)CC', 'CCCCC(C)C(C)C', 'CCCC(C)CC(C)C', 'CCCC(C)(C)CCC', 'CCCCCC(C)(C)C', 'CCCC(CC)C(C)C', 'CCC(CC)CC(C)C', 'CCCC(C)(CC)CC', 'CCC(C)C(CC)CC', 'CCCC(CC)CCC', 'CCCCC(CC)CC', 'CCCCCCC(C)C', 'CCCCCC(C)CC', 'CCCCC(C)CCC', 'CCC(C)(C)C(C)(C)C', 'CC(C)C(C)C(C)(C)C', 'CC(C)C(C)(C)C(C)C', 'CC(C)(C)CC(C)(C)C', 'CC(C)CCC(C)(C)C', 'CCCC(C)C(C)(C)C', 'CC(C)CC(C)C(C)C', 'CCC(C)C(C)(C)CC', 'CCC(C)C(C)C(C)C', 'CCC(C)(C)CC(C)C', 'CCCC(C)(C)C(C)C', 'CCC(C)CC(C)(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['3,3-diethylpentane', '2,2-dimethyl-3-ethylpentane', '2,4-dimethyl-3-ethylpentane', '2,3-dimethyl-3-ethylpentane', '2,5-dimethylheptane', '2,6-dimethylheptane', '3,4-dimethylheptane', '3,5-dimethylheptane', '2,3-dimethylheptane', '2,4-dimethylheptane', '4,4-dimethylheptane', '2,2-dimethylheptane', '3-ethyl-2-methylhexane', '4-ethyl-2-methylhexane', '3-ethyl-3-methylhexane', '3-ethyl-4-methylhexane', '4-ethylheptane', '3-ethylheptane', '2-methyloctane', '3-methyloctane', '4-methyloctane', '2,2,3,3-tetramethylpentane', '2,2,3,4-tetramethylpentane', '2,3,3,4-tetramethylpentane', '2,2,4,4-tetramethylpentane', '2,2,5-trimethylhexane', '2,2,3-trimethylhexane', '2,3,5-trimethylhexane', '3,3,4-trimethylhexane', '2,3,4-trimethylhexane', '2,4,4-trimethylhexane', '2,3,3-trimethylhexane', '2,2,4-trimethylhexane']

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
filename = "iC9.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C10 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCC(C)(C)CC', 'CCCC(C)C(C)CCC', 'CCC(C)CCC(C)CC', 'CCCC(C)CC(C)CC', 'CCCCC(C)(C)CCC', 'CCCCC(C)C(C)CC', 'CC(C)CCCCC(C)C', 'CCC(C)CCCC(C)C', 'CCCC(C)CCC(C)C', 'CCCCCC(C)C(C)C', 'CCCCCCC(C)(C)C', 'CCCCC(C)CC(C)C', 'CCC(C)(CC)C(C)(C)C', 'CCC(C(C)C)C(C)(C)C', 'CCC(CC)CC(C)(C)C', 'CCCC(CC)C(C)(C)C', 'CCC(C)(C(C)C)C(C)C', 'CCC(CC)C(C)C(C)C', 'CCCC(C)(CC)C(C)C', 'CCC(C)(CC)CC(C)C', 'CCC(C)C(CC)C(C)C', 'CCC(CC(C)C)C(C)C', 'CCCC(CC)CC(C)C', 'CCC(CC)CCC(C)C', 'CCCCC(CC)C(C)C', 'CCC(CC)C(C)(C)CC', 'CCC(C)C(C)(CC)CC', 'CCCC(CC)C(C)CC', 'CCCCC(C)(CC)CC', 'CCCC(C)C(CC)CC', 'CCCC(C)(CC)CCC', 'CCC(C)CC(CC)CC', 'CCCCC(CC)CCC', 'CCCCCC(CC)CC', 'CCCC(C(C)C)C(C)C', 'CCCC(CCC)C(C)C', 'CCCCCCC(C)CC', 'CCCCCCCC(C)C', 'CCCCCC(C)CCC', 'CCCCC(C)CCCC', 'CC(C(C)(C)C)C(C)(C)C', 'CC(C)C(C)(C)C(C)(C)C', 'CCCC(CCC)CCC', 'CC(C)C(C)CC(C)(C)C', 'CCC(C)(C)CC(C)(C)C', 'CCCC(C)(C)C(C)(C)C', 'CCC(C)C(C)C(C)(C)C', 'CCC(C)(C)C(C)(C)CC', 'CCC(C)(C)C(C)C(C)C', 'CC(C)CC(C)(C)C(C)C', 'CCC(C)C(C)(C)C(C)C', 'CC(C)C(C)C(C)C(C)C', 'CC(C)(C)CCC(C)(C)C', 'CC(C)CC(C)C(C)(C)C', 'CCC(C)CC(C)C(C)C', 'CCC(C)C(C)CC(C)C', 'CCC(C)C(C)C(C)CC', 'CCCC(C)(C)C(C)CC', 'CCC(C)CC(C)(C)CC', 'CCCC(C)C(C)(C)CC', 'CCCCC(C)(C)C(C)C', 'CC(C)CC(C)CC(C)C', 'CCCC(C)(C)CC(C)C', 'CC(C)CCC(C)C(C)C', 'CCCC(C)C(C)C(C)C', 'CC(C)CCCC(C)(C)C', 'CCC(C)CCC(C)(C)C', 'CCCCC(C)C(C)(C)C', 'CCC(C)(C)CCC(C)C', 'CCCC(C)CC(C)(C)C', 'CCC(CC)(CC)C(C)C', 'CCCC(CC)(CC)CC', 'CCC(CC)C(CC)CC', 'CC(C)C(C(C)C)C(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['3,3-dimethyloctane', '4,5-dimethyloctane', '3,6-dimethyloctane', '3,5-dimethyloctane', '4,4-dimethyloctane', '3,4-dimethyloctane', '2,7-dimethyloctane', '2,6-dimethyloctane', '2,5-dimethyloctane', '2,3-dimethyloctane', '2,2-dimethyloctane', '2,4-dimethyloctane', '3-ethyl-2,2,3-trimethylpentane', '3-ethyl-2,2,4-trimethylpentane', '4-ethyl-2,2-dimethylhexane', '3-ethyl-2,2-dimethylhexane', '3-ethyl-2,3,4-trimethylpentane', '4-ethyl-2,3-dimethylhexane', '3-ethyl-2,3-dimethylhexane', '4-ethyl-2,4-dimethylhexane', '3-ethyl-2,4-dimethylhexane', '3-ethyl-2,5-dimethylhexane', '4-ethyl-2-methylheptane', '5-ethyl-2-methylheptane', '3-ethyl-2-methylheptane', '4-ethyl-3,3-dimethylhexane', '3-ethyl-3,4-dimethylhexane', '4-ethyl-3-methylheptane', '3-ethyl-3-methylheptane', '3-ethyl-4-methylheptane', '4-ethyl-4-methylheptane', '3-ethyl-5-methylheptane', '4-ethyloctane', '3-ethyloctane', '3-isopropyl-2-methylhexane', '4-isopropylheptane', '3-methylnonane', '2-methylnonane', '4-methylnonane', '5-methylnonane', '2,2,3,4,4-pentamethylpentane', '2,2,3,3,4-pentamethylpentane', '4-propylheptane', '2,2,4,5-tetramethylhexane', '2,2,4,4-tetramethylhexane', '2,2,3,3-tetramethylhexane', '2,2,3,4-tetramethylhexane', '3,3,4,4-tetramethylhexane', '2,3,4,4-tetramethylhexane', '2,3,3,5-tetramethylhexane', '2,3,3,4-tetramethylhexane', '2,3,4,5-tetramethylhexane', '2,2,5,5-tetramethylhexane', '2,2,3,5-tetramethylhexane', '2,3,5-trimethylheptane', '2,4,5-trimethylheptane', '3,4,5-trimethylheptane', '3,4,4-trimethylheptane', '3,3,5-trimethylheptane', '3,3,4-trimethylheptane', '2,3,3-trimethylheptane', '2,4,6-trimethylheptane', '2,4,4-trimethylheptane', '2,3,6-trimethylheptane', '2,3,4-trimethylheptane', '2,2,6-trimethylheptane', '2,2,5-trimethylheptane', '2,2,3-trimethylheptane', '2,5,5-trimethylheptane', '2,2,4-trimethylheptane', '3,3-diethyl-2-methylpentane', '3,3-diethylhexane', '3,4-diethylhexane', '2,4-dimethyl-3-isopropylpentane']

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
filename = "iC10.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C11 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCC(CC)C(CC)CC', 'CCCCC(CC)(CC)CC', 'CCCC(CC)(CC)CCC', 'CCC(CC)CC(CC)CC', 'CCC(CC)(C(C)C)C(C)C', 'CCC(CC)(CC)C(C)(C)C', 'CCCC(C)(C)C(CC)CC', 'CCCC(C)C(CC)C(C)C', 'CCCC(C)C(C)(CC)CC', 'CCC(C)CC(C)(CC)CC', 'CCC(CCC(C)C)C(C)C', 'CCC(C)CC(CC)C(C)C', 'CCCCC(C)(CC)C(C)C', 'CCCCC(CC)C(C)(C)C', 'CCC(C)C(C(C)C)C(C)C', 'CCCC(C(C)C)C(C)(C)C', 'CCCC(C)(C(C)C)C(C)C', 'CC(C)CC(C(C)C)C(C)C', 'CCC(C)C(CC)C(C)CC', 'CCCC(CC)C(C)C(C)C', 'CCCC(C)(CC)CC(C)C', 'CCC(C)C(CC)CC(C)C', 'CCC(CC(C)C)CC(C)C', 'CCCC(CC)C(C)(C)CC', 'CCCC(C)(CC)C(C)CC', 'CCCC(CC)CC(C)(C)C', 'CCC(CC)CCC(C)(C)C', 'CCC(CC)CC(C)C(C)C', 'CCC(CC)C(C)CC(C)C', 'CCC(C)(CC)CCC(C)C', 'CCC(CC)CC(C)(C)CC', 'CCC(C)C(C)C(CC)CC', 'CCCCCC(C)(C)CCC', 'CCCCC(C)CCC(C)C', 'CCCCCC(C)C(C)CC', 'CCCCC(C)CC(C)CC', 'CCCC(C)CCC(C)CC', 'CCCC(C)CCCC(C)C', 'CCCCCCC(C)(C)CC', 'CCC(C)CCCC(C)CC', 'CCCCCC(C)CC(C)C', 'CCCCCCC(C)C(C)C', 'CCCCCCCC(C)(C)C', 'CCCCC(C)C(C)CCC', 'CCCC(C)CC(C)CCC', 'CCCCC(C)(C)CCCC', 'CCC(C)CCCCC(C)C', 'CC(C)CCCCCC(C)C', 'CCCCC(CC)CCCC', 'CCCCCC(CC)CCC', 'CCCCCCC(CC)CC', 'CC(C)(C)C(C)(C)C(C)(C)C', 'CCCCC(CCC)C(C)C', 'CCCC(CC)(CC)C(C)C', 'CCC(CC)C(CC)C(C)C', 'CCC(CC)C(C)(CC)CC', 'CCCCCC(CC)C(C)C', 'CCCCC(C)C(CC)CC', 'CCCCCC(C)(CC)CC', 'CCCCC(C(C)C)C(C)C', 'CCC(C)C(CC)(CC)CC', 'CCC(CC)(CC)CC(C)C', 'CCCCC(CC)C(C)CC', 'CCCCC(CC)CC(C)C', 'CCCCC(C)(CC)CCC', 'CCCC(C)(CCC)C(C)C', 'CCCC(CC(C)C)C(C)C', 'CCCC(C(C)C)C(C)CC', 'CCCC(CCC)C(C)CC', 'CCCC(CCC)CC(C)C', 'CCCC(C)(CCC)CCC', 'CCCC(CC)CCC(C)C', 'CCCC(C)C(CC)CCC', 'CCCC(CC)CC(C)CC', 'CCCC(C)CC(CC)CC', 'CCC(C)CCC(CC)CC', 'CCC(CC)CCCC(C)C', 'CCCCCCC(C)CCC', 'CCCCCCCCC(C)C', 'CCCCCCCC(C)CC', 'CCCCCC(C)CCCC', 'CC(C)C(C)C(C)(C)C(C)C', 'CC(C)C(C)C(C)C(C)(C)C', 'CC(C)CC(C)(C)C(C)(C)C', 'CCC(C)C(C)(C)C(C)(C)C', 'CC(CC(C)(C)C)C(C)(C)C', 'CC(C)C(C)(C)CC(C)(C)C', 'CCC(C)(C)C(C)C(C)(C)C', 'CCC(C)(C)C(C)(C)C(C)C', 'CCCCC(CCC)CCC', 'CCCC(CCC)C(C)(C)C', 'CCC(C(C)(C)C)C(C)(C)C', 'CCC(C)(C(C)C)C(C)(C)C', 'CC(C)(C)CCCC(C)(C)C', 'CCC(C)CC(C)(C)C(C)C', 'CCCC(C)C(C)(C)C(C)C', 'CCCC(C)C(C)C(C)(C)C', 'CC(C)C(C)CCC(C)(C)C', 'CCCC(C)(C)C(C)C(C)C', 'CCC(C)(C)CCC(C)(C)C', 'CCC(C)(C)C(C)CC(C)C', 'CC(C)CC(C)CC(C)(C)C', 'CC(C)CCC(C)(C)C(C)C', 'CCCC(C)(C)C(C)(C)CC', 'CCCCC(C)(C)C(C)(C)C', 'CCC(C)C(C)(C)C(C)CC', 'CC(C)CCC(C)C(C)(C)C', 'CCC(C)(C)CC(C)(C)CC', 'CCC(C)C(C)C(C)(C)CC', 'CCC(C)C(C)(C)CC(C)C', 'CCC(C)CC(C)C(C)(C)C', 'CCC(C)C(C)C(C)C(C)C', 'CC(C)CC(C)(C)CC(C)C', 'CC(C)C(C)CC(C)C(C)C', 'CCC(C)(C)CC(C)C(C)C', 'CCCC(C)(C)CC(C)(C)C', 'CC(C)CC(C)C(C)C(C)C', 'CCC(C)C(C)CC(C)(C)C', 'CCC(C)C(C)(CC)C(C)C', 'CCC(CC(C)C)C(C)(C)C', 'CCC(C)(CC(C)C)C(C)C', 'CCC(C(C)C)C(C)(C)CC', 'CCC(C)C(CC)C(C)(C)C', 'CCCC(C)(CC)C(C)(C)C', 'CC(C)C(C(C)C)C(C)(C)C', 'CC(C)C(C)(C(C)C)C(C)C', 'CCC(CC)C(C)(C)C(C)C', 'CCC(CC(C)(C)C)C(C)C', 'CCC(C)(CC)C(C)C(C)C', 'CCC(C)(CC)CC(C)(C)C', 'CCC(C)(C)C(C)(CC)CC', 'CCC(C(C)C)C(C)C(C)C', 'CCC(CC)C(C)C(C)(C)C', 'CCC(C)CCCC(C)(C)C', 'CCCCC(C)(C)CC(C)C', 'CC(C)CCCC(C)C(C)C', 'CCC(C)CCC(C)C(C)C', 'CCCC(C)C(C)CC(C)C', 'CCCCCC(C)C(C)(C)C', 'CCCC(C)CC(C)C(C)C', 'CCC(C)CC(C)CC(C)C', 'CCCC(C)CCC(C)(C)C', 'CCCCC(C)(C)C(C)CC', 'CCCCC(C)C(C)C(C)C', 'CCCCCC(C)(C)C(C)C', 'CCCCC(C)CC(C)(C)C', 'CCCCC(C)C(C)(C)CC', 'CC(C)CCCCC(C)(C)C', 'CCCC(C)C(C)(C)CCC', 'CCCC(C)(C)CC(C)CC', 'CCC(C)CC(C)C(C)CC', 'CCCC(C)CC(C)(C)CC', 'CCC(C)CCC(C)(C)CC', 'CC(C)CCC(C)CC(C)C', 'CCC(C)(C)CCCC(C)C', 'CCC(C)C(C)CCC(C)C', 'CCCC(C)(C)CCC(C)C', 'CCCC(C)C(C)C(C)CC']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['3,4-diethylheptane', '3,3-diethylheptane', '4,4-diethylheptane', '3,5-diethylheptane', '2,4-dimethyl-3,3-diethylpentane', '2,2-dimethyl-3,3-diethylpentane', '4,4-dimethyl-3-ethylheptane', '2,4-dimethyl-3-ethylheptane', '3,4-dimethyl-3-ethylheptane', '3,5-dimethyl-3-ethylheptane', '2,6-dimethyl-3-ethylheptane', '2,5-dimethyl-3-ethylheptane', '2,3-dimethyl-3-ethylheptane', '2,2-dimethyl-3-ethylheptane', '2,4-dimethyl-3-isopropylhexane', '2,2-dimethyl-3-isopropylhexane', '2,3-dimethyl-3-isopropylhexane', '2,5-dimethyl-3-isopropylhexane', '3,5-dimethyl-4-ethylheptane', '2,3-dimethyl-4-ethylheptane', '2,4-dimethyl-4-ethylheptane', '2,5-dimethyl-4-ethylheptane', '2,6-dimethyl-4-ethylheptane', '3,3-dimethyl-4-ethylheptane', '3,4-dimethyl-4-ethylheptane', '2,2-dimethyl-4-ethylheptane', '2,2-dimethyl-5-ethylheptane', '2,3-dimethyl-5-ethylheptane', '2,4-dimethyl-5-ethylheptane', '2,5-dimethyl-5-ethylheptane', '3,3-dimethyl-5-ethylheptane', '3,4-dimethyl-5-ethylheptane', '4,4-dimethylnonane', '2,5-dimethylnonane', '3,4-dimethylnonane', '3,5-dimethylnonane', '3,6-dimethylnonane', '2,6-dimethylnonane', '3,3-dimethylnonane', '3,7-dimethylnonane', '2,4-dimethylnonane', '2,3-dimethylnonane', '2,2-dimethylnonane', '4,5-dimethylnonane', '4,6-dimethylnonane', '5,5-dimethylnonane', '2,7-dimethylnonane', '2,8-dimethylnonane', '5-ethylnonane', '4-ethylnonane', '3-ethylnonane', '2,2,3,3,4,4-hexamethylpentane', '4-isopropyloctane', '2-methyl-3,3-diethylhexane', '2-methyl-3,4-diethylhexane', '3-methyl-3,4-diethylhexane', '2-methyl-3-ethyloctane', '4-methyl-3-ethyloctane', '3-methyl-3-ethyloctane', '2-methyl-3-isopropylheptane', '3-methyl-4,4-diethylhexane', '2-methyl-4,4-diethylhexane', '3-methyl-4-ethyloctane', '2-methyl-4-ethyloctane', '4-methyl-4-ethyloctane', '4-methyl-4-isopropylheptane', '2-methyl-4-isopropylheptane', '3-methyl-4-isopropylheptane', '3-methyl-4-propylheptane', '2-methyl-4-propylheptane', '4-methyl-4-propylheptane', '2-methyl-5-ethyloctane', '4-methyl-5-ethyloctane', '3-methyl-5-ethyloctane', '4-methyl-6-ethyloctane', '3-methyl-6-ethyloctane', '2-methyl-6-ethyloctane', '4-methyldecane', '2-methyldecane', '3-methyldecane', '5-methyldecane', '2,3,3,4,5-pentamethylhexane', '2,2,3,4,5-pentamethylhexane', '2,2,3,3,5-pentamethylhexane', '2,2,3,3,4-pentamethylhexane', '2,2,3,5,5-pentamethylhexane', '2,2,4,4,5-pentamethylhexane', '2,2,3,4,4-pentamethylhexane', '2,3,3,4,4-pentamethylhexane', '4-propyloctane', '4-tert-butylheptane', '2,2,4,4-tetramethyl-3-ethylpentane', '2,2,3,4-tetramethyl-3-ethylpentane', '2,2,6,6-tetramethylheptane', '2,3,3,5-tetramethylheptane', '2,3,3,4-tetramethylheptane', '2,2,3,4-tetramethylheptane', '2,2,5,6-tetramethylheptane', '2,3,4,4-tetramethylheptane', '2,2,5,5-tetramethylheptane', '2,4,5,5-tetramethylheptane', '2,2,4,6-tetramethylheptane', '2,3,3,6-tetramethylheptane', '3,3,4,4-tetramethylheptane', '2,2,3,3-tetramethylheptane', '3,4,4,5-tetramethylheptane', '2,2,3,6-tetramethylheptane', '3,3,5,5-tetramethylheptane', '3,3,4,5-tetramethylheptane', '2,4,4,5-tetramethylheptane', '2,2,3,5-tetramethylheptane', '2,3,4,5-tetramethylheptane', '2,4,4,6-tetramethylheptane', '2,3,5,6-tetramethylheptane', '2,3,5,5-tetramethylheptane', '2,2,4,4-tetramethylheptane', '2,3,4,6-tetramethylheptane', '2,2,4,5-tetramethylheptane', '2,3,4-trimethyl-3-ethylhexane', '2,2,5-trimethyl-3-ethylhexane', '2,3,5-trimethyl-3-ethylhexane', '2,4,4-trimethyl-3-ethylhexane', '2,2,4-trimethyl-3-ethylhexane', '2,2,3-trimethyl-3-ethylhexane', '2,2,4-trimethyl-3-isopropylpentane', '2,3,4-trimethyl-3-isopropylpentane', '2,3,3-trimethyl-4-ethylhexane', '2,2,5-trimethyl-4-ethylhexane', '2,3,4-trimethyl-4-ethylhexane', '2,2,4-trimethyl-4-ethylhexane', '3,3,4-trimethyl-4-ethylhexane', '2,3,5-trimethyl-4-ethylhexane', '2,2,3-trimethyl-4-ethylhexane', '2,2,6-trimethyloctane', '2,4,4-trimethyloctane', '2,3,7-trimethyloctane', '2,3,6-trimethyloctane', '2,4,5-trimethyloctane', '2,2,3-trimethyloctane', '2,3,5-trimethyloctane', '2,4,6-trimethyloctane', '2,2,5-trimethyloctane', '3,4,4-trimethyloctane', '2,3,4-trimethyloctane', '2,3,3-trimethyloctane', '2,2,4-trimethyloctane', '3,3,4-trimethyloctane', '2,2,7-trimethyloctane', '4,4,5-trimethyloctane', '3,5,5-trimethyloctane', '3,4,6-trimethyloctane', '3,3,5-trimethyloctane', '3,3,6-trimethyloctane', '2,4,7-trimethyloctane', '2,6,6-trimethyloctane', '2,5,6-trimethyloctane', '2,5,5-trimethyloctane', '3,4,5-trimethyloctane']

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
filename = "iC11.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C12 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCC(CC)(CC)CCC', 'CCCCCC(CC)(CC)CC', 'CCC(CC)CCC(CC)CC', 'CCCC(CC)CC(CC)CC', 'CCCC(CC)C(CC)CCC', 'CCCCC(CC)C(CC)CC', 'CCCC(CC)(CC)C(C)(C)C', 'CCC(CC)(CC(C)C)C(C)C', 'CCC(C)C(CC)(CC)C(C)C', 'CCC(CC)C(CC)C(C)(C)C', 'CCC(CC)C(C)(CC)C(C)C', 'CCC(C(C)C)C(C)(CC)CC', 'CCC(C(C)C)C(CC)C(C)C', 'CCC(C)(CC)C(C)(CC)CC', 'CCC(C(C)C)(C(C)C)C(C)C', 'CCCCC(C)C(C)(CC)CC', 'CCCCCC(CC)C(C)(C)C', 'CCCCCC(C)(CC)C(C)C', 'CCC(C)CCC(C)(CC)CC', 'CCCCC(C)(C)C(CC)CC', 'CCCC(C)C(C)C(CC)CC', 'CCCC(C)CC(CC)C(C)C', 'CCCC(C)CC(C)(CC)CC', 'CCC(C)CCC(CC)C(C)C', 'CCCCC(C)C(CC)C(C)C', 'CCC(CCCC(C)C)C(C)C', 'CCCC(C)C(C(C)C)C(C)C', 'CCCCC(C(C)C)C(C)(C)C', 'CCC(C)CC(C(C)C)C(C)C', 'CC(C)CCC(C(C)C)C(C)C', 'CCCCC(C)(C(C)C)C(C)C', 'CCCC(C(C)(C)C)C(C)(C)C', 'CCC(CC)(CC)CC(C)(C)C', 'CCC(CC)(CC)C(C)C(C)C', 'CCC(C)(C)C(CC)(CC)CC', 'CCCC(C)C(CC)CC(C)C', 'CCCCC(C)(CC)CC(C)C', 'CCC(C)CC(CC)CC(C)C', 'CCCCC(CC)C(C)C(C)C', 'CCC(CCC(C)C)CC(C)C', 'CCC(C)CC(CC)C(C)CC', 'CCCC(C)C(C)(CC)CCC', 'CCCCC(C)(CC)C(C)CC', 'CCCC(C)C(CC)C(C)CC', 'CCCCC(CC)CC(C)(C)C', 'CCCCC(CC)C(C)(C)CC', 'CCC(C)C(C(C)C)C(C)CC', 'CCCC(C(C)C)C(C)C(C)C', 'CCCC(C)(CC(C)C)C(C)C', 'CCC(C)C(CC(C)C)C(C)C', 'CCCC(C(C)C)C(C)(C)CC', 'CC(C)CC(CC(C)C)C(C)C', 'CCCC(C)(C(C)C)C(C)CC', 'CCCC(CC(C)(C)C)C(C)C', 'CCCC(CCC)C(C)(C)CC', 'CCCC(CC(C)C)C(C)CC', 'CCCC(CCC)C(C)C(C)C', 'CCCC(CCC)CC(C)(C)C', 'CCCC(CC(C)C)CC(C)C', 'CCCC(C)(CCC)C(C)CC', 'CCCC(C(C)CC)C(C)CC', 'CCCC(C)(CCC)CC(C)C', 'CCCC(C)(CC)CC(C)CC', 'CCCC(CC)C(C)(C)CCC', 'CCCC(CC)CC(C)(C)CC', 'CCCC(C)(CC)CCC(C)C', 'CCCC(CC)CCC(C)(C)C', 'CCCC(CC)CC(C)C(C)C', 'CCC(C)C(CC)CCC(C)C', 'CCCC(CC)C(C)CC(C)C', 'CCCC(CC)C(C)C(C)CC', 'CCC(CC)CCC(C)(C)CC', 'CCC(CC)CC(C)CC(C)C', 'CCC(CC)CC(C)C(C)CC', 'CCCC(C)(C)CC(CC)CC', 'CCC(C)CC(C)C(CC)CC', 'CCC(CC)CCCC(C)(C)C', 'CCC(CC)C(C)CCC(C)C', 'CCC(C)(CC)CCCC(C)C', 'CCC(CC)CCC(C)C(C)C', 'CCCCC(C)CC(C)CCC', 'CCCC(C)CCC(C)CCC', 'CCCCCC(C)(C)CCCC', 'CC(C)CCCCCCC(C)C', 'CCCCCCCCC(C)(C)C', 'CCCCCCCC(C)C(C)C', 'CCCCC(C)C(C)CCCC', 'CCCCCC(C)CCC(C)C', 'CCCCC(C)CCCC(C)C', 'CCCCCCC(C)CC(C)C', 'CCC(C)CCCCCC(C)C', 'CCCCCC(C)C(C)CCC', 'CCCCCCCC(C)(C)CC', 'CCCCCCC(C)C(C)CC', 'CCCCCC(C)CC(C)CC', 'CCCCC(C)CCC(C)CC', 'CCCC(C)CCCC(C)CC', 'CCC(C)CCCCC(C)CC', 'CCCCCCC(C)(C)CCC', 'CCCC(C)CCCCC(C)C', 'CCCC(CC)(CCC)C(C)C', 'CCCC(C(C)C)C(CC)CC', 'CCCC(CC)(CCC)CCC', 'CCCC(CCC)C(CC)CC', 'CCCCCCC(CC)CCC', 'CCCCCCCC(CC)CC', 'CCCCCC(CC)CCCC', 'CC(C)C(C)(C)C(C)(C)C(C)C', 'CCC(C)(C)C(C)(C)C(C)(C)C', 'CC(C)C(C)C(C)(C)C(C)(C)C', 'CC(C)C(C)(C)C(C)C(C)(C)C', 'CC(C)(C)CC(C)(C)C(C)(C)C', 'CC(C(C)C(C)(C)C)C(C)(C)C', 'CCCCC(CCCC)C(C)C', 'CCCCCC(CCC)C(C)C', 'CCCCC(CC)(CC)C(C)C', 'CCCC(C)C(CC)(CC)CC', 'CCCC(CC)C(C)(CC)CC', 'CCCC(CC)C(CC)C(C)C', 'CCCC(C)(CC)C(CC)CC', 'CCC(CC)CC(CC)C(C)C', 'CCC(CC)CC(C)(CC)CC', 'CCC(CC)C(C)C(CC)CC', 'CCCC(CC)(C(C)C)C(C)C', 'CCCCC(C)CC(CC)CC', 'CCCCCCC(C)(CC)CC', 'CCCCCCC(CC)C(C)C', 'CCCCCC(C)C(CC)CC', 'CCCCCC(C(C)C)C(C)C', 'CCCC(CC)(CC)CC(C)C', 'CCCC(CC)(CC)C(C)CC', 'CCC(CC)C(CC)CC(C)C', 'CCC(C)C(CC)C(CC)CC', 'CCC(CC)C(C(C)C)C(C)C', 'CCCCCC(CC)C(C)CC', 'CCCCC(C)C(CC)CCC', 'CCCCCC(C)(CC)CCC', 'CCCCCC(CC)CC(C)C', 'CCCCC(CC(C)C)C(C)C', 'CCCCC(C)(CCC)C(C)C', 'CCCCC(C(C)C)C(C)CC', 'CCCCC(C)(CCC)CCC', 'CCCCC(CCC)CC(C)C', 'CCCCC(CCC)C(C)CC', 'CCCC(C)(CCC)C(C)(C)C', 'CCCC(CC(C)C)C(C)(C)C', 'CCCC(C(C)CC)C(C)(C)C', 'CCC(CC)(CC)CCC(C)C', 'CCC(C)CC(CC)(CC)CC', 'CCCCC(CC)CC(C)CC', 'CCCCC(CC)CCC(C)C', 'CCCCC(CC)C(C)CCC', 'CCCCC(C)(CC)CCCC', 'CCCC(C)C(CCC)C(C)C', 'CCCC(CC(C)CC)C(C)C', 'CCCC(CCC(C)C)C(C)C', 'CCCC(C)C(CCC)CCC', 'CCCC(CCC)CC(C)CC', 'CCCC(CCC)CCC(C)C', 'CCCC(CC)CCCC(C)C', 'CCCC(C)CC(CC)CCC', 'CCCC(CC)CCC(C)CC', 'CCC(CC)CCCCC(C)C', 'CCC(C)CCCC(CC)CC', 'CCCC(C)CCC(CC)CC', 'CCCCCCC(C)CCCC', 'CCCCCC(C)CCCCC', 'CCCCCCCC(C)CCC', 'CCCCCCCCC(C)CC', 'CCCCCCCCCC(C)C', 'CCC(C)(C(C)(C)C)C(C)(C)C', 'CC(CC(C)(C)C)CC(C)(C)C', 'CC(C)CC(C)(C)CC(C)(C)C', 'CCC(C)C(C)(C)CC(C)(C)C', 'CCCC(C)(C)C(C)C(C)(C)C', 'CC(C)CCC(C)(C)C(C)(C)C', 'CCCC(C)C(C)(C)C(C)(C)C', 'CC(C)C(C)C(C)CC(C)(C)C', 'CCC(C)CC(C)(C)C(C)(C)C', 'CCC(C)(C)C(C)CC(C)(C)C', 'CC(C)C(C)C(C)C(C)C(C)C', 'CC(C)C(C)CC(C)C(C)(C)C', 'CC(C)CC(C)C(C)(C)C(C)C', 'CC(CCC(C)(C)C)C(C)(C)C', 'CCC(C)C(C)C(C)(C)C(C)C', 'CCCC(C)(C)C(C)(C)C(C)C', 'CC(C)C(C)(C)CCC(C)(C)C', 'CC(C)C(C)CC(C)(C)C(C)C', 'CCC(C)C(C)(C)C(C)C(C)C', 'CCC(C)(C)CC(C)C(C)(C)C', 'CC(C)CC(C)(C)C(C)C(C)C', 'CCC(C)(C)C(C)C(C)C(C)C', 'CCC(C)(C)CC(C)(C)C(C)C', 'CCC(C)(C)C(C)(C)CC(C)C', 'CCC(C)C(C)(C)C(C)(C)CC', 'CCC(C)(C)C(C)C(C)(C)CC', 'CCC(C)C(C)C(C)C(C)(C)C', 'CC(C)CC(C)C(C)C(C)(C)C', 'CCCCC(CCC)CCCC', 'CCCCCC(CCC)CCC', 'CCCCC(CCC)C(C)(C)C', 'CCC(C(C)(C)C)C(C)(C)CC', 'CCC(C)(CC(C)C)C(C)(C)C', 'CCC(C)(C(C)C)C(C)C(C)C', 'CCC(C)(C)C(C)(CC)C(C)C', 'CCC(CC(C)(C)C)C(C)(C)C', 'CCC(C)C(C)(CC)C(C)(C)C', 'CCC(C(C)C(C)C)C(C)(C)C', 'CC(C)C(C(C)(C)C)C(C)(C)C', 'CC(C)C(C)(C(C)C)C(C)(C)C', 'CCC(C)(CC)C(C)(C)C(C)C', 'CCC(C(C)C)C(C)C(C)(C)C', 'CCC(CC)C(C)(C)C(C)(C)C', 'CCC(C)(CC(C)(C)C)C(C)C', 'CCC(C)(CC)C(C)C(C)(C)C', 'CCC(C(C)C)C(C)(C)C(C)C', 'CCC(C)CC(C)CC(C)(C)C', 'CC(C)CC(C)CC(C)C(C)C', 'CCC(C)CCC(C)C(C)(C)C', 'CCCC(C)CC(C)C(C)(C)C', 'CCCCC(C)C(C)C(C)(C)C', 'CCCCCC(C)(C)C(C)(C)C', 'CC(C)CCCC(C)C(C)(C)C', 'CCCCC(C)(C)CC(C)(C)C', 'CCCC(C)C(C)CC(C)(C)C', 'CCCC(C)C(C)C(C)C(C)C', 'CC(C)C(C)CCCC(C)(C)C', 'CCCC(C)C(C)(C)CC(C)C', 'CC(C)C(C)CCC(C)C(C)C', 'CCC(C)(C)CCC(C)C(C)C', 'CCC(C)C(C)CC(C)C(C)C', 'CCCC(C)(C)CC(C)C(C)C', 'CC(C)CCC(C)(C)CC(C)C', 'CCC(C)CC(C)C(C)C(C)C', 'CCCC(C)(C)C(C)CC(C)C', 'CCCCC(C)(C)C(C)C(C)C', 'CC(C)CCCC(C)(C)C(C)C', 'CCC(C)CCC(C)(C)C(C)C', 'CCCC(C)CC(C)(C)C(C)C', 'CCCCC(C)C(C)(C)C(C)C', 'CC(C)(C)CCCCC(C)(C)C', 'CC(C)CCC(C)C(C)C(C)C', 'CCC(C)CC(C)C(C)(C)CC', 'CCCC(C)(C)C(C)(C)CCC', 'CCC(C)C(C)C(C)C(C)CC', 'CCCC(C)(C)C(C)C(C)CC', 'CCC(C)CC(C)(C)C(C)CC', 'CCCC(C)C(C)(C)C(C)CC', 'CCC(C)(C)CCC(C)(C)CC', 'CCC(C)CC(C)(C)CC(C)C', 'CCCC(C)(C)CC(C)(C)CC', 'CC(C)CC(C)C(C)CC(C)C', 'CCCC(C)C(C)C(C)(C)CC', 'CCCCC(C)(C)C(C)(C)CC', 'CCC(C)(C)C(C)CCC(C)C', 'CCC(C)C(C)(C)CCC(C)C', 'CCC(C)(C)CC(C)CC(C)C', 'CCC(C)C(C)C(C)CC(C)C', 'CCC(C)C(C)CC(C)(C)CC', 'CCC(C)(C)CCCC(C)(C)C', 'CC(C)CC(C)CCC(C)(C)C', 'CCC(C)C(C)CCC(C)(C)C', 'CC(C)CCC(C)CC(C)(C)C', 'CCCC(C)(C)CCC(C)(C)C', 'CCC(CC)C(CC)(CC)CC', 'CCC(CC)(C(C)C)C(C)(C)C', 'CCC(C)C(C)C(C)(CC)CC', 'CCC(C)CC(CC)C(C)(C)C', 'CCC(CCC(C)C)C(C)(C)C', 'CCC(C)(CCC(C)C)C(C)C', 'CCCC(C)(C)C(CC)C(C)C', 'CCC(CC(C)(C)CC)C(C)C', 'CCCC(C)C(CC)C(C)(C)C', 'CCCCC(C)(CC)C(C)(C)C', 'CCC(C)C(C)C(CC)C(C)C', 'CCCC(C)C(C)(CC)C(C)C', 'CCCC(C)(C)C(C)(CC)CC', 'CCC(C)CC(C)(CC)C(C)C', 'CCC(C(C)C)C(C)CC(C)C', 'CC(C)CC(C)(C(C)C)C(C)C', 'CCC(C)C(C)(C(C)C)C(C)C', 'CC(C)CC(C(C)C)C(C)(C)C', 'CCC(C)C(C(C)C)C(C)(C)C', 'CCCC(C)(C(C)C)C(C)(C)C', 'CCC(C)(C)C(C(C)C)C(C)C', 'CCCC(C)(CC)C(C)C(C)C', 'CCC(CC(C)C)C(C)C(C)C', 'CCC(C)C(CC)C(C)C(C)C', 'CCCC(CC)C(C)C(C)(C)C', 'CCCC(CC)C(C)(C)C(C)C', 'CCC(CC(C)C)C(C)(C)CC', 'CCC(C)C(C)(CC)CC(C)C', 'CCC(C)C(C)(CC)C(C)CC', 'CCC(C)C(CC)CC(C)(C)C', 'CCC(CC(C)C)CC(C)(C)C', 'CCC(C)(CC(C)C)CC(C)C', 'CCCC(C)(CC)C(C)(C)CC', 'CCC(C)C(CC)C(C)(C)CC', 'CCCC(C)(CC)CC(C)(C)C', 'CC(C)C(CC(C)(C)C)C(C)C', 'CC(C)C(C)C(C(C)C)C(C)C', 'CCC(CC)CC(C)C(C)(C)C', 'CCC(CCC(C)(C)C)C(C)C', 'CCC(CC)CC(C)(C)C(C)C', 'CCC(CC)C(C)CC(C)(C)C', 'CCC(C)(CC)CCC(C)(C)C', 'CCC(C)C(C)(C)C(CC)CC', 'CCC(CC)C(C)C(C)C(C)C', 'CCC(CC(C)C(C)C)C(C)C', 'CCC(CC)C(C)(C)CC(C)C', 'CCC(C)(C)CC(C)(CC)CC', 'CCC(C)(CC)C(C)CC(C)C', 'CCC(C)(CC)CC(C)C(C)C', 'CCC(CC)C(C)C(C)(C)CC', 'CCC(C)(C)CCCCC(C)C', 'CCCC(C)CC(C)C(C)CC', 'CCCC(C)CC(C)CC(C)C', 'CCC(C)CC(C)CC(C)CC', 'CCC(C)CCC(C)CC(C)C', 'CC(C)CCCC(C)CC(C)C', 'CCCC(C)C(C)CC(C)CC', 'CCCCC(C)(C)CCC(C)C', 'CCCCC(C)(C)CC(C)CC', 'CCCCC(C)C(C)CC(C)C', 'CCC(C)CCC(C)C(C)CC', 'CCCC(C)CCC(C)(C)CC', 'CCCCC(C)C(C)C(C)CC', 'CCCCCC(C)(C)C(C)CC', 'CCC(C)CCCC(C)(C)CC', 'CCCCCC(C)C(C)(C)CC', 'CCC(C)C(C)CCCC(C)C', 'CCCC(C)(C)CCCC(C)C', 'CC(C)CCC(C)CCC(C)C', 'CCC(C)CC(C)CCC(C)C', 'CCCC(C)C(C)CCC(C)C', 'CCCCCC(C)C(C)C(C)C', 'CCCCC(C)CC(C)(C)CC', 'CCCCC(C)CCC(C)(C)C', 'CCCC(C)CCCC(C)(C)C', 'CCC(C)CCCCC(C)(C)C', 'CCCCCC(C)CC(C)(C)C', 'CC(C)CCCCCC(C)(C)C', 'CCCC(C)(C)CCC(C)CC', 'CCCCCCC(C)(C)C(C)C', 'CCCCCCC(C)C(C)(C)C', 'CCCCC(C)CC(C)C(C)C', 'CCCC(C)C(C)C(C)CCC', 'CCCCC(C)(C)C(C)CCC', 'CCCC(C)CCC(C)C(C)C', 'CCC(C)CCCC(C)C(C)C', 'CCCC(C)CC(C)(C)CCC', 'CC(C)CCCCC(C)C(C)C', 'CCCCC(C)C(C)(C)CCC', 'CCCCCC(C)(C)CC(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['4,4-diethyloctane', '3,3-diethyloctane', '3,6-diethyloctane', '3,5-diethyloctane', '4,5-diethyloctane', '3,4-diethyloctane', '2,2-dimethyl,3,3-diethylhexane', '2,5-dimethyl-3,3-diethylhexane', '2,4-dimethyl-3,3-diethylhexane', '2,2-dimethyl-3,4-diethylhexane', '2,3-dimethyl-3,4-diethylhexane', '2,4-dimethyl-3,4-diethylhexane', '2,5-dimethyl-3,4-diethylhexane', '3,4-dimethyl-3,4-diethylhexane', '2,4-dimethyl-3-ethyl-3-isopropylpentane', '3,4-dimethyl-3-ethyloctane', '2,2-dimethyl-3-ethyloctane', '2,3-dimethyl-3-ethyloctane', '3,6-dimethyl-3-ethyloctane', '4,4-dimethyl-3-ethyloctane', '4,5-dimethyl-3-ethyloctane', '2,5-dimethyl-3-ethyloctane', '3,5-dimethyl-3-ethyloctane', '2,6-dimethyl-3-ethyloctane', '2,4-dimethyl-3-ethyloctane', '2,7-dimethyl-3-ethyloctane', '2,4-dimethyl-3-isopropylheptane', '2,2-dimethyl-3-isopropylheptane', '2,5-dimethyl-3-isopropylheptane', '2,6-dimethyl-3-isopropylheptane', '2,3-dimethyl-3-isopropylheptane', '2,2-dimethyl-3-tert-butylhexane', '2,2-dimethyl-4,4-diethylhexane', '2,3-dimethyl-4,4-diethylhexane', '3,3-dimethyl-4,4-diethylhexane', '2,5-dimethyl-4-ethyloctane', '2,4-dimethyl-4-ethyloctane', '2,6-dimethyl-4-ethyloctane', '2,3-dimethyl-4-ethyloctane', '2,7-dimethyl-4-ethyloctane', '3,6-dimethyl-4-ethyloctane', '4,5-dimethyl-4-ethyloctane', '3,4-dimethyl-4-ethyloctane', '3,5-dimethyl-4-ethyloctane', '2,2-dimethyl-4-ethyloctane', '3,3-dimethyl-4-ethyloctane', '3,5-dimethyl-4-isopropylheptane', '2,3-dimethyl-4-isopropylheptane', '2,4-dimethyl-4-isopropylheptane', '2,5-dimethyl-4-isopropylheptane', '3,3-dimethyl-4-isopropylheptane', '2,6-dimethyl-4-isopropylheptane', '3,4-dimethyl-4-isopropylheptane', '2,2-dimethyl-4-isopropylheptane', '3,3-dimethyl-4-propylheptane', '2,5-dimethyl-4-propylheptane', '2,3-dimethyl-4-propylheptane', '2,2-dimethyl-4-propylheptane', '2,6-dimethyl-4-propylheptane', '3,4-dimethyl-4-propylheptane', '3,5-dimethyl-4-propylheptane', '2,4-dimethyl-4-propylheptane', '3,5-dimethyl-5-ethyloctane', '4,4-dimethyl-5-ethyloctane', '3,3-dimethyl-5-ethyloctane', '2,5-dimethyl-5-ethyloctane', '2,2-dimethyl-5-ethyloctane', '2,3-dimethyl-5-ethyloctane', '2,6-dimethyl-5-ethyloctane', '2,4-dimethyl-5-ethyloctane', '3,4-dimethyl-5-ethyloctane', '3,3-dimethyl-6-ethyloctane', '2,4-dimethyl-6-ethyloctane', '3,4-dimethyl-6-ethyloctane', '4,4-dimethyl-6-ethyloctane', '3,5-dimethyl-6-ethyloctane', '2,2-dimethyl-6-ethyloctane', '2,5-dimethyl-6-ethyloctane', '2,6-dimethyl-6-ethyloctane', '2,3-dimethyl-6-ethyloctane', '4,6-dimethyldecane', '4,7-dimethyldecane', '5,5-dimethyldecane', '2,9-dimethyldecane', '2,2-dimethyldecane', '2,3-dimethyldecane', '5,6-dimethyldecane', '2,5-dimethyldecane', '2,6-dimethyldecane', '2,4-dimethyldecane', '2,8-dimethyldecane', '4,5-dimethyldecane', '3,3-dimethyldecane', '3,4-dimethyldecane', '3,5-dimethyldecane', '3,6-dimethyldecane', '3,7-dimethyldecane', '3,8-dimethyldecane', '4,4-dimethyldecane', '2,7-dimethyldecane', '4-ethyl-4-isopropylheptane', '3-ethyl-4-isopropylheptane', '4-ethyl-4-propylheptane', '3-ethyl-4-propylheptane', '4-ethyldecane', '3-ethyldecane', '5-ethyldecane', '2,3,3,4,4,5-hexamethylhexane', '2,2,3,3,4,4-hexamethylhexane', '2,2,3,3,4,5-hexamethylhexane', '2,2,3,4,4,5-hexamethylhexane', '2,2,3,3,5,5-hexamethylhexane', '2,2,3,4,5,5-hexamethylhexane', '5-isopropylnonane', '4-isopropylnonane', '2-methyl-3,3-diethylheptane', '4-methyl-3,3-diethylheptane', '3-methyl-3,4-diethylheptane', '2-methyl-3,4-diethylheptane', '4-methyl-3,4-diethylheptane', '2-methyl-3,5-diethylheptane', '3-methyl-3,5-diethylheptane', '4-methyl-3,5-diethylheptane', '2-methyl-3-ethyl-3-isopropylhexane', '5-methyl-3-ethylnonane', '3-methyl-3-ethylnonane', '2-methyl-3-ethylnonane', '4-methyl-3-ethylnonane', '2-methyl-3-isopropyloctane', '2-methyl-4,4-diethylheptane', '3-methyl-4,4-diethylheptane', '2-methyl-4,5-diethylheptane', '3-methyl-4,5-diethylheptane', '2-methyl-4-ethyl-3-isopropylhexane', '3-methyl-4-ethylnonane', '5-methyl-4-ethylnonane', '4-methyl-4-ethylnonane', '2-methyl-4-ethylnonane', '2-methyl-4-isopropyloctane', '4-methyl-4-isopropyloctane', '3-methyl-4-isopropyloctane', '4-methyl-4-propyloctane', '2-methyl-4-propyloctane', '3-methyl-4-propyloctane', '4-methyl-4-tert-butylheptane', '2-methyl-4-tert-butylheptane', '3-methyl-4-tert-butylheptane', '2-methyl-5,5-diethylheptane', '3-methyl-5,5-diethylheptane', '3-methyl-5-ethylnonane', '2-methyl-5-ethylnonane', '4-methyl-5-ethylnonane', '5-methyl-5-ethylnonane', '4-methyl-5-isopropyloctane', '3-methyl-5-isopropyloctane', '2-methyl-5-isopropyloctane', '4-methyl-5-propyloctane', '3-methyl-5-propyloctane', '2-methyl-5-propyloctane', '2-methyl-6-ethylnonane', '4-methyl-6-ethylnonane', '3-methyl-6-ethylnonane', '2-methyl-7-ethylnonane', '3-methyl-7-ethylnonane', '4-methyl-7-ethylnonane', '5-methylundecane', '6-methylundecane', '4-methylundecane', '3-methylundecane', '2-methylundecane', '2,2,3,4,4-pentamethyl-3-ethylpentane', '2,2,4,6,6-pentamethylheptane', '2,2,4,4,6-pentamethylheptane', '2,2,4,4,5-pentamethylheptane', '2,2,3,4,4-pentamethylheptane', '2,2,3,3,6-pentamethylheptane', '2,2,3,3,4-pentamethylheptane', '2,2,4,5,6-pentamethylheptane', '2,2,3,3,5-pentamethylheptane', '2,2,4,5,5-pentamethylheptane', '2,3,4,5,6-pentamethylheptane', '2,2,3,5,6-pentamethylheptane', '2,3,3,4,6-pentamethylheptane', '2,2,3,6,6-pentamethylheptane', '2,3,3,4,5-pentamethylheptane', '2,3,3,4,4-pentamethylheptane', '2,2,5,5,6-pentamethylheptane', '2,3,3,5,6-pentamethylheptane', '2,3,4,4,5-pentamethylheptane', '2,2,3,5,5-pentamethylheptane', '2,3,4,4,6-pentamethylheptane', '2,3,4,5,5-pentamethylheptane', '2,3,3,5,5-pentamethylheptane', '2,4,4,5,5-pentamethylheptane', '3,3,4,4,5-pentamethylheptane', '3,3,4,5,5-pentamethylheptane', '2,2,3,4,5-pentamethylheptane', '2,2,3,4,6-pentamethylheptane', '5-propylnonane', '4-propylnonane', '4-tert-butyloctane', '2,2,4,4-tetramethyl-3-ethylhexane', '2,2,3,5-tetramethyl-3-ethylhexane', '2,3,4,5-tetramethyl-3-ethylhexane', '2,3,4,4-tetramethyl-3-ethylhexane', '2,2,5,5-tetramethyl-3-ethylhexane', '2,2,3,4-tetramethyl-3-ethylhexane', '2,2,4,5-tetramethyl-3-ethylhexane', '2,2,4,4-tetramethyl-3-isopropylpentane', '2,2,3,4-tetramethyl-3-isopropylpentane', '2,3,3,4-tetramethyl-4-ethylhexane', '2,2,3,5-tetramethyl-4-ethylhexane', '2,2,3,3-tetramethyl-4-ethylhexane', '2,2,4,5-tetramethyl-4-ethylhexane', '2,2,3,4-tetramethyl-4-ethylhexane', '2,3,3,5-tetramethyl-4-ethylhexane', '2,2,4,6-tetramethyloctane', '2,3,5,7-tetramethyloctane', '2,2,3,6-tetramethyloctane', '2,2,3,5-tetramethyloctane', '2,2,3,4-tetramethyloctane', '2,2,3,3-tetramethyloctane', '2,2,3,7-tetramethyloctane', '2,2,4,4-tetramethyloctane', '2,2,4,5-tetramethyloctane', '2,3,4,5-tetramethyloctane', '2,2,6,7-tetramethyloctane', '2,4,4,5-tetramethyloctane', '2,3,6,7-tetramethyloctane', '2,3,6,6-tetramethyloctane', '2,3,5,6-tetramethyloctane', '2,3,5,5-tetramethyloctane', '2,4,4,7-tetramethyloctane', '2,3,4,6-tetramethyloctane', '2,4,5,5-tetramethyloctane', '2,3,4,4-tetramethyloctane', '2,3,3,7-tetramethyloctane', '2,3,3,6-tetramethyloctane', '2,3,3,5-tetramethyloctane', '2,3,3,4-tetramethyloctane', '2,2,7,7-tetramethyloctane', '2,3,4,7-tetramethyloctane', '3,3,4,6-tetramethyloctane', '4,4,5,5-tetramethyloctane', '3,4,5,6-tetramethyloctane', '3,4,5,5-tetramethyloctane', '3,4,4,6-tetramethyloctane', '3,4,4,5-tetramethyloctane', '3,3,6,6-tetramethyloctane', '2,4,4,6-tetramethyloctane', '3,3,5,5-tetramethyloctane', '2,4,5,7-tetramethyloctane', '3,3,4,5-tetramethyloctane', '3,3,4,4-tetramethyloctane', '2,5,6,6-tetramethyloctane', '2,5,5,6-tetramethyloctane', '2,4,6,6-tetramethyloctane', '2,4,5,6-tetramethyloctane', '3,3,5,6-tetramethyloctane', '2,2,6,6-tetramethyloctane', '2,2,5,7-tetramethyloctane', '2,2,5,6-tetramethyloctane', '2,2,4,7-tetramethyloctane', '2,2,5,5-tetramethyloctane', '3,3,4-triethylhexane', '2,2,4-trimethyl-3,3-diethylpentane', '3,4,5-trimethyl-3-ethylheptane', '2,2,5-trimethyl-3-ethylheptane', '2,2,6-trimethyl-3-ethylheptane', '2,3,6-trimethyl-3-ethylheptane', '2,4,4-trimethyl-3-ethylheptane', '2,5,5-trimethyl-3-ethylheptane', '2,2,4-trimethyl-3-ethylheptane', '2,2,3-trimethyl-3-ethylheptane', '2,4,5-trimethyl-3-ethylheptane', '2,3,4-trimethyl-3-ethylheptane', '3,4,4-trimethyl-3-ethylheptane', '2,3,5-trimethyl-3-ethylheptane', '2,4,6-trimethyl-3-ethylheptane', '2,3,5-trimethyl-3-isopropylhexane', '2,3,4-trimethyl-3-isopropylhexane', '2,2,5-trimethyl-3-isopropylhexane', '2,2,4-trimethyl-3-isopropylhexane', '2,2,3-trimethyl-3-isopropylhexane', '2,4,4-trimethyl-3-isopropylhexane', '2,3,4-trimethyl-4-ethylheptane', '2,3,6-trimethyl-4-ethylheptane', '2,3,5-trimethyl-4-ethylheptane', '2,2,3-trimethyl-4-ethylheptane', '2,3,3-trimethyl-4-ethylheptane', '2,5,5-trimethyl-4-ethylheptane', '2,4,5-trimethyl-4-ethylheptane', '3,4,5-trimethyl-4-ethylheptane', '2,2,5-trimethyl-4-ethylheptane', '2,2,6-trimethyl-4-ethylheptane', '2,4,6-trimethyl-4-ethylheptane', '3,3,4-trimethyl-4-ethylheptane', '3,3,5-trimethyl-4-ethylheptane', '2,2,4-trimethyl-4-ethylheptane', '2,2,5-trimethyl-4-isopropylhexane', '2,3,5-trimethyl-4-isopropylhexane', '2,2,3-trimethyl-5-ethylheptane', '2,2,6-trimethyl-5-ethylheptane', '2,3,3-trimethyl-5-ethylheptane', '2,2,4-trimethyl-5-ethylheptane', '2,2,5-trimethyl-5-ethylheptane', '3,4,4-trimethyl-5-ethylheptane', '2,3,4-trimethyl-5-ethylheptane', '2,3,6-trimethyl-5-ethylheptane', '2,4,4-trimethyl-5-ethylheptane', '3,3,5-trimethyl-5-ethylheptane', '2,4,5-trimethyl-5-ethylheptane', '2,3,5-trimethyl-5-ethylheptane', '3,3,4-trimethyl-5-ethylheptane', '2,7,7-trimethylnonane', '3,4,6-trimethylnonane', '2,4,6-trimethylnonane', '3,5,7-trimethylnonane', '2,4,7-trimethylnonane', '2,4,8-trimethylnonane', '3,5,6-trimethylnonane', '2,5,5-trimethylnonane', '3,5,5-trimethylnonane', '2,4,5-trimethylnonane', '3,4,7-trimethylnonane', '3,3,6-trimethylnonane', '3,4,5-trimethylnonane', '3,4,4-trimethylnonane', '3,3,7-trimethylnonane', '3,3,4-trimethylnonane', '2,6,7-trimethylnonane', '2,6,6-trimethylnonane', '2,5,8-trimethylnonane', '2,5,7-trimethylnonane', '2,5,6-trimethylnonane', '2,3,4-trimethylnonane', '3,3,5-trimethylnonane', '2,2,5-trimethylnonane', '2,2,6-trimethylnonane', '2,2,7-trimethylnonane', '2,2,4-trimethylnonane', '2,2,8-trimethylnonane', '3,6,6-trimethylnonane', '2,3,3-trimethylnonane', '2,2,3-trimethylnonane', '2,3,5-trimethylnonane', '4,5,6-trimethylnonane', '4,5,5-trimethylnonane', '2,3,6-trimethylnonane', '2,3,7-trimethylnonane', '4,4,6-trimethylnonane', '2,3,8-trimethylnonane', '4,4,5-trimethylnonane', '2,4,4-trimethylnonane']

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
filename = "iC12.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C13 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCC(C)CC(C)C', 'CCCCCCCCCC(C)(C)C', 'CCCCCCCCC(C)C(C)C', 'CCCCCCCCCCC(C)C', 'CCCCCCCCCC(C)CC']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2,4-dimethylundecane', '2,2-dimethylundecane', '2,3-dimethylundecane', '2-methyldodecane', '3-methyldodecane']

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
filename = "iC13.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C14 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCCCCC(C)(C)C', 'CCCCCCCCC(C)CC(C)C', 'CCCCCCCCCC(C)C(C)C', 'CCCCCCCCCCCC(C)C', 'CCCCCCCCCCC(C)CC']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2,2-dimethyldodecane', '2,4-dimethyldodecane', '2,3-dimethyldodecane', '2-methyltridecane', '3-methyltridecane']

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
filename = "iC14.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C15 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCCCCC(C)C(C)C', 'CCCCCCCCCCCC(C)(C)C', 'CCCCCCCCCC(C)CC(C)C', 'CCCCCCCCCCCCC(C)C', 'CCCCCCCCCCCC(C)CC']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2,3-dimethyltridecane', '2,2-dimethyltridecane', '2,4-dimethyltridecane', '2-methyltetradecane', '3-methyltetradecane']

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
filename = "iC15.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C16 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCCCCCCC(C)(C)C', 'CCCCCCCCCCC(C)CC(C)C', 'CCCCCCCCCCCC(C)C(C)C', 'CCCCCCCCCCCCC(C)CC', 'CCCCCCCCCCCCCC(C)C', 'CC(CC(C)(C)C)CC(C)(C)CC(C)(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2,2-dimethyltetradecane', '2,4-dimethyltetradecane', '2,3-dimethyltetradecane', '3-methylpentadecane', '2-methylpentadecane', '2,2,4,4,6,8,8-heptamethylnonane']

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
filename = "iC16.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C17 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCCCCCCC(C)C(C)C', 'CCCCCCCCCCCC(C)CC(C)C', 'CCCCCCCCCCCCCC(C)(C)C', 'CCCCCCCCCCCCCCC(C)C', 'CCCCCCCCCCCCCC(C)CC']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2,3-dimethylpentadecane', '2,4-dimethylpentadecane', '2,2-dimethylpentadecane', '2-methylhexadecane', '3-methylhexadecane']

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
filename = "iC17.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C18 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCCCCCCCCC(C)(C)C', 'CCCCCCCCCCCCCC(C)C(C)C', 'CCCCCCCCCCCCC(C)CC(C)C', 'CCCCCCCCCCCCCCCC(C)C', 'CCCCCCCCCCCCCCC(C)CC']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2,2-dimethylhexadecane', '2,3-dimethylhexadecane', '2,4-dimethylhexadecane', '2-methylheptadecane', '3-methylheptadecane']

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
filename = "iC18.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)

# ++++++ C19 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCCCCCCCCCC(C)(C)C', 'CCCCCCCCCCCCCC(C)CC(C)C', 'CCCCCCCCCCCCCCC(C)C(C)C', 'CCCCCCCCCCCCCCCC(C)CC', 'CCCCCCCCCCCCCCCCC(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['2,2-dimethylheptadecane', '2,4-dimethylheptadecane', '2,3-dimethylheptadecane', '3-methyloctadecane', '2-methyloctadecane']

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
filename = "iC19.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)
    
# ++++++ C20 iso-paraffins ++++++

# list of SMILES
smiles_list = ['CCCCCCCCCCCCCCC(C)CC(C)C', 'CCCCCCCCCCCCCCCC(C)C(C)C', 'CCCCCCCCCCCCCCCCC(C)(C)C', 'CCCCCCCCCCCCCCCCC(C)CC', 'CCCCCCCCCCCCCCCCCC(C)C']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

species = ['2,4-dimethyloctadecane', '2,3-dimethyloctadecane', '2,2-dimethyloctadecane', '3-methylnonadecane', '2-methylnonadecane']

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
filename = "iC20.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_star)