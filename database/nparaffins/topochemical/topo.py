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

# list of SMILES
smiles_list = ['CCCCCCC', 'CCCCCCCC', 'CCCCCCCCC', 'CCCCCCCCCC', 'CCCCCCCCCCC', 'CCCCCCCCCCCC', 'CCCCCCCCCCCCC', 'CCCCCCCCCCCCCC', 'CCCCCCCCCCCCCCC', 'CCCCCCCCCCCCCCCC', 'CCCCCCCCCCCCCCCCC', 'CCCCCCCCCCCCCCCCCC', 'CCCCCCCCCCCCCCCCCCC', 'CCCCCCCCCCCCCCCCCCCC']

# calculate molecular descriptors
descriptors = from_smiles(smiles_list)

# list of chemical species names
species = ['heptane', 'octane', 'nonane', 'decane', 'undecane', 'dodecane', 'tridecane', 'tetradecane', 'pentadecane', 'hexadecane', 'heptadecane', 'octadecane', 'nonadecane', 'eicosane']

eta_B_star = []

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
				
	eta_B_star.append(eta_B*vol)

if has_duplicates(eta_B_star) == True:
	
	indexes = []
	
	for i, v in enumerate(eta_B_star):
		if eta_B_star.count(v) > 1 and i not in indexes:
			indexes.append(i)
			
	print('The list of duplicate elements is :  ' + str(indexes))
	
# name of csv file
filename = "nC7-nC20.csv"

# writing to csv file
with open(filename, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
 
    # writing the fields
    csvwriter.writerow(species)
 
    # writing the data rows
    csvwriter.writerow(eta_B_star)