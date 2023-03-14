# prepare fingerprint XML
import glob

xml_files = glob.glob('fingerprints_xml/*.xml')
xml_files.sort()
xml_files

FP_list = ['AtomPairs2DCount',
 'AtomPairs2D',
 'EState',
 'CDKextended',
 'CDK',
 'CDKgraphonly',
 'KlekotaRothCount',
 'KlekotaRoth',
 'MACCS',
 'PubChem',
 'SubstructureCount',
 'Substructure']

# create a dictionary
fp = dict(zip(FP_list, xml_files))

# prepare data subset as input to PaDEL 
import pandas as pd

df3 = pd.read_csv('acetylcholinesterase_04_bioactivity_data_3clase_pIC50.csv')
df3

selection = ['canonical_smiles', 'molecule_chembl_id']
df3_selection = df3[selection]
df3_selection.to_csv('molecule.smi', sep='\t',index=False, header=False)

# calculate descriptors
# there are 12 fingerprint types in PaDEL. to calculate all 12, make sure to make 
# adjustments to the descriptortypes input augument to any of the ones in the fp dictionary
# variable as shown above, eg, substructurefingerprintcount.xml
fp
from padelpy import padeldescriptor

fingerprint = 'Substructure'

fingerprint_output_file = ''.join([fingerprint, '.csv']) # Substructure.csv
fingerprint_descriptortypes = fp[fingerprint]

padeldescriptor(mol_dir = 'molecule.smi',
                d_file=fingerprint_output_file, #Substructure.csv
                # descriptortypes = 'SubstructureFingerprint.xml',
                descriptortypes = fingerprint_descriptortypes,
                detectaromaticity=True,
                standardizenitro = True,
                standardizetautomers = True,
                threads=2,
                removesalt = True,
                log=True,
                fingerprints=True)

# 测试另一种fingerprint
fingerprint = 'PubChem'

fingerprint_output_file = ''.join([fingerprint, '.csv']) # Substructure.csv
fingerprint_descriptortypes = fp[fingerprint]

padeldescriptor(mol_dir = 'molecule.smi',
                d_file=fingerprint_output_file, #Substructure.csv
                # descriptortypes = 'SubstructureFingerprint.xml',
                descriptortypes = fingerprint_descriptortypes,
                detectaromaticity=True,
                standardizenitro = True,
                standardizetautomers = True,
                threads=2,
                removesalt = True,
                log=True,
                fingerprints=True)

# 利用循环把每一个fingerprint都做了
for i in FP_list:
    fingerprint = i

    fingerprint_output_file = ''.join([fingerprint, '.csv']) # Substructure.csv
    fingerprint_descriptortypes = fp[fingerprint]

    padeldescriptor(mol_dir = 'molecule.smi',
                    d_file=fingerprint_output_file, #Substructure.csv
                    # descriptortypes = 'SubstructureFingerprint.xml',
                    descriptortypes = fingerprint_descriptortypes,
                    detectaromaticity=True,
                    standardizenitro = True,
                    standardizetautomers = True,
                    threads=2,
                    removesalt = True,
                    log=True,
                    fingerprints=True)

# X data matrix
df3_x = pd.read_csv('Substructure.csv')
df3_x
df3_x = df3_x.drop(columns=['Name'])
df3_x

# Y variable
df3_y = df3['pIC50']
df3_y

# combining X and Y variable
dataset3 = pd.concat([df3_x, df3_y], axis=1)
dataset3

dataset3.to_csv('acetylcholinesterase_06_bioactivity_data_3class_pIC50_substructure_fp.csv', index=False)



