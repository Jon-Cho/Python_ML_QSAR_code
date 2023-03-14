import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

# load bioactivity data
df = pd.read_csv('bioactivity_preprocessed_data.csv')

# calculate lipinski descriptors
'''The Lipinski's Rule stated the following:

Molecular weight < 500 Dalton
Octanol-water partition coefficient (LogP) < 5
Hydrogen bond donors < 5
Hydrogen bond acceptors < 10'''

def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
    
    baseData = np.arange(1,1)
    i = 0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])
        
        if (i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors

df_lipinski = lipinski(df.canonical_smiles)

# combine the 2 dataframe
df_combined = pd.concat([df, df_lipinski], axis=1)
df_combined

# convert IC50 to pIC50
# to allow IC50 data to be more uniformly distributed, we will convert IC50 to the negative logarithmic scale which is essentially -log10(IC50)
def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # convert nM to M
        pIC50.append(-np.log10(molar))
    
    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)

    return x

# point to note: values greater than 100,000,000 will be fixed at 100,000,000 otherwise the negative logarithmic value will become negative
df_combined.standard_value.describe()

def norm_value(input):
    norm = []
    
    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000
        norm.append(i) 

    input['standard_value_norm'] = norm   
    x = input.drop('standard_value', 1)

    return x

# first, apply the norm_value() function so that the values in the standard_value is normalized
df_norm = norm_value(df_combined)
df_norm
df_norm.standard_value_norm.describe()        

df_final = pIC50(df_norm)
df_final
df_final.pIC50.describe()

# remove the 'intermediat' bioactivity class
df_2class = df_final[df_final.bioactivity_class != 'intermediate']
df_2class

# exploaratory data analysis (chemical space analysis) via Lipinski descriptors
# Frequency plot of the 2 bioactivity classes
plt.figure(figsize=(5.5, 5.5))
sns.countplot(x = 'bioactivity_class', data = df_2class, edgecolor = 'k')

plt.xlabel('Bioactivity class', fontsize=14, fontweight = 'bold')
plt.ylabel('Frequency', fontsize = 14, fontweight = 'bold')

plt.savefig('plot_bioactivity_class.pdf')

# scatter plot of MW versus LogP
plt.figure(figsize=(5.5, 5.5))

sns.scatterplot(x='MW', y='LogP', data=df_2class, hue = 'bioactivity_class', size='pIC50', edgecolor = 'k', alpha=0.7)

plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.tight_layout()

plt.savefig('plot_MW_vs_LogP.pdf')

# box plots of pIC50 values
plt.figure(figsize=(5.5,5.5))

sns.boxplot(x='bioactivity_class', y='pIC50', data=df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')

plt.savefig('plot_ic50.pdf')

# statistical analysis (Mann-Whitney U Test) 它假设两个样本分别来自除了总体均值以外完全相同的两个总体，目的是检验这两个总体的均值是否有显著的差别
def mannwhitney (descriptor, verbose=False):
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu

    # seed the random number generator
    seed (1)

    # activates and inactives
    selection = [descriptor, 'bioactivity_class']
    df = df_2class[selection]
    active = df[df.bioactivity_class == 'active']
    active = active[descriptor]

    selection = [descriptor, 'bioactivity_class']
    df = df_2class[selection]
    inactive = df[df.bioactivity_class == 'inactive']
    inactive = inactive[descriptor]

    # compare samples
    stat, p = mannwhitneyu(active, inactive)

    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'
    
    results = pd.DataFrame({'Descriptor': descriptor,
                            'Statistics':stat,
                            'p':p,
                            'alpha':alpha,
                            'Interpretation':interpretation}, index=[0])
    filename = 'mannwhitneyu_'+descriptor+'.csv'
    results.to_csv(filename)

    return results

mannwhitney('pIC50')

# plot MW vs bioactivity
plt.figure(figsize=(5.5,5.5))

sns.boxplot(x='bioactivity_class', y='MW', data=df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')

plt.savefig('plot_MW.pdf')

mannwhitney('MW')

# plot LogP vs bioactivity
plt.figure(figsize=(5.5,5.5))

sns.boxplot(x='bioactivity_class', y='LogP', data=df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')

plt.savefig('plot_LogP.pdf')

mannwhitney('LogP')

# plot NumHDonors vs bioactivity
plt.figure(figsize=(5.5,5.5))

sns.boxplot(x='bioactivity_class', y='NumHDonors', data=df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')

plt.savefig('plot_NumHDonors.pdf')

mannwhitney('NumHDonors')

# plot NumHAcceptors vs bioactivity
plt.figure(figsize=(5.5,5.5))

sns.boxplot(x='bioactivity_class', y='NumHAcceptors', data=df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')

plt.savefig('plot_NumHAcceptors.pdf')

mannwhitney('NumHAcceptors')

'''Interpretation of Statistical Results
Box Plots
pIC50 values
Taking a look at pIC50 values, the actives and inactives displayed 
*statistically significant difference*, which is to be expected since 
threshold values (IC50 < 1,000 nM = Actives while IC50 > 10,000 nM = Inactives, 
corresponding to pIC50 > 6 = Actives and pIC50 < 5 = Inactives) were used to 
define actives and inactives.

Lipinski's descriptors
Of the 4 Lipinski's descriptors (MW, LogP, NumHDonors and NumHAcceptors), 
only LogP exhibited *no difference* between the actives and inactives 
while the other 3 descriptors (MW, NumHDonors and NumHAcceptors) shows 
*statistically significant difference* between actives and inactives.'''

