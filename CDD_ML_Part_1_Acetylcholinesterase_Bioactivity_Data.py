import pandas as pd
from chembl_webresource_client.new_client import new_client

# target search for acetylcholinesterase
target = new_client.target
target_query = target.search('acetylcholinesterase')
targets = pd.DataFrame.from_dict(target_query)
targets

# select and retrieve bioactivity data for human acetylcholinesterase (first entry)
selected_target = targets.target_chembl_id[0]
selected_target

# retrieve only bioactivity data for human acetylcholinesterase (CHEMBL220) that are reported as pCHEMBL values
activity = new_client.activity
res = activity.filter(target_chembl_id = selected_target).filter(standard_type='IC50')

df = pd.DataFrame.from_dict(res)
df

df.to_csv('acetylcholinesterase_01_bioactivity_data_raw.csv', index=False)

# handling missing data
df2 = df[df.standard_value.notna()]
df2 = df2[df2.canonical_smiles.notna()]
df2

len(df2.canonical_smiles.unique()) #查看唯一分子式的数量

df2_nr = df2.drop_duplicates(['canonical_smiles']) #去掉重复的分子
df2_nr

# data pre-processing of the bioactivity data
# combine the 3 columns (molecule_chembl_id, canonical_smiles, standard_value) and the bioactivity_class into a dataframe
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2_nr[selection]
df3

df3.to_csv('acetylcholinesterase_02_bioactivity_data_preprocessed.csv', index=False)

# labeling compounds as either being active, inactive or intermediate
# the bioactivity data is in the IC50 unit. compounds having values of less than 1000nM will be considered to be active
# while those greater than 10,000nM will be considered to be inactive,
# those in between 1000 and 10000 nM will be referred to as intermediate

df4 = pd.read_csv('acetylcholinesterase_02_bioactivity_data_preprocessed.csv')

bioactivity_threshold = []
for i in df4.standard_value:
    if float(i) >= 10000:
        bioactivity_threshold.append('inactive')
    elif float(i) <= 1000:
        bioactivity_threshold.append('active')
    else:
        bioactivity_threshold.append('intermediate')

bioactivity_class = pd.Series(bioactivity_threshold, name='class')
df5 = pd.concat([df4, bioactivity_class], axis=1)
df5
    
df5.to_csv('acetylcholinesterase_03_bioactivity_data_curated.csv', index=False)

