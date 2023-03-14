import pandas as pd
from chembl_webresource_client.new_client import new_client

#target search for coronavirus
target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame.from_dict(target_query)
targets

# select and retrieve bioactivity data for SARS coronavirus 3C-like proteinase (fifth entry)
selected_target = targets.target_chembl_id[4]
selected_target

# retrieve only bioactivity data for coronavirus 3C-like proteinase (CHEMBL3927) that are reported as IC50 values in nM unit
activity = new_client.activity
res = activity.filter(target_chembl_id = selected_target).filter(standard_type='IC50')
df = pd.DataFrame.from_dict(res)
df.head(3)
df.standard_type.unique() #查看包含类型

# save data
df.to_csv('bioactivity_data.csv', index=False)

#handling missing data
#if any compounds has missing value for the standard_value column then drop it
df2 = df[df.standard_value.notna()]
pd.set_option('display.max_columns', None)
df2

# data pre-processing of the bioactivity data
# labeling compounds as either being active, inactive or intermediate
# the bioactivity data is in the IC50 unit. compounds having values of less than 1000nM will be consider to be active,
# while those greater than 10,000 nM will be considered to be inactive. as those values in between 1000 and 10,000 nM will be referred to as intermediate

bioactivity_class = []
for i in df2.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append('inactive')
    elif float(i) <= 1000:
        bioactivity_class.append('active')
    else:
        bioactivity_class.append('intermediate')

# iterate the molecule_chembl_id to a list
mol_cid = []
for i in df2.molecule_chembl_id:
    mol_cid.append(i)

# iterate canonical_smiles to a list
canonical_smiles = []
for i in df2.canonical_smiles:
    canonical_smiles.append(i)

# iterate standard_value to a list
standard_value = []
for i in df2.standard_value:
    standard_value.append(i)

# combine the 4 lists into a dataframe
data_tuples = list(zip(mol_cid, canonical_smiles, bioactivity_class, standard_value))
df3 = pd.DataFrame(data_tuples, columns=['molecule_chembl_id', 'canonical_smiles', 'bioactivity_class', 'standard_value'])
df3

# alternative method (直接从数据里取出对应的列)
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2[selection]
df3

pd.concat([df3, pd.Series(bioactivity_class)], axis=1)

# save dataframe to csv file
df3.to_csv('bioactivity_preprocessed_data.csv', index=False)
