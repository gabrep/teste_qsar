import pandas as pd
from padelpy import padeldescriptor

#Criar descritores para as moléculas com atividade sob MMP1
#Os descritores podem ser criados com Mordred ou com PaDEL
#Nesse caso, criarei com PaDEL para aprender a utilizar a ferramenta e o pacote padeldescriptor

df = pd.read_csv('../data/MMP1_bioactivity_3class_pIC50.csv')
df.head()
selection = ['canonical_smiles', 'molecule_chembl_id']
df_selection = df[selection]

#Criar arquivo .smi contendo o SMILES e o chembl_id de cada molécula, sem header, para input no PaDEL descriptors
df_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)

padeldescriptor(mol_dir='molecule.smi', d_file='../descriptors.csv', removesalt=True, fingerprints=True)

df3_x = pd.read_csv('../descriptors.csv')

df3_x = df3_x.drop(columns=['Name'])
df['pIC50']

df3 = pd.concat([df3_x, df['pIC50']], axis=1)     

df3.to_csv('../data/MMP1_3class_pIC50_pubchem_fp.csv')
