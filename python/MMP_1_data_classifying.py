import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import Draw

#Criar classificadores de Lipinski para as moléculas e trabalhar a unidade de IC50

df = pd.read_csv('../data/MMP1_bioactivity_processed_data.csv', index_col=0)
df = df.dropna()

##LIPINSKI DESCRIPTORS
#Descritores que indicam a atividade do composto
#baseado nas propriedades cineticas farmacologicas:
#absorção, distribuição, metabolismo e excreção

#Molecular weight < 500 Dalton
#Octanol-water partition coefficient (LogP) < 5
#Hydrogen bond donors < 5
#Hydrogen bond acceptors < 10


#====Lipinski Descriptors====#
# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
#Fonte: https://www.youtube.com/watch?v=jBlTQjcKuaY

def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

#Desenho bacana das moléculas com base na estrutura SMILES
Draw.MolToImage(Chem.MolFromSmiles(df.canonical_smiles[2]))

df_lipinski = lipinski(df.canonical_smiles)


df_combined = pd.concat([df,df_lipinski], axis=1)
df_combined = df_combined.dropna()

df.shape
df_lipinski.shape

#=========Transformação IC50=======
#Transoformar IC50 para pIC50
#Logaritmo negativo, para corrigir a distribuição do IC50

# https://github.com/chaninlab/estrogen-receptor-alpha-qsar/blob/master/02_ER_alpha_RO5.ipynb
def pIC50(input):
    pIC50 = []

    for i in input['corrected_value']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('corrected_value', axis = 1)
        
    return x

df_combined.standard_value.describe()

#Se valore de IC50 estiverem acima de 10^8, a transformação em log negativo retorna valores abaixo de 0, que dificulta a interpretação
#pode se criar um limite máximo de valores para substituir os valores acima disso
(df_combined.standard_value > 100000000).sum()
#Como nenhum valor é maior que isto, nao precisa normalizar para limite maximo

df_final = pIC50(df_combined)
df_final.pIC50.describe()
df_final.to_csv('../data/MMP1_bioactivity_3class_pIC50.csv')

df_2class = df_final[df_final['class'] != 'intermediate']
df_2class.to_csv('../data/MMP1_bioactivity_2class_pIC50.csv')
     