import pandas as pd
from chembl_webresource_client.new_client import new_client

#Neste projeto, usaremos compostos que possuem ação inibitória contra metaloproteinase-1 para posterior identificação
#da atividade inibidora de compostos naturais contidos em um subset da base de dados LOTUS

#Este projeto foi realizado para aprendizado do uso de:
    #1) Programação em python voltada as ferramentas de análise de compostos quimicos e análise do tipo QSAR
        #a) rdkit, sklearn
    #2) Entendimento do uso das bases chembl ou outras que contém dados de estruturas químicas, uso de códigos CAS e estrutura SMILES para identificação de moléculas
    #3) Testes de regressão e machine learning para classificação das moléculas
 
#======Análise inicial======#
    
#Obter dados de moléculas que possuem a MMP1 como alvo
    
target = new_client.target
target_query = target.search('MMP1')
targets = pd.DataFrame.from_dict(target_query)

     
#Capturar o primeiro resultado, de espécie humana
selected_target = targets.target_chembl_id[0]


#Com o indicador correto da molécula MMP1, obter os dados de atividade dos alvos, mantendo o dado de IC50
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type = 'IC50')

df = pd.DataFrame.from_dict(res)
df.to_csv('../data/MMP1_bioactivity_raw_data.csv', index=False)


#Filtrar para remover NAs e dplicados
df2 = df[df.standard_value.notna()]
df2 = df2[df2.canonical_smiles.notna()]

df2_nr = df2.drop_duplicates(['canonical_smiles'])
df2_nr.head()

df2_nr['units'].unique()


#Pré-processamento
df2_nr.columns

selection = ['molecule_chembl_id','canonical_smiles','standard_value', 'units']
df3 = df2_nr[selection]
df3.to_csv('MMP1_bioactivity_pre_processed.csv', index=False)

df3.columns
df3['units'].unique()

#Remover unidades estranhas ou ausentes
df4 = df3[~df3['units'].isin([None, 'None', 'ug ml-1'])]
df4['units'].unique()

#Converter unidade dos valores de IC50 para todos ficarem nM

def para_nM(valor, unidade):
    if unidade == "nM":
        return float(valor)
    if unidade == "uM":
        return float(valor) * 1000
    if unidade == "10^4nM":
        return float(valor) * 10**4

df4['corrected_value'] = []

df4['corrected_value'] = df4.apply(lambda row: para_nM(row['standard_value'], row['units']), axis=1)

bioactivity_threshold = []
for i in df4.corrected_value:
  if float(i) >= 10000:
    bioactivity_threshold.append("inactive")
  elif float(i) <= 1000:
    bioactivity_threshold.append("active")
  else:
    bioactivity_threshold.append("intermediate")
     
bioactivity_class = pd.Series(bioactivity_threshold, name='class')
df5 = pd.concat([df4, bioactivity_class], axis=1)
df5.to_csv('../data/MMP1_bioactivity_processed_data.csv')
