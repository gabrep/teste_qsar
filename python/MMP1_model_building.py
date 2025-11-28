import pandas as pd
import seaborn as sns
import matplotlib as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold

#Utilizando os descritores criados, criar um modelo para predição de IC50

df = pd.read_csv('../data/MMP1_3class_pIC50_pubchem_fp.csv', index_col=0)

X = df.drop('pIC50', axis=1)
X

Y = df.pIC50
Y

X.shape
Y.shape

#Remover baixa variancia
#Parecido com o principio para calculo de DEGs
#Remover o que nao classifica as moleculas para manter os descritores mais relevantes
selection = VarianceThreshold(threshold=(.8 * (1 - .8)))    
selection.fit(X)
mask = selection.get_support()
X = X.loc[:, mask] 
##Dividir o dataset para criar o dataset de treino e de validação
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)


X_train.shape, Y_train.shape
X_test.shape, Y_test.shape


##===== Random Forest ======#
import numpy as np
np.random.seed(111)
model = RandomForestRegressor(n_estimators=200)

#Lembrando que X = descritores, Y = pIC50
model.fit(X_train, Y_train)
r2 = model.score(X_test, Y_test)
r2
#r quadrado de 0.75 nos dados de teste


Y_pred = model.predict(X_test)


sns.set(color_codes=True)
sns.set_style("white")

ax = sns.regplot(x=Y_test,y =  Y_pred, scatter_kws={'alpha':0.4})
ax.set_xlabel('Experimental pIC50', fontsize='large', fontweight='bold')
ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.figure.set_size_inches(5, 5)
plt.show


#====Teste de outros models====#
import lazypredict
from lazypredict.Supervised import LazyRegressor

clf = LazyRegressor(verbose=0,ignore_warnings=True, custom_metric=None)
models_train,predictions_train = clf.fit(X_train, X_train, Y_train, Y_train)
models_test,predictions_test = clf.fit(X_train, X_test, Y_train, Y_test)

predictions_test.to_csv('../data/Prediction_models_test.csv')
predictions_train.to_csv('../data/Prediction_models_train.csv')

#===========================================#
#=====Rodar o modelo com os dados Lotus=====#
#===========================================#
from padelpy import padeldescriptor
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold

#Dados baixados da plataforma LOTUS
#Ler o lotus.smi mas apenas as primeiras 15000 entradas, para reduzir o tamanho do arquivo
df_lotus = pd.read_csv('../LOTUS_DB.smi', sep="\t", header=None, nrows=15000)
df_lotus.to_csv('../LOTUS_DB_reduzido.smi', sep='\t', index=False, header=False)

padeldescriptor(mol_dir='../LOTUS_DB_reduzido.smi', d_file='../descriptors_LOTUS.csv', removesalt=True, fingerprints=True)

#Carregar os descritores criados para as moleculas (subset) da base LOTUS
lotus = pd.read_csv('../descriptors_LOTUS.csv')

lotus= pd.read_csv('../LOTUS_DB_reduzido.smi', sep='\t', header=None, )
lotus[0]
lotus_desc = pd.read_csv('../Lotus_descritores.csv', index_col=0)

#Manter para as moleculas LOTUS apenas os descritores contidos no modelo de treinamento
lotus_aligned = lotus_desc[X_test.columns]

#Predição de pIC50
lotus_pred = model.predict(lotus_aligned)

from rdkit.Chem import Descriptors, Lipinski
from rdkit import Chem

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

#Obter descritores de Lipinski para as moleculas LOTUS
lotus_lipinski = lipinski(lotus[0])

lotus_p = pd.DataFrame()
lotus_p['pred_pIC50']= []
lotus_p['pred_pIC50']=lotus_pred


lotus_final = pd.concat([lotus_p, lotus_lipinski], axis=1)
#Adicionar dados estruturais das moleculas
lotus_final = pd.concat([lotus_final,lotus[0]], axis=1)

#Converter de pIC50 para IC50
lotus_final['IC50'] = 10**(-lotus_final.pred_pIC50)

from rdkit.Chem import Draw
Draw.MolToImage(Chem.MolFromSmiles(lotus.loc[11134,0]))
lotus.loc[1,1]

10**(-lotus_final.pred_pIC50)

lotus_final.to_csv('../data/Lotus_predito_final.csv')

#Desenha as 10 moleculas com menor IC50
lotus_top10 = lotus_final.sort_values(by='IC50').iloc[0:10]
lotus_top10.to_csv('../data/Lotus_predito_top10.csv')
img = [Draw.MolToImage(Chem.MolFromSmiles(sm)) for sm in lotus_top10.iloc[0:10,7]]

#Chat criou o loop para exportar as figuras de cada uma das top 10, nomeadas com o codigo lotus
from pathlib import Path
out = Path("../figures")

for i in range(10):
    smiles = lotus_top10.iloc[i, 7]        # coluna SMILES
    name   = str(lotus_top10.iloc[i, 6])   # nome da molécula

    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    filename = f"{i+1}_{name}.png"   # garante nomes únicos

    img.save(out / filename)
