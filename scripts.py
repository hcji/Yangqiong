# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:27:22 2019

@author: hcji
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

from chemopy.molproperty import GetMolecularProperty
from chemopy.connectivity import GetConnectivity
from chemopy.charge import GetCharge

# load data
subs = pd.read_excel('data/data1.xlsx')
reac = pd.read_excel('data/data2.xlsx')

# define functions
maccs = [MACCSkeys.GenMACCSKeys(Chem.MolFromInchi(i)) for i in subs['InChI']]
rdkfp = [Chem.rdmolops.RDKFingerprint(Chem.MolFromInchi(i)) for i in subs['InChI']]
prope = [GetMolecularProperty(Chem.MolFromInchi(i)) for i in subs['InChI']]
conne = [GetConnectivity(Chem.MolFromInchi(i)) for i in subs['InChI']]
charg = [GetCharge(Chem.MolFromInchi(i)) for i in subs['InChI']]

def condition_encoder(conditions):
    encoder = preprocessing.LabelEncoder()
    encoder.fit(conditions)
    hot = encoder.transform(conditions)
    le = len(encoder.classes_)
    output = []
    for j in range(len(hot)):
        i = np.zeros(le)
        i[hot[j]] = 1
        output.append(i)
    return np.array(output)

def fp_to_data(reactants, fps):
    output = []
    for r in reactants:
        wh = np.where(subs['reactant']==r)[0][0]
        output.append(np.array(fps[wh]))
    return np.array(output)    

def dp_to_data(reactants, dps):
    output = []
    for r in reactants:
        wh = np.where(subs['reactant']==r)[0][0]
        output.append(list(dps[wh].values()))
    return np.array(output)   

# get descriptors
cata_code = condition_encoder(reac['potocatalyst'])
base_code = condition_encoder(reac['base'])
solv_code = condition_encoder(reac['solvents'])
react_maccs = fp_to_data(reac['reactant'], maccs)
react_rdkfp = fp_to_data(reac['reactant'], rdkfp)
react_prope = dp_to_data(reac['reactant'], prope)
react_conne = dp_to_data(reac['reactant'], conne)
react_charg = dp_to_data(reac['reactant'], charg)

# concatenate all descriptors
Yield = np.array(reac['yield'])
Xdata = np.concatenate((cata_code, base_code, solv_code, react_maccs, react_rdkfp, react_prope, react_conne, react_charg), axis=1)

# split data
X_train, X_test, y_train, y_test = train_test_split(Xdata, Yield, test_size=0.2)

# scale data
scaler = MinMaxScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

# train random forest regression model
w = []
for y in y_train:
    if y == 0:
        w.append(0.5)
    else:
        w.append(1)
rf_regr = RandomForestRegressor(max_depth=20, n_estimators=500)
rf_regr.fit(X_train, y_train, np.array(w))
y_pred = rf_regr.predict(X_test)

# print output
r2 = r2_score(y_test, y_pred)
print('r2:' + str(round(r2, 4)))
plt.figure(figsize=(6, 6))
plt.plot(y_test, y_pred, 'bo')
plt.plot([0, max(y_test) + 2], [0,max(y_test) + 2], color ='red')
plt.xlabel('yield')
plt.ylabel('prediction')
plt.show()