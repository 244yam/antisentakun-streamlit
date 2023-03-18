! pip install biopython
! pip install seqfold
import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction as gf
from Bio.SeqUtils import MeltingTemp as mt
import base64
from seqfold import dg
import itertools
from sklearn.model_selection import train_test_split
from sklearn import tree
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

st.header('アンチ選太くん２')
df =st.sidebar.file_uploader("ファイルアップロード (LLLdddddLLLなどのリスト)", type='csv')
if df is not None:
    df = pd.read_csv(df)
    st.dataframe(df)

gfv = []
tmv = []
dgv = []

for j in range(len(df)):
  gfv.append(gf(df['aso'].loc[j]))
  tmv.append(mt.Tm_NN(df['aso'].loc[j]))
  dgv.append(dg(df['aso'].loc[j],temp = 37.0))
  ddf = pd.DataFrame({'GC%': gfv, 'Tm':tmv, 'deltaG':dgv})
  dddf = df.join(ddf)
dddf = dddf.replace([np.inf, -np.inf], np.nan)
dddf = dddf.fillna(100) # 0でnanを置換


#トリプレット生成
base = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']
tlist = []

for trp in itertools.product(base, repeat=3):
      tlist.append(trp[0]+trp[1]+trp[2])

#print(tlist)
tf = [] 
list_tf = pd.DataFrame()

for i in tlist:
  for j in range(len(df)):
    list_tf.loc[j,i] = i in str(df['aso'].loc[j])
#df2 = pd.concat([df, list_tf], axis =1)
fdf = dddf.join(list_tf)

x = fdf.loc[:, 'GC%':'CCC']

ml_menu = st.selectbox("実施する予測モデルを選択してください", ["DecisionTree","RandomForest"])
if ml_menu == "DecisionTree":
    st.markdown("毒性を予測します")
    execute = st.button("実行")
    
    if execute:
        with open('yoshidaBurdick_ki_model.pkl', 'rb') as f:
            ktg = pickle.load(f)
        pred = ktg.predict(x)
    tox = pd.DataFrame({'tox_ki':pred})
    alldf = fdf.join(tox)
    csv1 = alldf.to_csv(index=False) 
    st.download_button(label="Download data as CSV",
    data=csv1,
    file_name='predict_tree.csv',
    mime='text/csv')

elif ml_menu == "RandomForest":
    st.markdown("毒性を予測します")
    execute = st.button("実行")
  
    if execute:
        with open('yoshidaBurdick_mori_model.pkl', 'rb') as f:
            ktg2 = pickle.load(f)  
        pred2 = ktg2.predict(x)
    tox2 = pd.DataFrame({'tox_mori':pred2})
    alldf = fdf.join(tox2)
    csv2 = alldf.to_csv(index=False) 
    st.download_button(label="Download data as CSV",
    data=csv2,
    file_name='predict_forest.csv',
    mime='text/csv')
st.dataframe(alldf)


