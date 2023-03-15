import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
import base64
from seqfold import dg
import itertools
from sklearn.model_selection import train_test_split
from sklearn import tree
import numpy as np
import pickle

st.header('アンチ選太くん')
id1 = st.sidebar.text_input('遺伝子名','Gapdh')
id2 = st.sidebar.text_input('Suffix（例：ASO(14)）','ASO(14)')
seq1 = st.sidebar.text_input('遺伝子配列を入力（例：AAATGGT...）','')
seq2 = st.sidebar.text_input('比較配列1（Ref1）を入力（例：AAATGGT...）','')
seq3 = st.sidebar.text_input('比較配列2（Ref2）を入力（例：AAATGGT...）','')

seq1 = Seq(seq1.upper())
seq2 = Seq(seq2.upper())
seq3 = Seq(seq3.upper())
numr = st.sidebar.slider('左翼の数', 0, 5,3)
gap = st.sidebar.slider('ギャップの数',0,15,8)
numl = st.sidebar.slider('右翼の数',0,5,3)
tail = st.sidebar.slider('テイルの数',0,5,0)
motif = st.sidebar.slider('トリプレット',0,5,3)

lenaso = int(numr)+int(gap)+int(numl)+int(tail)
mRNA = st.sidebar.slider('mRNA snippet',0,50,int(lenaso))

st.sidebar.write('ASO鎖長',lenaso)

st.write('query配列長:',len(seq1))



seq1r=str(seq1.reverse_complement())
seq1f=str(seq1)
seq2r=str(seq2.reverse_complement())
seq3r=str(seq3.reverse_complement())

#emptyリスト
list_df = pd.DataFrame( columns=["No",'ID','5wing',"ASO（5'to3'）",'3wing','GC%','Tm','deltaG', 'cpg in gap','hom vs ref1', 'hom vs ref2'] )
#リスト追加（Seqクラスはstrに直してから使う）
for i in range(len(seq1)-int(lenaso)+1):
    tmp_se = pd.Series( 
    [int(len(seq1))-(i+1)-(lenaso)+2,
    id1+"-"+str(int(len(seq1))-(i+1)-(lenaso)+2)+"-"+id2 ,  seq1r[i:i+int(numr)],
    seq1r[i:i+int(numr)]
    +seq1r[i+int(numr): i+int(numr)+int(gap)].lower()
    +seq1r[i+ int(numr)+ int(gap):i+ int(numr)+ int(gap)+ int(numl)]
    +seq1r[i+ int(numr)+ int(gap)+ int(numl):i+ int(numr)+ int(gap)+ int(numl)+int(tail)].lower(),
    seq1r[i+ int(numr)+ int(gap):i+ int(numr)+ int(gap)+ int(numl)],
    GC(seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))]),
    mt.Tm_NN(seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))]),
    dg(seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))],temp = 37.0),
     "cg" in seq1r[i+int(numr): i+int(numr)+int(gap)].lower(),
     seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))] in seq2r,
     seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))] in seq3r]
    
    ,index=list_df.columns )
    list_df = list_df.append( tmp_se, ignore_index=True )
#=========================機械学習用===========
gfv = []
tmv = []
dgv = []

for j in range(len(list_df)):
  gfv.append(gf(list_df["ASO（5'to3'）"].loc[j]))
  tmv.append(mt.Tm_NN(list_df["ASO（5'to3'）"].loc[j]))
  dgv.append(dg(list_df["ASO（5'to3'）"].loc[j],temp = 37.0))
  ddf = pd.DataFrame({'GC%': gfv, 'Tm':tmv, 'deltaG':dgv})
  dddf = list_df.join(ddf)

#トリプレット生成
base = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']
tlist = []

for trp in itertools.permutations(base, 3):
      tlist.append(trp[0]+trp[1]+trp[2])

#print(tlist)
tf = [] 
list_tf = pd.DataFrame()
dddf = pd.DataFrame()
for i in tlist:
  for j in range(len(list_df)):
    list_tf.loc[j,i] = i in str(list_df["ASO（5'to3'）"].loc[j])
#df2 = pd.concat([list_df, list_tf], axis =1)
fdf = dddf.join(list_tf)

x = fdf.loc[:, 'atg':'CGT']
t = fdf['ast']

#決定木
clf = pickle.load(open('yoshidamodel.pkl', 'rb'))

pred = clf.predict(x)
#tox = pd.DataFrame({'tox':pred})
#fdf2 = fdf.join(tox)

#=========================機械学習用===========ここまで

st.dataframe(fdf2.sort_values("No"))

list2_df = pd.DataFrame( columns=["No","snippet", "rev_compl"])

for i in range(len(seq1)-int(mRNA)+1):
    tmp2_se = pd.Series(
        [i+1, seq1f[i:i+int(mRNA)], str(seq1[i:i+int(mRNA)].reverse_complement())]
    ,index=list2_df.columns
    )
    list2_df = list2_df.append( tmp2_se, ignore_index=True )

st.dataframe(list2_df.sort_values("No"))

csv1 = df2.sort_values("No").to_csv(index=False) 
st.download_button('Download top table', csv1, 'text/csv')

csv2 = list2_df.sort_values("No").to_csv(index=False) 
st.download_button('Download bottom table', csv2, 'text/csv') 

#st.write("query配列（5'to3'）:",seq1)
#st.write("query逆相補鎖配列（5'to3'）: ",seq1.reverse_complement())

    
    

