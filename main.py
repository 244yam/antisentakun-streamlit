import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import GC

st.header('アンチ選択くん')
id1 = st.sidebar.text_input('遺伝子名','Gapdh')
seq1 = st.sidebar.text_input('主配列を入力','AAAAGGATATATATAGATAGAGATAGA')
seq2 = st.sidebar.text_input('参照配列1（Ref1）を入力','ATAGA')
seq3 = st.sidebar.text_input('参照配列2（Ref2）を入力','ATATAGATAGAGATAGA')

seq1 = Seq(seq1.upper())
seq2 = Seq(seq2.upper())
seq3 = Seq(seq3.upper())
numr = st.sidebar.text_input('右翼',3)
gap = st.sidebar.text_input('右翼',8)
numl = st.sidebar.text_input('左翼',3)
tail = st.sidebar.text_input('テイル',0)

st.sidebar.write('ASO鎖長',int(numr)+int(gap)+int(numl)+int(tail))
st.write('配列長:',len(seq1))
st.write('相補鎖:',seq1.reverse_complement())

seq1r=str(seq1.reverse_complement())
seq2r=str(seq2.reverse_complement())
seq3r=str(seq3.reverse_complement())

#emptyリスト
list_df = pd.DataFrame( columns=['ID','ASO','GC%','cpg in gap','hom vs ref1', 'hom vs ref2'] )
#リスト追加（Seqクラスはstrに直してから使う）
for i in range(len(seq1)-(int(numr)+int(gap)+int(numl)+int(tail)+1)):
    tmp_se = pd.Series( 
    [id1+"-"+str(i+1)+"-ASO" , 
    seq1r[i:i+int(numr)]
    +seq1r[i+int(numr): i+int(numr)+int(gap)].lower()
    +seq1r[i+ int(numr)+ int(gap):i+ int(numr)+ int(gap)+ int(numl)]
    +seq1r[i+ int(numr)+ int(gap)+ int(numl):i+ int(numr)+ int(gap)+ int(numl)+int(tail)].lower(),
    GC(seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))]),
     "cg" in seq1r[i+int(numr): i+int(numr)+int(gap)].lower(),
     seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))] in seq2r,
     seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))] in seq3r
    ]
    ,index=list_df.columns )
    list_df = list_df.append( tmp_se, ignore_index=True )

st.dataframe(list_df)
