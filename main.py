import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import base64

st.header('アンチ選択くん')
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
mRNA = st.sidebar.slider('mRNA snippet',0,50,25)
lenaso = int(numr)+int(gap)+int(numl)+int(tail)

st.sidebar.write('ASO鎖長',lenaso)
st.write('query配列長:',len(seq1))

seq1r=str(seq1.reverse_complement())
seq1f=str(seq1.reverse())
seq2r=str(seq2.reverse_complement())
seq3r=str(seq3.reverse_complement())

#emptyリスト
list_df = pd.DataFrame( columns=["No",'ID','5wing',"ASO（5'to3'）",'3wing','GC%','cpg in gap','hom vs ref1', 'hom vs ref2'] )
#リスト追加（Seqクラスはstrに直してから使う）
for i in range(len(seq1)-(int(numr)+int(gap)+int(numl)+int(tail))+1):
    tmp_se = pd.Series( 
    [int(len(seq1))-(i+1)-(lenaso)+2,
    id1+"-"+str(int(len(seq1))-(i+1)-(lenaso)+2)+"-"+id2 ,  seq1r[i:i+int(numr)],
    seq1r[i:i+int(numr)]
    +seq1r[i+int(numr): i+int(numr)+int(gap)].lower()
    +seq1r[i+ int(numr)+ int(gap):i+ int(numr)+ int(gap)+ int(numl)]
    +seq1r[i+ int(numr)+ int(gap)+ int(numl):i+ int(numr)+ int(gap)+ int(numl)+int(tail)].lower(),
    seq1r[i+ int(numr)+ int(gap):i+ int(numr)+ int(gap)+ int(numl)],
    GC(seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))]),
     "cg" in seq1r[i+int(numr): i+int(numr)+int(gap)].lower(),
     seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))] in seq2r,
     seq1r[i:i+(int(numr)+int(gap)+int(numl)+int(tail))] in seq3r]
    ,index=list_df.columns )
    list_df = list_df.append( tmp_se, ignore_index=True )

st.dataframe(list_df.sort_values("No"))

list2_df = pd.DataFrame( columns=["No","snippet"])

for i in range(len(seq1)-int(mRNA)+1)
    tmp2_se = pd.Series(
        [i+1, seq1f[i:i+int(mRNA)]]
    ,index=list2_df.columns
    )
    list2_df = list2_df.append( tmp2_se, ignore_index=True )

st.dataframe(list2_df.sort_values("No"))

csv = list_df.sort_values("No").to_csv(index=False)  
b64 = base64.b64encode(csv.encode()).decode()
href = f'<a href="data:application/octet-stream;base64,{b64}" download="result.csv">download</a>'
st.markdown(f" CSVとしてダウンロードする： {href}", unsafe_allow_html=True)

csv2 = list2_df.sort_values("No").to_csv(index=False)  
b65 = base64.b64encode(csv.encode()).decode()
href = f'<a href="data:application/octet-stream;base64,{b65}" download="result2.csv">download</a>'
st.markdown(f" CSVとしてダウンロードする： {href}", unsafe_allow_html=True)

st.write("query配列（5'to3'）:",seq1)
st.write("query逆相補鎖配列（5'to3'）: ",seq1.reverse_complement())
