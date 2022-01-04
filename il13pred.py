###################################################################################
# IL13Pred is developed for predicting, desigining and scanning interleukin-13    #
# inducing peptides. It is developed by Prof G. P. S. Raghava's group.            #
# Please cite: IL13Pred; available at https://webs.iiitd.edu.in/raghava/il13pred/ #
##################################################################################
import argparse
import warnings
import subprocess
import pkg_resources
import os
import sys
import numpy as np
import pandas as pd
import math
import itertools
from collections import Counter
import pickle
import re
import glob
import time
from time import sleep
from tqdm import tqdm
import xgboost
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Please provide following arguments') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3], help="Job Type: 1:predict, 2:design and 3:scan, by default 1")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.06")
parser.add_argument("-w","--winleng", type=int, choices =range(8, 36), help="Window Length: 8 to 35 (scan mode only), by default 9")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Interleukin-13 inducing peptide, 2: All peptides, by default 1")
args = parser.parse_args()

# Function for generating all possible mutants
def seq_mutants(aa_seq):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    aa_seq = aa_seq.upper()
    mut_seq=[]
    for pos in range(0,len(aa_seq)):
        for aa in std:
            mut_seq += [aa_seq[:pos] + aa + aa_seq[pos+1:]]
    return mut_seq
# Function for generating pattern of a given length
def seq_pattern(aa_seq,win_len):
    aa_seq == aa_seq.upper
    seq_pat=[]
    for i1 in range(0, (len(aa_seq) + 1 - win_len)):
        i2 = i1 + int(win_len)
        seq_pat += [aa_seq[i1:i2]]
    return seq_pat

def feature_gen(file):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    sf = ['BTC_T','SOC1_G1','ATC_N','AAC_L','CeTD_HB1','AAC_V','AAC_F','AAC_G','CeTD_SA3','AAC_A']
    def aac_comp(file):
        filename, file_extension = os.path.splitext(file)
        f = open('sam_allcomp.aac', 'w')
        sys.stdout = f
        df = pd.read_csv(file, header = None)
        zz = df.iloc[:,0]
        print("AAC_A,AAC_C,AAC_D,AAC_E,AAC_F,AAC_G,AAC_H,AAC_I,AAC_K,AAC_L,AAC_M,AAC_N,AAC_P,AAC_Q,AAC_R,AAC_S,AAC_T,AAC_V,AAC_W,AAC_Y,")
        for j in zz:
            for i in std:
                count = 0
                for k in j:
                    temp1 = k
                    if temp1 == i:
                        count += 1
                    composition = (count/len(j))*100
                print("%.2f"%composition, end = ",")
            print("")
        f.truncate()
    ###########################atom###############
    def atc(file):
        filename,file_ext = os.path.splitext(file)
        atom=pd.read_csv("Data/atom.csv",header=None)
        at=pd.DataFrame()
        i = 0
        C_atom = []
        H_atom = []
        N_atom = []
        O_atom = []
        S_atom = []

        while i < len(atom):
            C_atom.append(atom.iloc[i,1].count("C"))
            H_atom.append(atom.iloc[i,1].count("H"))
            N_atom.append(atom.iloc[i,1].count("N"))
            O_atom.append(atom.iloc[i,1].count("O"))
            S_atom.append(atom.iloc[i,1].count("S"))
            i += 1
        atom["C_atom"]=C_atom
        atom["O_atom"]=O_atom
        atom["H_atom"]=H_atom
        atom["N_atom"]=N_atom
        atom["S_atom"]=S_atom
    ##############read file ##########
        test1 = pd.read_csv(file,header=None)
        dd = []
        for i in range(0, len(test1)):
            dd.append(test1[0][i].upper())
        test = pd.DataFrame(dd)
        count_C = 0
        count_H = 0
        count_N = 0
        count_O = 0
        count_S = 0
        count = 0
        i1 = 0
        j = 0
        k = 0
        C_ct = []
        H_ct = []
        N_ct = []
        O_ct = []
        S_ct = []
        while i1 < len(test) :
            while j < len(test[0][i1]) :
                while k < len(atom) :
                    if test.iloc[i1,0][j]==atom.iloc[k,0].replace(" ","") :
                        count_C = count_C + atom.iloc[k,2]
                        count_H = count_H + atom.iloc[k,3]
                        count_N = count_N + atom.iloc[k,4]
                        count_O = count_O + atom.iloc[k,5]
                        count_S = count_S + atom.iloc[k,6]
                    #count = count_C + count_H + count_S + count_N + count_O
                    k += 1
                k = 0
                j += 1
            C_ct.append(count_C)
            H_ct.append(count_H)
            N_ct.append(count_N)
            O_ct.append(count_O)
            S_ct.append(count_S)
            count_C = 0
            count_H = 0
            count_N = 0
            count_O = 0
            count_S = 0
            j = 0
            i1 += 1
        test["C_count"]=C_ct
        test["H_count"]=H_ct
        test["N_count"]=N_ct
        test["O_count"]=O_ct
        test["S_count"]=S_ct

        ct_total = []
        m = 0
        while m < len(test) :
            ct_total.append(test.iloc[m,1] + test.iloc[m,2] + test.iloc[m,3] + test.iloc[m,4] + test.iloc[m,5])
            m += 1
        test["count"]=ct_total
    ##########final output#####
        final = pd.DataFrame()
        n = 0
        p = 0
        C_p = []
        H_p = []
        N_p = []
        O_p = []
        S_p = []
        while n < len(test):
            C_p.append((test.iloc[n,1]/test.iloc[n,6])*100)
            H_p.append((test.iloc[n,2]/test.iloc[n,6])*100)
            N_p.append((test.iloc[n,3]/test.iloc[n,6])*100)
            O_p.append((test.iloc[n,4]/test.iloc[n,6])*100)
            S_p.append((test.iloc[n,5]/test.iloc[n,6])*100)
            n += 1
        final["ATC_C"] = C_p
        final["ATC_H"] = H_p
        final["ATC_N"] = N_p
        final["ATC_O"] = O_p
        final["ATC_S"] = S_p

        (final.round(2)).to_csv('sam_allcomp.atc', index = None, encoding = 'utf-8')
    ########################################atom#################
    def bond(file) :
        tota = []
        hy = []
        Si = []
        Du = []
        b1 = []
        b2 = []
        b3 = []
        b4 = []
        bb = pd.DataFrame()
        filename, file_extension = os.path.splitext(file)
        df = pd.read_csv(file, header = None)
        bonds=pd.read_csv("Data/bonds.csv", sep = ",")
        for i in range(0,len(df)) :
            tot = 0
            h = 0
            S = 0
            D = 0
            tota.append([i])
            hy.append([i])
            Si.append([i])
            Du.append([i])
            for j in range(0,len(df[0][i])) :
                temp = df[0][i][j]
                for k in range(0,len(bonds)) :
                    if bonds.iloc[:,0][k] == temp :
                        tot = tot + bonds.iloc[:,1][k]
                        h = h + bonds.iloc[:,2][k]
                        S = S + bonds.iloc[:,3][k]
                        D = D + bonds.iloc[:,4][k]
            tota[i].append(tot)
            hy[i].append(h)
            Si[i].append(S)
            Du[i].append(D)
        for m in range(0,len(df)) :
            b1.append(tota[m][1])
            b2.append(hy[m][1])
            b3.append(Si[m][1])
            b4.append(Du[m][1])

        bb["BTC_T"] = b1
        bb["BTC_H"] = b2
        bb["BTC_S"] = b3
        bb["BTC_D"] = b4

        bb.to_csv('sam_allcomp.btc', index=None, encoding="utf-8")
    ######################################CTD###################################
    def ctd(file):
        attr=pd.read_csv("Data/aa_attr_group.csv", sep="\t")
        filename, file_extension = os.path.splitext(file)
        df1 = pd.read_csv(file, header = None)
        df = pd.DataFrame(df1[0].str.upper())
        n = 0
        stt1 = []
        m = 1
        for i in range(0,len(attr)) :
            st =[]
            stt1.append([])
            for j in range(0,len(df)) :
                stt1[i].append([])
                for k in range(0,len(df[0][j])) :
                    while m < 4 :
                        while n < len(attr.iloc[i,m]) :
                            if df[0][j][k] == attr.iloc[i,m][n] :
                                st.append(m)
                                stt1[i][j].append(m)
                            n += 2
                        n = 0
                        m += 1
                    m = 1
    #####################Composition######################
        f = open("compout_1", 'w')
        sys.stdout = f
        std = [1,2,3]
        print("1,2,3,")
        for p in range (0,len(df)) :
            for ii in range(0,len(stt1)) :
                #for jj in stt1[ii][p]:
                for pp in std :
                    count = 0
                    for kk in stt1[ii][p] :
                        temp1 = kk
                        if temp1 == pp :
                            count += 1
                        composition = (count/len(stt1[ii][p]))*100
                    print("%.2f"%composition, end = ",")
                print("")
        f.truncate()

    #################################Transition#############
        tt = []
        tr=[]
        kk =0
        for ii in range(0,len(stt1)) :
            tt = []
            tr.append([])
            for p in range (0,len(df)) :
                tr[ii].append([])
                while kk < len(stt1[ii][p]) :
                    if kk+1 <len(stt1[ii][p]):
                    #if  stt1[ii][p][kk] < stt1[ii][p][kk+1] or stt1[ii][p][kk] > stt1[ii][p][kk+1]: # condition for adjacent values
                        tt.append(stt1[ii][p][kk])
                        tt.append(stt1[ii][p][kk+1])
                        tr[ii][p].append(stt1[ii][p][kk])
                        tr[ii][p].append(stt1[ii][p][kk+1])

                    kk += 1
                kk = 0

        pp = 0
        xx = []
        xxx = []
        for mm in range(0,len(tr)) :
            xx = []
            xxx.append([])
            for nn in range(0,len(tr[mm])):
                xxx[mm].append([])
                while pp < len(tr[mm][nn]) :
                    xx .append(tr[mm][nn][pp:pp+2])
                    xxx[mm][nn].append(tr[mm][nn][pp:pp+2])
                    pp+=2
                pp = 0

        f1 = open("compout_2", 'w')
        sys.stdout = f1
        std1 = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
        print("1->1,1->2,1->3,2->1,2->2,2->3,3->1,3->2,3->3,")
        for rr in range(0,len(df)) :
            for qq in range(0,len(xxx)):
                for tt in std1 :
                    count = 0
                    for ss in xxx[qq][rr] :
                        temp2 = ss
                        if temp2 == tt :
                            count += 1
                    print(count, end = ",")
                print("")
        f1.truncate()

        #################################Distribution#############
        c_11 = []
        c_22 = []
        c_33 = []
        zz = []
        #print("0% 25% 50% 75% 100%")
        for x in range(0,len(stt1)) :
            #c_11.append([])
            c_22.append([])
            #c_33.append([])
            yy_c_1 = []
            yy_c_2 = []
            yy_c_3 = []
            ccc = []

            k = 0
            j = 0
            for y in range(0,len(stt1[x])):
                #c_11[x].append([])
                c_22[x].append([])
                for i in range(1,4) :
                    cc = []
                    c1 = [index for index,value in enumerate(stt1[x][y]) if value == i]
                    c_22[x][y].append(c1)
        cc = []
        for ss in range(0,len(df)):
            for uu in range(0,len(c_22)):
                for mm in range(0,3):
                    for ee in range(0,101,25):
                        k = (ee*(len(c_22[uu][ss][mm])))/100
                        cc.append(math.floor(k))
        f2 = open('compout_3', 'w')
        sys.stdout = f2
        print("0% 25% 50% 75% 100%")
        for i in range (0,len(cc),5):
            print(*cc[i:i+5])
        f2.truncate()
        head = []
        header1 = ['CeTD_HB','CeTD_VW','CeTD_PO','CeTD_PZ','CeTD_CH','CeTD_SS','CeTD_SA']
        for i in header1:
            for j in range(1,4):
                head.append(i+str(j))
        df11 = pd.read_csv("compout_1")
        df_1 = df11.iloc[:,:-1]
        zz = pd.DataFrame()
        for i in range(0,len(df_1),7):
            zz = zz.append(pd.DataFrame(pd.concat([df_1.loc[i],df_1.loc[i+1],df_1.loc[i+2],df_1.loc[i+3],df_1.loc[i+4],df_1.loc[i+5],df_1.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
        zz.columns = head
        #zz.to_csv(filename+".ctd_comp", index=None, encoding='utf-8')
        head2 = []
        header2 = ['CeTD_11','CeTD_12','CeTD_13','CeTD_21','CeTD_22','CeTD_23','CeTD_31','CeTD_32','CeTD_33']
        for i in header2:
            for j in ('HB','VW','PO','PZ','CH','SS','SA'):
                head2.append(i+'_'+str(j))
        df12 = pd.read_csv("compout_2")
        df_2 = df12.iloc[:,:-1]
        ss = pd.DataFrame()
        for i in range(0,len(df_2),7):
            ss = ss.append(pd.DataFrame(pd.concat([df_2.loc[i],df_2.loc[i+1],df_2.loc[i+2],df_2.loc[i+3],df_2.loc[i+4],df_2.loc[i+5],df_2.loc[i+6]],axis=0)).transpose()).reset_index(drop=True)
        ss.columns = head2
        #ss.to_csv(filename+".ctd_trans", index=None, encoding='utf-8')
        head3 = []
        header3 = ['CeTD_0_p','CeTD_25_p','CeTD_50_p','CeTD_75_p','CeTD_100_p']
        header4 = ['HB','VW','PO','PZ','CH','SS','SA']
        for j in range(1,4):
            for k in header4:
                for i in header3:
                    head3.append(i+'_'+k+str(j))
        df_3 = pd.read_csv("compout_3", sep=" ")
        rr = pd.DataFrame()
        for i in range(0,len(df_3),21):
            rr = rr.append(pd.DataFrame(pd.concat([df_3.loc[i],df_3.loc[i+1],df_3.loc[i+2],df_3.loc[i+3],df_3.loc[i+4],df_3.loc[i+5],df_3.loc[i+6],df_3.loc[i+7],df_3.loc[i+8],df_3.loc[i+9],df_3.loc[i+10],df_3.loc[i+11],df_3.loc[i+12],df_3.loc[i+13],df_3.loc[i+14],df_3.loc[i+15],df_3.loc[i+16],df_3.loc[i+17],df_3.loc[i+18],df_3.loc[i+19],df_3.loc[i+20]],axis=0)).transpose()).reset_index(drop=True)
        rr.columns = head3
        cotrdi= pd.concat([zz,ss,rr],axis=1)
        cotrdi.to_csv('sam_allcomp.ctd', index=None, encoding='utf-8')
    ##################################SOCN#################################
    def soc(file):
        gap = 1
        ff = []
        filename, file_extension = os.path.splitext(file)
        df = pd.read_csv(file, header = None)
        df2 = pd.DataFrame(df[0].str.upper())
        for i in range(0,len(df2)):
            ff.append(len(df2[0][i]))
        if min(ff) < gap:
            print("Error: All sequences' length should be higher than :", gap)
            return 0
        mat1 = pd.read_csv("Data/Schneider-Wrede.csv", index_col = 'Name')
        mat2 = pd.read_csv("Data/Grantham.csv", index_col = 'Name')
        h1 = []
        h2 = []
        for n in range(1, gap+1):
            h1.append('SC' + str(n))
        for n in range(1, gap + 1):
            h2.append('G' + str(n))
        h1 = ['SOC'+str(gap)+'_'+sam for sam in h1]
        h2 = ['SOC'+str(gap)+'_'+sam for sam in h2]
        s1 = []
        s2 = []
        for i in range(0,len(df2)):
            for n in range(1, gap+1):
                sum = 0
                sum1 =0
                sum2 =0
                sum3 =0
                for j in range(0,(len(df2[0][i])-n)):
                    sum = sum + (mat1[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum1 = sum/(len(df2[0][i])-n)
                    sum2 = sum2 + (mat2[df2[0][i][j]][df2[0][i][j+n]])**2
                    sum3 = sum2/(len(df2[0][i])-n)
                s1.append(sum1)
                s2.append(sum3)
        zz = np.array(s1).reshape(len(df2),gap)
        zz2 = np.array(s2).reshape(len(df2),gap)
        zz3 = round(pd.concat([pd.DataFrame(zz, columns = h1),pd.DataFrame(zz2,columns = h2)], axis = 1),4)
        zz3.to_csv('sam_allcomp.soc', index = None, encoding = 'utf-8')
    ####################Function Calling#####################################
    aac_comp(file)
    atc(file)
    bond(file)
    ctd(file)
    soc(file)
    ##############################################################################
    df1 = pd.read_csv("sam_allcomp.aac")
    df4 = pd.read_csv("sam_allcomp.atc")
    df5 = pd.read_csv("sam_allcomp.btc")
    df14 = pd.read_csv("sam_allcomp.ctd")
    df15 = pd.read_csv("sam_allcomp.soc")
    df19 = pd.concat([df1.iloc[:,:-1],df4,df5,df14,df15],axis=1)
    filelist=glob.glob("sam_allcomp*")
    for file_2 in filelist:
        os.remove(file_2)
    os.remove('compout_1')
    os.remove('compout_2')
    os.remove('compout_3')
    df4 = df19[sf]
    df4.columns = ['BTC_T','SOC1_G1','ATC_N','AAC_L','CeTD_HB1','AAC_V','AAC_F','AAC_G','CeTD_SA3','AAC_A']
    return df4

def adjusted_classes(y_scores, t):
    return [1 if y >= t else 0 for y in y_scores]

def Perform_testing(clf,name,X,t):
    Y_pred = clf.predict(X)
    Y_scores=[]
    Y_scores=clf.predict_proba(X)[:,-1]
    Y_pred = adjusted_classes(Y_scores,t)
    return Y_pred,Y_scores

print('####################################################################################')
print('# This program IL13Pred is developed for predicting, desigining and scanning       #')
print('# interleukin-13 inducing peptides, developed by Prof G. P. S. Raghava group.      #')
print('# Please cite: IL13Pred; available at https://webs.iiitd.edu.in/raghava/il13pred/  #')
print('####################################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Threshold 
if args.threshold == None:
        Threshold = 0.11
else:
        Threshold= float(args.threshold)
# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)
# Window Length 
if args.winleng == None:
        Win_len = int(9)
else:
        Win_len = int(args.winleng)

# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)


if Job==2:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job)
    print('Output File: ',result_filename,'; Window Length: ',Win_len,'; Display: ',dplay)
    print('##############################################################################')
else:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job)
    print('Output File: ',result_filename,'; Display: ',dplay)
    print('# ############################################################################')
#------------------ Read input file ---------------------
def load_model(path):
    clf = pickle.load(open(path,'rb'))
    return clf

f=open(Sequence,"r")
len1 = f.read().count('>')
f.close()

f=open(Sequence,"r")
seqs=[]
seqid=[]
str1='';
header=0
line=0
if len1 >= 1: # read fasta file
    for l in f:
        if l.startswith('>'):
            if header != 0:
                seqs += [str1]
                str1 = '';
            header = 1
            line+=1
            seqid += [l.rstrip()]
        else:
            str1 += l.rstrip()
    seqs += [str1]
else: # read single line file    
    for l in f:
        if len(l) >= 8 and len(l)<= 25:
            seqs+=[l.strip()]
            seqid += ['>Seq_' + str(line)]
        line+=1
f.close()

fout= open(result_filename,"w+")

i1 = 0
#======================= Prediction Module start from here =====================
if Job == 1:
    print('\n======= Thanks for using Predict module of IL13Pred. Your results will be stored in file :',result_filename,' =====\n')
    fout.write('# Sequence_ID, Sequence, Prediction, Score\n')
    clf=load_model('XGB_model')
    for Sequence in tqdm(seqs):
        header = seqid[i1]
        if len(Sequence) >= 8: 
            if len(Sequence) >= 35:
                Sequence = Sequence[0:35]
            ss = []
            ss.append(Sequence)
            pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
            X = feature_gen(Sequence)
            Y_pred,Y_score=Perform_testing(clf,'XGB',X,Threshold)
            Y_score = np.round(Y_score,2)
            flag=""
            if Y_pred[0]==1:
                flag='IL-13 inducer'
            else:
                flag='IL-13 non-inducer'
            if dplay == 1:
                if Y_pred[0]==1:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
            else:
                fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
        os.remove(Sequence)
        i1 = i1 +1
#===================== Design Model Start from Here ======================
elif Job == 2:
    print('\n======= Thanks for using Design module of IL13Pred. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Designing Peptides: Processing sequences please wait ...')
    fout.write('# Sequence_ID, pattern, Prediction, Score\n')
    i1 = 0
    for Sequence1 in tqdm(seqs):
        pat_seq=[]
        header = seqid[i1]
#        print('Sequence #: ',i1,'# under process ...')
        if len(Sequence1) >= 8: 
            if len(Sequence1) >= 35:
                Sequence1 = Sequence1[0:35]
            pat_seq = seq_mutants(Sequence1)
            for Sequence in pat_seq:
                clf=load_model('XGB_model')
                ss = []
                ss.append(Sequence)
                pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
                X = feature_gen(Sequence)
                Y_pred,Y_score=Perform_testing(clf,'XGB',X,Threshold)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='IL-13 inducer'
                else:
                    flag='IL-13 non-inducer'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                os.remove(Sequence)
        i1 = i1 + 1
#=============== Scan Model start from here ==================
else:
    print('\n======= Thanks for using Scan module of IL13Pred. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Scanning Peptides: Processing sequences please wait ...')
    i1 = 0
    fout.write('# Sequence_ID, pattern, Prediction, Score\n')
    for Sequence1 in tqdm(seqs):
        pat_seq=[]
        header = seqid[i1]
        if len(Sequence1) >= Win_len: 
            pat_seq = seq_pattern(Sequence1,Win_len)
            for Sequence in pat_seq:
                clf=load_model('XGB_model')
                ss = []
                ss.append(Sequence)
                pd.DataFrame(ss).to_csv(Sequence,index=None,header=False)
                X = feature_gen(Sequence)
                Y_pred,Y_score=Perform_testing(clf,'XGB',X,Threshold)
                Y_score = np.round(Y_score,2)
                flag=""
                if Y_pred[0]==1:
                    flag='IL-13 inducer'
                else:
                    flag='IL-13 non-inducer'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%.2f\n" % (header,Sequence,flag,Y_score[0]))
                os.remove(Sequence)
        i1 = i1 +1
fout.close()
