# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 22:51:42 2015

@author: mark klick
@MIDTERM EXAM BMI550
"""
import regex as re
def unique_short_subs(s):
    for j in range(len(s)):
        x = []
        for i in range(len(s)):
            #if len(re.findall(s[i:i+j],s)) == 1:#regex
            if s.find(s[i:i+j]) == s.rfind(s[i:i+j]):#first index matches last one i.e. unique.
                #print i,j,s[i:i+j]
                x.append(s[i:i+j])
        if x != []:
            return x
#print unique_short_subs(b)

#BLOSUM62 file retrieved from ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62
#read in blosum matrix into a dictionary of dictionaries to represent the rows and columns of BLOSUM62 for easy lookup give AA x AA
with open('blosum62.txt') as f1:
    AAs = f1.readline().rstrip().split()#grab the AA column names
    blos_mat = {}#tried to use dict comprehension but couldn't get it to work out :/
    for line in f1:
        row = line.split()
        #grab AA from row
        AA = row.pop(0)#good ol c++ vector method forgot it worked on .py lists but perfect for this dict creation e.g.pop head off vector
        blos_mat[AA] = {}
        #dict with each column name as the key and value = the BLOSUM62 substitution score for that [row][col] AA x AA sub score
        for col in AAs:
            blos_mat[AA][col] = row.pop(0) 

def NW_global_align(seq1,seq2):
    #gap_pen = -1
    traceback = [[0 for i in range(len(seq1)+1)] for i in range(len(seq2)+1)]
    for i in range(len(seq1)+1):
        traceback[i][0] = int(blos_mat[seq1[0]]['*'])*i
    for i in range(len(seq2)+1):
        traceback[0][i] = int(blos_mat['*'][seq1[0]])*i
    
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            diag = int(traceback[i-1][j-1]) + int(blos_mat[seq1[i-1]][seq2[j-1]])
            up = int(traceback[i-1][j]) + int(blos_mat[seq1[i-1]]['*'])
            #up = int(traceback[i-1][j]) + gap_pen
            left = int(traceback[i][j-1]) + int(blos_mat['*'][seq2[j-1]])
            #left = int(traceback[i][j-1]) + gap_pen
            traceback[i][j] = max(diag,up,left)
            #print traceback[i][j]
    with open('NW_traceback_matrix1.txt','w') as f2:
        for i in range(len(seq1)+1):
            f2.write('\n')
            for j in range(len(seq2)+1):
                f2.write(str(traceback[i][j])+'\t')
    align1 = ''
    align2 = ''
    i = len(seq1)
    j = len(seq2)
    while(i>0 or j>0):
        if i>0 and j>0 and traceback[i][j]==int(traceback[i-1][j-1]) + int(blos_mat[seq1[i-1]][seq2[j-1]]):
            align1 += seq1[i-1] 
            align2 += seq2[j-1] 
            i = i-1
            j = j-1
        elif i>0 and traceback[i][j] ==int(traceback[i-1][j]) + int(blos_mat[seq1[i-1]]['*']):
        #elif i>0 and traceback[i][j] ==int(traceback[i-1][j]) + gap_pen:
            align1 += seq1[i-1]
            align2 += '-' 
            i = i-1
        else:
            align1 += '-' 
            align2 += seq2[j-1] 
            j = j-1
    align1 = align1[::-1]
    align2 = align2[::-1]
    with open ('NW_global_alignment1.txt','w') as f3:
        f3.write(align1 + '\n' + align2)

def SW_local_align(seq1,seq2):
    traceback = [[1 for i in range(len(seq1)+1)] for i in range(len(seq2)+1)]
    #path = [[1 for i in range(len(seq1)+1)] for i in range(len(seq2)+1)]
    for i in range(len(seq1)+1):
        traceback[i][0] = 0
    for i in range(len(seq2)+1):
        traceback[0][i] = 0
    maxx = 0
    maxx_cood = (0,0)
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            diag = int(traceback[i-1][j-1]) + int(blos_mat[seq1[i-1]][seq2[j-1]])
            up = int(traceback[i-1][j]) + int(blos_mat[seq1[i-1]]['*'])
            #up = int(traceback[i-1][j]) + gap_pen
            left = int(traceback[i][j-1]) + int(blos_mat['*'][seq2[j-1]])
            #left = int(traceback[i][j-1]) + gap_pen
            traceback[i][j] = max(0,diag,up,left)
            if traceback[i][j] >= maxx:
                maxx = traceback[i][j]
                maxx_cood = (i,j)
    with open('sW_traceback_matrix1.txt','w') as f2:
        for i in range(len(seq1)+1):
            f2.write('\n')
            for j in range(len(seq2)+1):
                f2.write(str(traceback[i][j])+'\t')    
    align1 = ''
    align2 = ''
    i = maxx_cood[0]
    j = maxx_cood[1]
    while(i>0 or j>0):
        if traceback[i][j]==0:
            break
        if i>0 and j>0 and traceback[i][j]==int(traceback[i-1][j-1]) + int(blos_mat[seq1[i-1]][seq2[j-1]]):
            align1 += seq1[i-1] 
            align2 += seq2[j-1] 
            i = i-1
            j = j-1
        elif i>0 and traceback[i][j] ==int(traceback[i-1][j]) + int(blos_mat[seq1[i-1]]['*']):
        #elif i>0 and traceback[i][j] ==int(traceback[i-1][j]) + gap_pen:
            align1 += seq1[i-1]
            align2 += '-' 
            i = i-1
        else:
            align1 += '-' 
            align2 += seq2[j-1] 
            j = j-1
    align1 = align1[::-1]
    align2 = align2[::-1]
    with open ('SW_local_alignment.txt','w') as f5:
        f5.write(align1 + '\n' + align2)
        
#seq1 = 'GCATWYACLA' 
#seq2 = 'AVKHYWTAGL'
#NW_global_align(seq1,seq2)   
#SW_local_align(seq1,seq2) 

import random
AAs = ['W','C','H','Y','P','F','D','N','G','M','K','R','Q','E','T','V','L','I','A','S']

def identity(seq1,seq2):
    identity = 0
    for i in range(len(seq1)):
        if seq1[i]==seq2[i]:
            identity = identity + 1
    identity = float(identity)/float(len(seq1))
    return identity
    
def disimilar(seq1,seq2):
    blos = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            if int(blos_mat[seq1[i]][seq2[i]]) <= -2:
                blos = blos+1
    return blos
    
def similar(seq1,seq2):
    blos = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            if int(blos_mat[seq1[i]][seq2[i]]) >= 3:
                blos = blos+1
    return blos  
    
identity_threshold = 0.3
dissim_threshold = 4 #at least 4 residues are very dissimilar 
sim_threshold = 4 #at least 4 residues are very similar
#intialize our peptides
pep1 =[AAs[random.randint(0,19)] for i in range(10) ]
pep2 =[AAs[random.randint(0,19)] for i in range(10) ]
pep3 =[AAs[random.randint(0,19)] for i in range(10) ]
pep4 =[AAs[random.randint(0,19)] for i in range(10) ]
#filter to meet similarity criteria
#while identity(pep1,pep2) <=  identity_threshold:
#    pep1 =[AAs[random.randint(0,19)] for i in range(10) ]
#    pep2 =[AAs[random.randint(0,19)] for i in range(10) ]
while disimilar(pep1,pep2) <=  dissim_threshold:
    pep1 =[AAs[random.randint(0,19)] for i in range(10) ]
    pep2 =[AAs[random.randint(0,19)] for i in range(10) ]
while similar(pep3,pep4) <= sim_threshold:
    pep3 =[AAs[random.randint(0,19)] for i in range(10) ]
    pep4 =[AAs[random.randint(0,19)] for i in range(10) ]

dis_seq1 = ''.join(pep1)
dis_seq2 = ''.join(pep2)
sim_seq1 = ''.join(pep3)
sim_seq2 = ''.join(pep4)

NW_global_align(dis_seq1,dis_seq2)
SW_local_align(dis_seq1,dis_seq2)

NW_global_align(sim_seq1,sim_seq2)
SW_local_align(sim_seq1,sim_seq2)