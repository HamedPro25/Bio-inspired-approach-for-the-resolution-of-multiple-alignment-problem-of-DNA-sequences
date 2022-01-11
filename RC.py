import xml.etree.ElementTree as ET
from Bio.pairwise2 import *
from Bio.Align.substitution_matrices import *
from Bio.SubsMat.MatrixInfo import *
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from numpy.random import default_rng
import numpy as np
import math
import solution
from solution import Solution
import time
#original_sequences=["WGKVNVDEVGGEAL","WDKVNEEEVGGEAL","WGKVGAHAGEYGAEAL","WSKVGGHAGEYGAEAL"]

import random
import copy
def seqence_align(deb,sequence):
    #fin_gap=""
    deb_gap=""
    if deb!=0:
        for i in range(deb):
            deb_gap+="-"

    return (str(deb_gap)+str(sequence)) #+str(fin_gap))  
def sequence_init(offset,index,original_sequences):
    sequence=original_sequences[index]
    while offset>0:
        position=random.randint(0,len(sequence))
        gap_nb=random.randint(0,offset)
        sequence_part1=sequence[0:position]
        sequence_part2=seqence_align(gap_nb,sequence[position:])
        sequence=str(sequence_part1)+str(sequence_part2)
        offset=offset-gap_nb
    #print(sequence)
    return sequence
        

def initialisation(po_nb,original_sequences):
    max_seq=max(original_sequences,key=len)
    max_lent=len(max_seq)+2
    solutions=[]
    offset_Array=[]
    for i in range(0,len(original_sequences)):
        offset_Array.append(max_lent-len(original_sequences[i]))
    for i in range(0,po_nb):
        sequence=[]
        for j in range(0,len(original_sequences)):
            sequence.append(sequence_init(offset_Array[j],j,original_sequences))
        #print("solution1"+str(sequence))
        solutions.append(Solution(sequence))
    return solutions



#methode de calcul l'aligment en utilsant le Recuit simule
def AG(seuil,t,meuill):
    begin=time.clock()
    bestscore=[]
    i=0
    iterat=[]
    Scour=Solution(meuill)
    Scour.calculate_score()
    while t>=seuil:

        Si=copy.deepcopy(Scour)
        #Creer nouvelle voisin de la solution courante
        Si.voisinage()  
        Si.calculate_score()     
        operator_choice=random.random()
        #Condition ou la nouvelle solution est mieux que la courante
        if Si.score<=Scour.score :
            print("voisinage operator")
            Scour.sequence=Si.sequence
            Scour.calculate_score()  
       #Cas ou la nouvelle solution est pire que la courante
        else:
            print("pass to next generation")
        #Ajouter meuilleur solution
        bestscore.append(Scour.score)
        iterat.append(t)
        t-=1
        i+=1
    return (bestscore,iterat,time.clock()-begin)