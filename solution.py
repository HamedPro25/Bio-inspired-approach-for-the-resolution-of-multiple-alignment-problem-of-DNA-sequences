#Liste des acide amines
acid_amines=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
import random
class Solution:
    score=0
    sequence=[]
    def __init__(self,seq_array):
        self.sequence=seq_array
    
    #methode de calculer le score entre deux sequences
    def score_pairwise(self,seq1, seq2, gap_s, gap_e, gap = True):
        sum=0
        for A,B in zip(seq1, seq2):
            diag = ('-'==A) or ('-'==B)
            if diag:
                if gap:
                    sum+=gap_s
                else:
                    sum+=gap_e
            elif A!=B:
                sum+=1
        return sum
    #methode calculer score d'un solution
    def calculate_score(self):
        sc=0
        for i in range(1,len(self.sequence)):
            for j in range(0,i-1):
                sc+=self.score_pairwise(self.sequence[i], self.sequence[j], +3,+2)
        self.score=sc
        #print(self.score)
    def get_sequencess(self):
        return self.sequence
    #croissement entre deux sequence
    def signle_cross(self,seq1,seq2,ran,index):
        s1part1=seq1[0:ran]
        s1part2=seq1[ran:]
        s2part1=seq2[0:ran]
        s2part2=seq2[ran:]
        seq1=str(s1part1)+str(s2part2)
        seq2=str(s2part1)+str(s1part2)
        return (seq1,seq2)
    #croissmsnt avec autre solution
    def crossover(self,ran,seq):
        result1=[]
        result2=[]
        for i in range(len(self.sequence)):
            result1.append(self.signle_cross(self.sequence[i],seq[i],ran,i)[0])
            result2.append(self.signle_cross(self.sequence[i],seq[i],ran,i)[1])
        return (result1,result2) 
    
    #methode de mutation
    def mutation(self):
        seqnb=random.randint(0,len(self.sequence)-1)
        s=list(self.sequence[seqnb])
        position=random.randint(0,len(self.sequence[seqnb])-1)
        if position!=len(self.sequence[seqnb])-1:
            s[position],s[position+1]=s[position+1],s[position]
        else:
            s[position],s[0]=s[0],s[position]
        self.sequence[seqnb]="".join(s)
    

    #methode creer voisinage d'un solution
    def voisinage(self):
        seqnb=random.randint(0,len(self.sequence)-1)
        s=list(self.sequence[seqnb])
        position=random.randint(0,len(self.sequence[seqnb])-1)
        random_caracter = random.choice(acid_amines)
        s[position]=str(random_caracter)
        self.sequence[seqnb]="".join(s)
