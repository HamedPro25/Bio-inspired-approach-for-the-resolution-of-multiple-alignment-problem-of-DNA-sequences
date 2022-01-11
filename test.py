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


formatters = {
    'int': lambda x: '%4d' % x,
    'float': lambda x: '%.02f' % x
}


def print_report(population, fitness, wheel, selectors, selected_individuals):#taffichie le resulta ta3 methode de roulete
    with np.printoptions(formatter=formatters):
        print('          Population:', population)
        print('             Fitness:', fitness)
        print('      Roulette wheel:', wheel)       # fitness cumulative sum
        print('      Sampled values:', selectors)   # roulette "extractions"
        print('Selected individuals:', selected_individuals)


# This should be equivalent to np.choice(population, size, weights=fitness)
#methode ta3 selection de roulette
#param2: list des individus,param3: list ta3 les scores ta3 kol individu param4: nombre des individus
def roulette_wheel_selection(rng: np.random.Generator, 
                             population: np.ndarray,
                             fitness: np.ndarray,
                             size: int,indice) -> np.ndarray:
    if size > len(population):
        raise ValueError
    fitness_cumsum = fitness.cumsum()   # the "roulette wheel"
    
    fitness_sum = fitness_cumsum[-1]    # sum of all fitness values (size of the wheel)
    sampled_values = rng.random(size) * fitness_sum
    
    # For each sampled value, get the corresponding roulette wheel slot
    selected = np.searchsorted(fitness_cumsum, sampled_values)
    #print_report(population, fitness, fitness_cumsum, sampled_values, selected)
    return selected


#methode Alignment en utilsant l'algorithme genetique
def AG(Nb_Iter,original_sequences,crois,mut,pop):
    begin=time.clock()
    max_seq=max(original_sequences,key=len)
    max_lent=len(max_seq)+2
    solutionss=pop
    bestalign=[]
    bestscore=[]
    scores_init=[]
    #evaluation de primiere population
    for sol in solutionss:
            sol.calculate_score()
            scores_init.append(sol.score) #remplir list des scores
    bestscore.append(min(scores_init))
    index=scores_init.index(bestscore[-1])
    bestalign.append(solutionss[index].sequence)
    #Parcourir les iterations
    for j in range(Nb_Iter):
        #list des scores
        scores=[]
        solutions=copy.copy(solutionss)
        i=0
        #Evaluation de generation
        for sol in solutions:
            sol.calculate_score()
            scores.append(sol.score) #remplir list des scores
            i+=1
        muttrv=False
        croistrv=False
        mini=min(scores)
        indice=scores.index(mini)
        select_numbers=int(len(solutions)/2)#numer de selection des individus 50%
        rng=default_rng()
        population_nb=np.arange(len(solutions))
        fitness=np.array(scores)
        #selection des individus en utilsant la methode de roulette
        select_individulas=roulette_wheel_selection(rng, population_nb, fitness, select_numbers,indice)#selection par roullette
        
        while indice in select_individulas:
            select_individulas=roulette_wheel_selection(rng, population_nb, fitness, select_numbers,indice)#selection par roullette
        operator_choice=random.random()
        #Condioton de faire l'opearation de mutation
        if operator_choice<mut:
            muttrv=True
            print("mutation operator")
            for i in range(len(select_individulas)):
                solutions[select_individulas[i]].mutation()
                solutions[select_individulas[i]].calculate_score()
        #Condioton de faire l'opearation de croissment
        if operator_choice<crois:
            croistrv=True
            print("croosover operator")
            for i in range(0,len(select_individulas),2):
                position=random.randint(2,max_lent)
                if i==len(select_individulas)-1:
                    results=(solutions[select_individulas[i]].crossover(position,solutions[select_individulas[0]].sequence))
                    solutions[select_individulas[i]].sequence=results[0]
                    solutions[select_individulas[0]].sequence=results[1]
                    solutions[select_individulas[i]].calculate_score()
                    solutions[select_individulas[0]].calculate_score()    
                else:
                    results=(solutions[select_individulas[i]].crossover(position,solutions[select_individulas[i+1]].sequence))
                    solutions[select_individulas[i]].sequence=results[0]
                    solutions[select_individulas[i+1]].sequence=results[1]
                    solutions[select_individulas[i]].calculate_score()
                    solutions[select_individulas[i+1]].calculate_score()
        #Condition de pass au gen suivante
        if not muttrv and not croistrv:
            print("pass to next generation")
        scores=[]
        for sol in solutions:
            scores.append(sol.score) #remplir list des scores
        #Ajouter la meuilleur solution
        bestscore.append(min(scores))
        index=scores.index(bestscore[-1])
        bestalign.append(solutions[index].sequence)
    iterat=[]
    for i in range(Nb_Iter+1):
        iterat.append(i)
    minimu=max(bestscore)
    ind=0
    for i in range(len(bestscore)):
        if bestscore[i]<=minimu:
            minimu=bestscore[i]
            ind=i
    meuilleur=bestalign[ind]
    return (bestscore,iterat,time.clock()-begin,meuilleur)

