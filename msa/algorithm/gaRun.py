import sys
import random
sys.path.append("./")

from msa.algorithm.crossover import clean_msa, generate_offspring
from msa.algorithm.scoreSumPairs import score_msa_sum_of_pairs

#Tamanho minimo e máximo que uma sequence pode ter a mais que o tamanho da maior sequencia existente
MIN_ADDED_SIZE  = 3
MAX_ADDED_SIZE  = 10

COLUMN_TO_SPLIT_ON_BREEDING = 3


def run_ga (seqs, pop_amount, subs_matrix_path, generation_amount):
    msa_list = initialize_population(seqs, pop_amount)
    joined_score_msa = []

    for _ in range(generation_amount):
        
        #Resetar a list de scores
        score = []
        for msa in msa_list:
            score.append(score_msa_sum_of_pairs(msa, subs_matrix_path))
        
        #Juntar cada msa com o seu score
        joined_score_msa = list(zip(msa_list,score))

        #Dar sort de acordo com o objeto na segunda posição (score)
        # Reverse a true para colocar os scores por ordem descendente
        # os objetos na primeira posição são respetivamente a melhor msa, e o melhor score
        joined_score_msa.sort(reverse=True, key= lambda x:x[1])
    
        #Preparar a próxima geração
        msa_list = nextGeneration(joined_score_msa)

    best_msa, best_score = joined_score_msa[0]
    return best_msa, best_score


#Verificar o comprimento da maior sequencia e adicionar um numero random de gaps (3-10)
#   Ir a cada sequencia e adicionar gaps no inicio até que tenham o tamanho max_sequence_lenght + random(3-10)
def initialize_population(sequences, size_population):
    msa_length              = len(sequences)
    max_sequence_length     = 0

    msas_list = []

    #Verificar o size da maior sequencia do file
    for seq in sequences:
        if len(seq) > max_sequence_length:
            max_sequence_length = len(seq)
    
    #Size universal para dar diversidade às MSA's
    final_sequence_length = max_sequence_length + random.randint(MIN_ADDED_SIZE,MAX_ADDED_SIZE)

    #Cada MSA
    for pop in range(size_population):
        msa_list = []
        #Dentro de cada MSA
        for current_sequence in range(msa_length):
            seq_size = len(sequences[current_sequence])

            #Acresentar à esquerda da sequencia um número random entre 0 e o máximo que 
            # uma sequencia pode ter de tamanho menos o actual size da sequencia.
            left_gaps   = random.randint(0,final_sequence_length - seq_size)
            #Acrescentar à direita da sequencia o que restar do espaço que não foi utilizado
            # pelos gaps da esquerda, matendo assim o mesmo tamanho de sequencia dentro de cada MSA 
            right_gaps  = final_sequence_length - seq_size - left_gaps

            sec_padded = "-" * left_gaps + sequences[current_sequence] + "-" * right_gaps
            msa_list.append(sec_padded)

        #Antes de dar append, verificar as colunas que apenas contêm gaps
        msa_cleaned = clean_msa(msa_list)
        msas_list.append(msa_cleaned)
    
    return msas_list


def nextGeneration(joined_score_msa):
    msa_list, score_list = zip(*joined_score_msa)
    msa_list    = list(msa_list)
    score_list  = list(score_list)
    half = len(msa_list)//2
    other_half = len(msa_list) - half
    new_generation = []

    elite_msa   = msa_list[:half]
    elite_msa_list = list(elite_msa)

    for seq in elite_msa_list:
        new_generation.append(seq)

    #breeding_msa    = msa_list[half:]
    #breeding_score  = score_list[half:]


    #Normalizar todos os scores para que não hajam valores negativos
    min_score = min(score_list)
    for index, score in enumerate(score_list):
        if min_score < 0:
            score_list[index] = score - min_score

    #Criar uma list com prob de cada MSA ser escolhido baseado no score
    sum_score = sum(score_list)
    prob = []
    if sum_score == 0:
        prob = [1 / len(score_list)] * len(score_list)
    else :
        for score in score_list:
            prob.append(score / sum_score)

    #Criar uma lista com as probabilidades comulativas para facilitar a escolha do MSA 
    comulative_prob = []
    p_sum = 0
    for p in prob:
        p_sum += p
        comulative_prob.append(p_sum)

    #Escolher os respetivos pais/ Escolher a mutação
    for index in range(other_half, len(msa_list)):
        r = random.random()
        if r > 0.5 :
            new_child = getChild(comulative_prob,msa_list)
            new_generation.append(new_child)
        else:
            new_protein_mutation = mutation(msa_list[index])
            new_generation.append(new_protein_mutation)
    
    return new_generation

def getChild(comulative_prob, breeding_msa):
    
    parents = []
    attempts = 0
    MAX_ATTEMPTS = 1000
    #Até dois parents
    while len(parents) < 2 and attempts < MAX_ATTEMPTS:
        attempts += 1
        r = random.random()
        for index, prob in enumerate(comulative_prob):
            if prob > r:
                possible_parent = breeding_msa[index]
                #Só adicionar se o parent escolhido não estiver nos parents escolhidos
                if possible_parent not in parents:
                    parents.append(possible_parent)
                break
    #Infinite loop problem, se a lista só conter sequencias MSA iguais e/ou se a chance de um MSA for bastante grande
    #  então escolhemos dois pais ao calhas
    if len(parents) < 2:
        parents = random.sample(breeding_msa,2)


    parent1, parent2 = parents
    return generate_offspring(parent1, parent2, COLUMN_TO_SPLIT_ON_BREEDING)

#Pegar numa sequencia dentro do msa e eliminar todos os gaps,
# depois colocar o número de gaps eliminados no inicio da sequência
def mutation(msa):
    mutated_msa = msa
    r = random.randint(0, len(mutated_msa) - 1)
    sequence_to_mutate = mutated_msa.pop(r)

    sequence_no_gaps = sequence_to_mutate.replace("-","")
    sequence_mutated = "-" * (len(sequence_to_mutate) - len(sequence_no_gaps)) + sequence_no_gaps

    mutated_msa.append(sequence_mutated)
    
    return mutated_msa



