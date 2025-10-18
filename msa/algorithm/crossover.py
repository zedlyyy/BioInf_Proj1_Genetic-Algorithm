import sys
sys.path.append("./")

from msa.multiplealign.myalign import MyAlign
from msa.algorithm.scoreSumPairs import score_msa_sum_of_pairs

SUBS_MATRIX_PATH    = "msa/multiplealign/blosum62.mat"

def align_Msa (sequence):
    return MyAlign(sequence, "protein")
    
def count_Residues(sequence):
    return len(sequence) - sequence.count("-")

def index_At_Residues(sequence, nr_Residues):
    residues_Found = 0
    index = 0
    
    while residues_Found != nr_Residues:
        if sequence[index] == "-":
            index += 1
        else: 
            residues_Found += 1
            index += 1

    return index 

def residues_to_indexes (sequences, residues):
    indexes = []

    for seq,residue in zip(sequences, residues):
        indexes.append(index_At_Residues(seq, residue))

    return indexes

def split_at (sequences, residues_row):
    
    #Se for só um inteiro, fazer uma lista com esse inteiro de size sequences.size
    if isinstance(residues_row, int):
        residues_row = [residues_row] * len(sequences)

    left_Side   = []
    right_side  = []
    residues    = []

    #Em cada seq, dar split no index indicado
    index = 0
    residues_left = 0
    for seq in sequences:
        left_Side.append(seq[:residues_row[index]])
        right_side.append(seq[residues_row[index]:])
        residues_left = count_Residues(left_Side[index])
        residues.append(residues_left)
        index += 1

    return left_Side, right_side, residues

def pad_aligment(sequences, side):

    max_length = 0
    padded_sequences = [] 

    for seq in sequences:
        seq_length = len(seq)
        if max_length < seq_length:
            max_length = seq_length

    #Lado direito e esquerdo da sequência
    for seq in sequences:
        seq_length = len(seq)
        match side:
            case "r":
                seq = seq + ("-") * (max_length - seq_length)
            case "l":
                seq = ("-") * (max_length - seq_length) + seq
            case _:
                print("Wrong side")
        padded_sequences.append(seq)

    return padded_sequences

def merge(sequences1, sequences2):

    seq_merged = []
    for seq1,seq2 in zip(sequences1, sequences2):
        seq_merged.append(seq1 + seq2)
    return seq_merged

def is_all_gaps(sequences, column):
    is_gap = True
    for seq in sequences:
        if seq[column] != "-":
            is_gap = False
    return is_gap

def clean_msa(sequences):
    clean_msa = sequences
    msa_length = len(clean_msa[0])

    column = 0
    #Para cada coluna, começamos do fim para o inicio devido à remoção de colunas
    for column in range(msa_length - 1, -1, -1):
        if is_all_gaps(clean_msa, column):
            #Para cada sequência dentro do clean_msa
            for i_sequence in range(len(clean_msa)):
                seq = clean_msa[i_sequence]
                clean_msa[i_sequence] = seq[:column] + seq[column+1:]

    return clean_msa

def generate_offspring(sequence1, sequence2, collum):

    align1 = align_Msa(sequence1)
    align2 = align_Msa(sequence2)

    l1, r1, residues = split_at(align1, collum)
    seq2_indexes = residues_to_indexes(align2, residues)

    l2, r2, _ = split_at(align2, seq2_indexes)

    l2_pad = pad_aligment(l2,"r")
    r2_pad = pad_aligment(r2,"l")

    l1r2_padded = merge(l1,r2_pad)
    l2r1_padded = merge(l2_pad,r1)

    l1r2_padded_clean = clean_msa(l1r2_padded)
    l2r1_padded_clean = clean_msa(l2r1_padded)

    score_child1 = score_msa_sum_of_pairs(l1r2_padded_clean,SUBS_MATRIX_PATH)
    score_child2 = score_msa_sum_of_pairs(l2r1_padded_clean,SUBS_MATRIX_PATH)

    #Dá return na melhor child
    if score_child1 > score_child2:
        return l1r2_padded_clean
    else :
        return l2r1_padded_clean

