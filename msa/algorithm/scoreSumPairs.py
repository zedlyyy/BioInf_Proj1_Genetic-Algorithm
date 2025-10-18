import sys
import random
sys.path.append("./")

from msa.multiplealign.substmatrix import SubstMatrix
from msa.multiplealign.myalign import MyAlign

GAP_PENALTY = -1

def score_msa_sum_of_pairs(msa, subs_matrix_path, gap_penalty = GAP_PENALTY):
    score = 0
    column_size = len(msa[0])
    msa_len = len(msa)
    msa_my_align = MyAlign(msa)

    #Dar load na BLOSUM62 submatrix 
    sub_matrix = SubstMatrix()
    sub_matrix.read_submat_file(subs_matrix_path)

    #Fazer tres loops, o primeiro a passar por cada coluna individual do MSA
    # o segundo a passar por todas as x sequências e o terceiro a passar pelas x+1 sequencias
    for column in range(column_size):

        column_residues = msa_my_align.column(column)

        for seq_line in range(msa_len - 1):
            for seq_line_next in range(seq_line+1, msa_len):
                if  column_residues[seq_line] == "-" or column_residues[seq_line_next] == "-":
                    score += gap_penalty
                else:
                    score += sub_matrix.sm[column_residues[seq_line] + column_residues[seq_line_next]]
                    
    return score