import sys
sys.path.append("./")

from msa.algorithm.gaRun import run_ga 
from Bio import SeqIO

FILE_PATH           = "data/cytochromes.fa"
SUBS_MATRIX_PATH    = "msa/multiplealign/blosum62.mat"
POP_SIZE          = 20
GENERATION_SIZE   = 50

def main():
    seqs = read_fasta()
    best_msa, best_score = run_ga(seqs, POP_SIZE, SUBS_MATRIX_PATH, GENERATION_SIZE)
    print(f"Best MSA: {best_msa} \n Best score : {best_score}")

def read_fasta():
    cytochromes = list(SeqIO.parse((FILE_PATH), "fasta-blast")) #Alteração no fasta para compatibilidade
    cytochromes_seq = []
    for seq in cytochromes:
        cytochromes_seq.append(str(seq.seq))
    return cytochromes_seq

if __name__ == "__main__":
    main()