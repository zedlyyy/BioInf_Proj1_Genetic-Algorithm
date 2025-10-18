import sys
sys.path.append("./")

from msa.algorithm.gaRun import mutation

align_exemple = [
    "ATC-G-G-TT",
    "A-CCGC-ATC",
    "A-CCG--ATC"
]

print(mutation(align_exemple))