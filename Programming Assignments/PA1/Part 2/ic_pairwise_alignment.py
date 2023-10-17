# coding: utf-8
#NAME:  ic_pairwise_alignment.py
#DESCRIPTION: This python script will read in the
# protein sequences of the files: seq1.txt and seq2.txt  
# and both globally and locally aligns them using the BLOSUM62 matrix

"""
AUTHOR: Ian Chavez

   Unpublished-rights reserved under the copyright laws of the United States.

   This data and information is proprietary to, and a valuable trade secret
   of, Leonard P. Wesley and Ian Chavez. It is given in confidence by Leonard
   P. Wesley and Ian Chavez. Its use, duplication, or disclosure is subject to
   the restrictions set forth in the License Agreement under which it has been
   distributed.

      Unpublished Copyright © 2023 Leonard P. Wesley and Ian Chavez
      All Rights Reserved

========================== MODIFICATION HISTORY ==============================
10/17/23:
    MOD:     Creation of file for part 2
    AUTHOR:  Ian Chavez
    COMMENT:
    - Copy file from Part 1, removed local alignment (unnecessary part 2)

    - Experimenting with different parameters to improve global score

====================== END OF MODIFICATION HISTORY ============================
"""
# Read in the sequences from the files: seq1.txt and seq2.txt
# and store them in variables
with open("seq1.txt", "r") as seq1_file:
    lines = seq1_file.readlines()
    seq1 = lines[1].strip()

with open("seq2.txt", "r") as seq2_file:
    lines = seq2_file.readlines()
    seq2 = lines[1].strip()

# Import the pairwise2 module from Biopython
from Bio import pairwise2
from Bio.Align import substitution_matrices
print("Imports successful")

# Create a BLOSUM62 matrix
matrix = substitution_matrices.load("BLOSUM62")
print("Matrix created")

# Global alignment using BLOSUM62 matrix
# and print the results to the console (using format_alignment function)
print("Global alignment starting====================")
for a in pairwise2.align.globalmx(seq1, seq2, 6, -1):
    print(pairwise2.format_alignment(*a))
print("Print global successful====================")

# Close the files
seq1_file.close()
seq2_file.close()
