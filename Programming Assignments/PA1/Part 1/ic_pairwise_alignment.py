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
09/26/23:
    MOD:     Creation of file and initial organization
    AUTHOR:  Ian Chavez
    COMMENT:
    - Changed package from pairwise2 to Align based on error:
    "BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated, and we intend to remove it in a future release of Biopython. As an alternative, please consider using Bio.Align.PairwiseAligner as a replacement, and contact the Biopython developers if you still need the Bio.pairwise2 module."

    - Addressed following overflow error by using documentation found at: http://biopython.org/DIST/docs/tutorial/Tutorial.html
    "OverflowError: number of optimal alignments is larger than 9223372036854775807"

    - Getting ValueError from running
    "ValueError: sequence contains letters not in the alphabet"
    
10/17/23:
    MOD:     Part 1 done
    AUTHOR:  Ian Chavez
    COMMENT:
    - Changed back to pairwise2 to keep consistent with assignment. 
    Not trying to lose points because of a warning.

    - Addressed multiple errors and warnings by using documentation found at: http://biopython.org/DIST/docs/tutorial/Tutorial.html

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
for a in pairwise2.align.globaldx(seq1, seq2, matrix):
    print(pairwise2.format_alignment(*a))
print("Print global successful====================")


# Local alignment using BLOSUM62 matrix
# and print the results to the console (using format_alignment function)
print("Local alignment starting====================")
alignments = pairwise2.align.localdx(seq1, seq2, matrix)
for a in alignments:
    print(pairwise2.format_alignment(*a))
print("Local alignment successful====================")

# Close the files
seq1_file.close()
seq2_file.close()
