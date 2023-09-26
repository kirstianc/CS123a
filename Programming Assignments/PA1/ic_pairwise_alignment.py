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
    
====================== END OF MODIFICATION HISTORY ============================
"""
# Based on how http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html   
#   describes pairwise local and global alignment using Biopython
from Bio import Align
from Bio.Align import substitution_matrices
aligner = Align.PairwiseAligner()

# Read in the sequences from the files: seq1.txt and seq2.txt
#   and store them in variables
seq1 = open("seq1.txt", "r")
seq2 = open("seq2.txt", "r")
seq1 = seq1.read()
seq2 = seq2.read()

#TEST Print the sequences to the console
#print("Sequence 1: ", seq1)
#print("Sequence 2: ", seq2)

# Global alignment using BLOSUM62 matrix
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
alignment = aligner.align(seq1, seq2)
for alignment in sorted(alignment):
    print("Score = %.1f:" % alignment.score)
    print(alignment)

# Local alignment using BLOSUM62 matrix
aligner.mode = 'local'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
alignment = aligner.align(seq1, seq2)
for alignment in sorted(alignment):
    print("Score = %.1f:" % alignment.score)
    print(alignment)


""" BASED ON PAIRWISE2 DOCUMENTATION BUT DEPRECATED
# Global alignment using BLOSUM62 matrix
#   and print the results to the console (using format_alignment function)
alignments = pairwise2.align.globaldx(seq1, seq2, matrix)
for a in alignments:
    print(format_alignment(*a))

# Local alignment using BLOSUM62 matrix
#   and print the results to the console (using format_alignment function)
alignments = pairwise2.align.localdx(seq1, seq2, matrix)
for a in alignments:
    print(format_alignment(*a))
"""

# Close the files
seq1.close()
seq2.close()
