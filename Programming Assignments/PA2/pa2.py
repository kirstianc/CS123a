# coding: utf-8
#NAME:  pa2.py
#DESCRIPTION: This python script will obtain FASTA from NCBI and save them to a file. Then it will continue the Programming Assignment 2.

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
10/30/23:
    MOD:     Creation of file and initial organization
    AUTHOR:  Ian Chavez
    COMMENT:
        - Created file and added initial organization
        - Base functionaltiy
11/5/23:
    MOD:     Finalization of code and function + print statements
    AUTHOR:  Ian Chavez
    COMMENT:
        - Modification title^
====================== END OF MODIFICATION HISTORY ============================
"""
# Imports
from Bio import Entrez
import time

accession_numbers = []
species_gc_count_file = {}
species_gc_count_fasta = {}
fasta_ranking = {}
ranking_medians = {}

def obtain_accession_nums():
    with open('cs123a_list_of_acc_nums.txt', 'r') as file:
        for line in file:
            accession_numbers.append(line.strip())

def save_fasta(accession_numbers):
    Entrez.email = ""

    for accession_num in accession_numbers:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            with open("cs123a_fasta_files.txt", "a") as f:
                f.write(fasta_data + "\n")
            handle.close()
            time.sleep(1)
        except Exception as e:
            print(f"An error occurred: {e}")

def obtain_species_gc_count_from_file():
    # store species and gc count in a dictionary from cs123a_list_of_species_gc_count.txt
    with open('cs123a_list_of_species_gc_count.txt', 'r') as file:
        for line in file:
            data = line.split()
            key = ' '.join(data[:-2])
            
            # convert percentages to floats
            percentages = [float(percent.strip('%')) for percent in data[-2:]]
            species_gc_count_file[key] = percentages

def obtain_gc_count_from_fasta():
    # counts g and c and calculates gc count, stores in a dictionary
    with open('cs123a_fasta_files.txt', 'r') as file:
        for line in file:
            if line.startswith('>'):
                species = line.split()[1]
                species_gc_count_fasta[species] = [0,0,0]
            else:
                for char in line:
                    # count g and c
                    if char == 'G' or char == 'C':
                        species_gc_count_fasta[species][0] += 1
                    # count total number
                    species_gc_count_fasta[species][1] += 1

def compute_gc_content():
    # computes gc content and prints out the results
    for species in species_gc_count_fasta:
        gc_content = (species_gc_count_fasta[species][0] / species_gc_count_fasta[species][1]) * 100
        species_gc_count_fasta[species][2] = gc_content
        print(species + " " + str("{:.2f}%".format(gc_content)))

def ranking():
    # for species from fasta files
    curr_fasta = 0
    for fasta_species in species_gc_count_fasta:
        fasta_ranking[fasta_species] = ["N/A"]
        
        # for species given in assignment
        for file_species in species_gc_count_file:
            ranking_medians[file_species] = []
            
            if fasta_species in species_gc_count_fasta:
                # check if gc content is within range
                if species_gc_count_fasta[fasta_species][2] >= species_gc_count_file[file_species][0]:
                    if species_gc_count_fasta[fasta_species][2] <= species_gc_count_file[file_species][1]:
                        # within range, calculate median and compare
                        median = (species_gc_count_file[file_species][0] + species_gc_count_file[file_species][1]) / 2
                        ranking_medians[file_species] = abs(species_gc_count_fasta[fasta_species][2] - median)
                    else: 
                        # too large for range, -1 denotes this species to not be ranked
                        ranking_medians[file_species] = -1
                else:
                    # too small for range, -1 denotes this species to not be ranked
                    ranking_medians[file_species] = -1   
        
        # rank species, lower median difference is ranked higher than higher median difference
        for file_species in ranking_medians:
            if ranking_medians[file_species] != -1:
                if fasta_ranking[fasta_species][0] == "N/A":
                    fasta_ranking[fasta_species][0] = file_species
                # if median difference is smaller than current ranking, replace and push old ranking down
                elif ranking_medians[file_species] < ranking_medians[fasta_ranking[fasta_species][0]]:
                    temp = fasta_ranking[fasta_species][0]
                    fasta_ranking[fasta_species][0] = file_species
                    fasta_ranking[fasta_species].append(temp)
        
        # print out ranking for each species like as
        # ACC: <accession number>
        # RANK: <ranked species>
        print("ACC: " + accession_numbers[curr_fasta])
        curr_fasta += 1
        
        curr_rank = 1
        print("RANK: ")
        for rank in fasta_ranking[fasta_species]:
            print(str(curr_rank) + ". " + rank)
            curr_rank += 1
        print("\n")                

if __name__ == '__main__':
    
    print("Starting program...\n")
    
    print("Obtaining accession nums...")
    obtain_accession_nums()
    print(accession_numbers)
    print("\n")
    
    print("Saving FASTA...")
    save_fasta(accession_numbers)
    print("\n")
    
    print("Obtaining GC count from file...")
    obtain_species_gc_count_from_file()
    print(species_gc_count_file)
    print("\n")
    
    print("Obtaining GC count from FASTA...")    
    obtain_gc_count_from_fasta()
    print(species_gc_count_fasta)
    print("\n")

    print("Computing GC content of FASTA...")
    compute_gc_content()
    print("\n")

    print("Ranking organisms per entry...")
    ranking()