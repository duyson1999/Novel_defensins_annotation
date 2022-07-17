####### 
# PURPOSE OF THIS FILE: 
# 1. Because there are huge number of hits (result) after running tblastn -> need to narrow it 
# 2. Transfer all output into new txt file under tabular format
#
#
#######################################     # From New_CDS.txt, retrieve all sequences into one fasta file
import Bio
from Bio import Entrez
import pandas as pd
from pandas import DataFrame
#path = input('link to the file: ')
#path_to_newCDS = path
stitle = []
location = []
path_to_newCDS = '/Users/sonvo/Documents/1.All_for_Genetic_Project/ForThesis/defensin_BSF/New_CDS.txt'
newCDS = pd.read_csv(path_to_newCDS, sep = '\t')


chromosome_lst = newCDS['chromosome'].unique().tolist()
#for row in newCDS.loc[:,'chromosome']:
#    if row not in chromosome_lst:
#        chromosome_lst.append(row)
accession_lst = newCDS['stitle'].unique().tolist()
#for row2 in newCDS.loc[:,'stitle']:
#    if row2 not in accession_lst:
#        accession_lst.append(row2)
        

#######################################     # Getting hits in chromosome1 from result_file
def retrieve_data_from_one_chromosome (number):
    chromosome = newCDS[(newCDS['chromosome'] == number)]
    locationlist_chromosome = []
    indexlst = chromosome.index
    for K in indexlst:
        A = chromosome.loc[K,'sstart']
        B = chromosome.loc[K,'send']
        if A <= B:
            locationlist_chromosome.append([A-1000,B+1000])
        else:
            locationlist_chromosome.append([B-1000,A+1000])
    locationlist_chromosome.sort(key = lambda x: x[0])
    return locationlist_chromosome


#######################################      # merge lists with same starting location (supplement)
def merge_list(list_with_same_start):      
    sublst = [X[1] for X in list_with_same_start]
    return [list_with_same_start[0][0], max(sublst)]

#######################################     # Join location if they are belong to the other (supplement)
def overlap(lst1,lst2):
    if lst1[0] <= lst2[0] <= lst1[1] or lst1[0] <= lst2[1] <= lst1[1] or lst2[0] <= lst1[0] <= lst1[1] <= lst2[1]:
        return True
    else:
        return False
def concatenate_list(biglst):     
    x = 0
    for i in range (len(biglst)-1):
        if overlap(biglst[i], biglst[i+1]):
            x += 1
            sublst= [min(biglst[i][0], biglst[i+1][0]), max (biglst[i][1], biglst[i+1][1])]
            anotherlst = [sublst] + [X for X in biglst if X !=biglst[i] if X != biglst[i+1]]
            break
    if x > 0:
        anotherlst.sort(key= lambda z: z[0])
        return concatenate_list(anotherlst)
    else:
        return biglst

#########################################    # Main code
def edit(lst):
    final_result = []
    result = []
    for loc in lst:         # remove duplicates                 
        if loc not in result:
            result.append(loc)      
    subset1 = []            # Merge location
    for loc2 in result:
        if loc2[0] not in subset1:
            subset1.append(loc2[0])
    for element in subset1:
        sublst = [X for X in result if X[0] == element]
        final_result.append(merge_list(sublst))
    return concatenate_list(final_result)

#########################################          # Create new dataframe to export into new file
sstart = []
send = []
chromosome  = []
accession = []
for number in chromosome_lst:
    final_concatenate_list = edit(retrieve_data_from_one_chromosome(number))
    sstart += [X[0] for X in final_concatenate_list]
    send += [X[1] for X in final_concatenate_list]
    chromosome += [number for i in range (len(final_concatenate_list))]
    accession += [accession_lst[chromosome_lst.index(number)] for j in range (len(final_concatenate_list))]

newdf = DataFrame({'sstart': sstart, 'send':send, 'chromosome':chromosome, 'stitle' : accession })
#print(newdf)
#print(newdf)
newdf.to_csv('Join_hits.tsv', sep='\t', index=False)

#########################################          # Retrieve sequences, store in 'fasta_seqs_newCDS.txt'
fasta_seqs = ""
Entrez.email = 'duyson1999@gmail.com'
for number in range(0,len(sstart)):
    record = Entrez.efetch(db = 'nucleotide', id = accession[number], rettype = 'fasta', retmode = 'text', seq_start = sstart[number], seq_stop = send[number])
    fasta_seqs += record.read()
    fasta_seqs += '\n'

open('fasta_seqs_newCDS.txt', 'w').write(fasta_seqs)



