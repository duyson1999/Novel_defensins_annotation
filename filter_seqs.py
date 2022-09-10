#############################################################################################
# PURPOSE OF THIS FILE:
# 1. Because there are huge number of hits (result) after running tblastn -> need to narrow it
# 2. Transfer all output into new txt file under tabular format
#############################################################################################

import pandas as pd
from Bio import Entrez

from constant import OUTPUT_DIR


def _gene_location_in_specific_chromosome(df, chromosome_number):
    chromosome = df[(df['chromosome'] == chromosome_number)]
    locationlist_chromosome = []
    for i in range(len(chromosome)):
        A = chromosome.iloc[i].loc['sstart']
        B = chromosome.iloc[i].loc['send']
        if A <= B:
            locationlist_chromosome.append([A-1000, B+1000])
        else:
            locationlist_chromosome.append([B-1000, A+1000])
    locationlist_chromosome.sort(key = lambda x: x[0])
    return locationlist_chromosome


def _merge_list(list_with_same_start):
    sublst = [X[1] for X in list_with_same_start]
    return [list_with_same_start[0][0], max(sublst)]


def _overlap(lst1, lst2):
	a, b = lst1
	c, d = lst2
	if a <= c <= b or a <= d <= b or c <= a <= b <= d:
		return True
	else:
		return False


def _concatenate_list(biglst):
    x = 0
    for i in range (len(biglst)-1):
        if _overlap(biglst[i], biglst[i+1]):
            x += 1
            sublst= [min(biglst[i][0], biglst[i+1][0]), max(biglst[i][1], biglst[i+1][1])]
            anotherlst = [sublst] + [y for y in biglst if y != biglst[i] if y != biglst[i+1]]
            break
    if x > 0:
        anotherlst.sort(key=lambda z: z[0])
        return _concatenate_list(anotherlst)
    else:
        return biglst


def _edit(lst):
    final_result = []
    result = []
    for loc in lst:
        if loc not in result:
            result.append(loc)
    subset1 = []
    for loc2 in result:
        if loc2[0] not in subset1:
            subset1.append(loc2[0])
    for element in subset1:
        sublst = [X for X in result if X[0] == element]
        final_result.append(_merge_list(sublst))
    return _concatenate_list(final_result)


def join_hit(CDS_df):
	chromosome_lst = CDS_df['chromosome'].unique().tolist()
	accession_lst = CDS_df['stitle'].unique().tolist()

	sstart = []
	send = []
	chromosome  = []
	accession = []
	for chr_number in chromosome_lst:
		location_lst = _gene_location_in_specific_chromosome(CDS_df, chr_number)
		print (1)
		final_concatenate_list = _edit(location_lst)
		print (2)
		sstart += [X[0] for X in final_concatenate_list]
		send += [X[1] for X in final_concatenate_list]
		chromosome += [chr_number for i in range (len(final_concatenate_list))]
		accession += [accession_lst[chromosome_lst.index(chr_number)] for j in range (len(final_concatenate_list))]
	return accession, sstart, send, chromosome


def retrieve_protein_seqs(accession, sstart, send):
	fasta_seqs = ''
	Entrez.email = 'duyson1999@gmail.com'
	for i in range(0,len(sstart)):
		record = Entrez.efetch(
			db = 'nucleotide', id = accession[i],
			rettype = 'fasta', retmode = 'text',
			seq_start = sstart[i], seq_stop = send[i]
		)
		fasta_seqs += record.read() + '\n'
	return fasta_seqs


if __name__ == '__main__':
	new_cds_df = pd.read_csv('{}/new_cds.txt'.format(OUTPUT_DIR), sep='\t', index_col=0)
	old_cds_df = pd.read_csv('{}/old_cds.txt'.format(OUTPUT_DIR), sep='\t', index_col=0)

	accession, sstart, send, chromosome = join_hit(new_cds_df)
	df = pd.DataFrame({
			'accession': accession,
			'sstart': sstart,
			'send': send,
			'chromosome': chromosome
		})
	df.to_csv('{}/new_cds_join_hits.tsv'.format(OUTPUT_DIR), sep='\t', index=False)
	fasta_seqs = retrieve_protein_seqs(accession, sstart, send)
	open('new_cds_seqs.txt', 'w').write(fasta_seqs)

	accession, sstart, send, chromosome = join_hit(old_cds_df)
	df = pd.DataFrame({
			'accession': accession,
			'sstart': sstart,
			'send': send,
			'chromosome': chromosome
		})
	df.to_csv('{}/old_cds_join_hits.tsv'.format(OUTPUT_DIR), sep='\t', index=False)
	fasta_seqs = retrieve_protein_seqs(accession, sstart, send)
	open('old_cds_seqs.txt', 'w').write(fasta_seqs)