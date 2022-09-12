#############################################################################################
# PURPOSE OF THIS CODE:
# 1. Retrieve protein queries from source (ncbi, uniprot) using for tBLASTn
# 2. Perform tBLASTn
# 3. Seperate output file into newCDS and oldCDS - Based on 'gene' database in ncbi
#
#
# Note: the chromosome name format is different between species in ncbi.
#############################################################################################


import os
import re
import common

import pandas as pd
import multiprocessing as mp

from Bio import Entrez
from bioservices import UniProt

from constant import AMP, ROOT_DIR, WHOLE_GENOME_FASTA
from constant import OUTPUT_DIR
from constant import SPECIES_NAME
from constant import RELATIVE_SPECIES
from constant import UNIPROT_QUERY_LIMIT


def _format_species_query_ncbi(species_lst):
	lst = ['"{}"[Organism]'.format(x) for x in species_lst]
	return ' OR '.join(lst)


def _format_species_query_uniprot(species_lst):
	return ' OR '.join(species_lst)


def _parse_list2string(lst):
	res = ''
	for sublst in lst:
		res += sublst
	return res


def _get_chunk(max, n_per_chunk):
	chunks = [[x, x + n_per_chunk] for x in range(0, max-n_per_chunk, n_per_chunk)]
	chunks.extend([[chunks[-1][-1], max]])
	return chunks


def _uniprot_organisms_chunks(chunks, uniprot_organisms):
	return [' OR '.join(
		uniprot_organisms.split(' OR ')[x[0] : x[1]]
	) for x in chunks]


def _collect_reference_proteins_ncbi(
		amp_name,
		ncbi_organisms
	):
	Entrez.email = 'duyson1999@gmail.com'
	search_term = '{}[All Fields] AND ({})'.format(
		amp_name, ncbi_organisms
	)
	search = Entrez.esearch(
			db='protein', term=search_term,
			retmax=200, usehistory='y', idtype='acc'
		)
	search_result = Entrez.read(search)
	search.close()
	idlist = ','.join(search_result['IdList'])
	ncbi_data = Entrez.efetch(
			db='protein', id=idlist,
			rettype='fasta', retmode='text'
		).read()
	# print ('ncbi_data', ncbi_data)
	return ncbi_data


def _query_uniprot(params):
	amp_name = params['amp_name']
	uniprot_organisms = params['uniprot_organisms']
	u = UniProt()
	res = u.search(
		'{} ({})'.format(amp_name, uniprot_organisms),
		limit =200,  frmt = 'list'
	)
	uniprotdata = ''
	if not res:
		return uniprotdata
	for id in res.rstrip('\n').split('\n'):
		fastafile = u.retrieve(id, frmt='fasta')
		uniprotdata += fastafile
	# print ('uniprotdata', uniprotdata)
	return uniprotdata


def _collect_reference_proteins_uniprot(
		amp_name,
		uniprot_organisms
	):
	org_len = len(uniprot_organisms.split('OR'))
	chunks = _get_chunk(org_len, UNIPROT_QUERY_LIMIT)
	org_chunks = _uniprot_organisms_chunks(chunks, uniprot_organisms)
	with mp.Pool(5) as pool:
		results = pool.map_async(
			_query_uniprot, [{
				'amp_name': amp_name,
				'uniprot_organisms': org_chunk
			} for org_chunk in org_chunks]
		).get()
	return _parse_list2string(results)


def _restructure_fasta(fasta):
	res = ''
	for line in fasta.split('\n'):
		if line.startswith('>'):
			line = line.replace('>', '\n>')
			res += line + '\n'
		else:
			line = line.strip()
			res += line
	res = res.lstrip('\n')
	print (1111, res.split('\n')[:10])
	return res


def _restructure_tblasn_output(output_path):
	df = pd.read_csv(output_path, sep='\t|\s+', header=None)
	df.drop(df.iloc[:,10:13], inplace=True, axis=1)
	df.drop(df.iloc[:,11:], inplace=True, axis=1)
	df.iloc[:,-1] = df.iloc[:,-1].map(lambda x: x.rstrip(','))
	HeaderList = ['qseqid', 'sseqid', 'sstart', 'send', 'evalue', 'bitscore', 'length', 'pident', 'nident', 'stitle', 'chromosome']
	df.columns = HeaderList
	return df


def _location_list(dataframe, position1, position2):
        gene_result = []
        for i in range(len(dataframe)):
            a = dataframe.iloc[i][position1]
            b = dataframe.iloc[i][position2]
            list1 = [a,b]
            gene_result.append(list1)
        return gene_result


def _location_of_CDS(file_path, position1, position2):
	df = pd.read_csv(file_path, sep = '\t')
	chromosomes_index = set(df.loc[:,'chromosome'])
	gene_result_dict = {}
	for chromosome_index in chromosomes_index:
		name = 'chromosome_{}'.format(chromosome_index)
		gene_result = _location_list(
				df.loc[df['chromosome'] == chromosome_index],
				position1, position2
			)
		tmp = {name : gene_result}
		gene_result_dict.update(tmp)
	return gene_result_dict


def _check_outer_region(location_from_tblastn, locations_from_reference):
    a, b = location_from_tblastn
    for ref in locations_from_reference:
        c, d = ref
        if c < a < d or c < b < d or a < c < d < b or b < c < d < a:
            return False
    return True


def _filter_reference_proteins(fasta, amp_name):
	fasta_arr = fasta.split('\n')
	fasta_arr = [x for x in fasta_arr if x != '']
	length = len(fasta_arr)
	new_fasta = ''
	for i in range (0, length, 2):
		if fasta_arr[i].count(amp_name.lower()) > 0\
		or fasta_arr[i].count(amp_name.upper()) > 0\
		or fasta_arr[i].count(amp_name.title()) > 0:
			if not re.search(fasta_arr[i+1], new_fasta) \
			and not re.search ('dorsal', fasta_arr[i]):
				new_fasta += fasta_arr[i] + '\n'
				new_fasta += fasta_arr[i+1] + '\n'
	path = '{}/reference_proteins.txt'.format(
		OUTPUT_DIR
	)
	common.write_txt(path, new_fasta.rstrip('\n').lstrip('\n'))
	return path


def retrieve_known_gene_location(amp_name, species):
	Entrez.email = 'duyson1999@gmail.com'
	search_term = '{} AND {}'.format(amp_name, species)
	with Entrez.esearch(db='gene', term=search_term, retmax=1000) as tmp:
		res = Entrez.read(tmp)
	handle_first_record = Entrez.efetch(
		db="gene",id=res["IdList"],rettype="tabular",retmode = 'text'
	).read()
	result_path = '{}/retrieve_known_gene_location'.format(OUTPUT_DIR)
	common.write_txt(result_path, handle_first_record)
	return result_path


def collect_reference_proteins(
		amp_name,
		relative_species
	):
	ncbi_organisms = _format_species_query_ncbi(relative_species)
	uniprot_oragnisms  = _format_species_query_uniprot(relative_species)
	ncbi_data = _collect_reference_proteins_ncbi(
			amp_name,
			ncbi_organisms
		)
	uniprotdata = _collect_reference_proteins_uniprot(
			amp_name,
			uniprot_oragnisms
		)
	res = ncbi_data + '\n' + uniprotdata
	fasta = _restructure_fasta(res)
	path = _filter_reference_proteins(fasta, amp_name)
	return path


def index_genome_db(wgs_fa_path):
	saving_dir = '{}/wgs_db'.format(ROOT_DIR)
	common.mkdir(saving_dir)
	cmd = 'makeblastdb -in {} -dbtype nucl -out {}/DB'.format(wgs_fa_path, saving_dir)
	os.system(cmd)
	return '{}/DB'.format(saving_dir)


def query_tblastn(wgs_db_dir, reference_path):
	output_path = '{}/tblastn_result.txt'.format(OUTPUT_DIR)
	cmd = 'tblastn -db {} -query {} -out {} -word_size 3 \
		-xdrop_gap 32 -xdrop_ungap 16 -evalue 0.0001 \
		-outfmt "6 qseqid sseqid sstart send evalue bitscore length pident nident stitle"'.format(
			wgs_db_dir, reference_path, output_path
		)
	os.system(cmd)
	df = _restructure_tblasn_output(output_path)
	df.to_csv(output_path, sep='\t', index=False)
	return output_path


def separate_old_and_new_result(tblastn_output_path, known_genes_location_path):
	tBLASTn_dict = _location_of_CDS(
		tblastn_output_path, 'sstart', 'send'
	)
	known_location_dict = _location_of_CDS(
		known_genes_location_path,
		'start_position_on_the_genomic_accession',
		'end_position_on_the_genomic_accession'
	)
	new_CDS = {}
	old_CDS = {}
	for tblastn_chromosome in tBLASTn_dict:
		for reference_chromosome in known_location_dict:
			if list(known_location_dict.keys()).count(tblastn_chromosome) == 0:
				new_CDS[tblastn_chromosome] = tBLASTn_dict[tblastn_chromosome]
			else:
				if tblastn_chromosome == reference_chromosome:
					for CDS_location in tBLASTn_dict[tblastn_chromosome]:
						if _check_outer_region(CDS_location, known_location_dict[reference_chromosome]):
							list_CDS = new_CDS.get(tblastn_chromosome, [])
							list_CDS.append(CDS_location)
							new_CDS[tblastn_chromosome] = list_CDS
						else:
							list_CDS = old_CDS.get(tblastn_chromosome, [])
							list_CDS.append(CDS_location)
							old_CDS[tblastn_chromosome] = list_CDS
	return new_CDS, old_CDS


def retrieve_back_tBLASTn_result (CDS_dict, original_file, new_file_name):
	if len(CDS_dict) == 0:
		print ('no result in ' + new_file_name)
		open (new_file_name, 'a')
	else:
		original_df = pd.read_csv(original_file, sep='\t')
		print (original_df)
		final_df = pd.DataFrame({x: [] for x in original_df.columns})
		for chromosome in CDS_dict:
			for location in CDS_dict[chromosome]:
				temp_df = original_df.loc[
							(original_df['chromosome'] == int(chromosome.replace('chromosome_','')))
							& (original_df['sstart'] == location[0])
							& (original_df['send'] == location[1])
						]
				final_df = pd.concat([final_df, temp_df])
		final_df.to_csv(new_file_name, sep='\t', index=False)
	return


if __name__ == '__main__':
	wgs_db_dir = index_genome_db(WHOLE_GENOME_FASTA)
	known_genes_location_path = retrieve_known_gene_location(AMP, SPECIES_NAME)
	ref_path = collect_reference_proteins(AMP, RELATIVE_SPECIES)
	tblastn_path = query_tblastn(wgs_db_dir, '{}/reference_proteins.txt'.format(OUTPUT_DIR))
	new_CDS, old_CDS = separate_old_and_new_result(tblastn_path, known_genes_location_path)
	retrieve_back_tBLASTn_result(new_CDS, tblastn_path, '{}/{}'.format(OUTPUT_DIR, 'new_cds.txt'))
	retrieve_back_tBLASTn_result(old_CDS, tblastn_path, '{}/{}'.format(OUTPUT_DIR, 'old_cds.txt'))





