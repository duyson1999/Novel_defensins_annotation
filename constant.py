AMP = 'defensin'
SPECIES_NAME = 'Hermetia Illucens'

ORGANISM_NCBI = open('organisms_ncbi.txt', 'r').read()
ORGANISM_UNIPROT = open('organisms_uniprot.txt', 'r').read()
ROOT_DIR = '/Users/sonvo/Documents/1.All_for_Genetic_Project/test_newcode'
GENOME_INDEX_DIR = '/Users/sonvo/Documents/1.All_for_Genetic_Project/AMPs/Hermetia_illucens_WGS/genome_assemblies_genome_fasta/ncbi-genomes-2021-10-18/database/DB'
OUTPUT_DIR = '{}/{}'.format(ROOT_DIR, AMP)
UNIPROT_QUERY_LIMIT = 34