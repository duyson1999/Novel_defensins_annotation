AMP = 'defensin'
SPECIES_NAME = 'Hermetia Illucens'

ORGANISM_NCBI = open('organisms_ncbi.txt', 'r').read()
ORGANISM_UNIPROT = open('organisms_uniprot.txt', 'r').read()

ROOT_DIR = '/Users/sonvo/Documents/1.All_for_Genetic_Project/Hermetia_Illucens'
WHOLE_GENOME_FASTA = '{}/whole_genome.fna'.format(ROOT_DIR)
OUTPUT_DIR = '{}/{}'.format(ROOT_DIR, AMP)
UNIPROT_QUERY_LIMIT = 34