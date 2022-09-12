AMP = 'defensin'
SPECIES_NAME = 'Hermetia Illucens'

RELATIVE_SPECIES = open('relative_species.txt').read().splitlines()

ROOT_DIR = '/Users/sonvo/Documents/1.All_for_Genetic_Project/Hermetia_Illucens'
WHOLE_GENOME_FASTA = '{}/whole_genome.fna'.format(ROOT_DIR)
OUTPUT_DIR = '{}/{}'.format(ROOT_DIR, AMP)
UNIPROT_QUERY_LIMIT = 34