########
# PURPOSE OF THIS CODE:
# 1. Retrieve protein queries from source (ncbi, uniprot) using for tBLASTn
# 2. Perform tBLASTn 
# 3. Seperate output file into newCDS and oldCDS - Based on 'gene' database in ncbi
#
#
# Note: the ncbi name of chromomome is different between species. 
# Edit the code to change chromosome name of result2.txt and reference_result.txt before running the code
 
#############################################################################################
AMP = str(input('AMP name: '))
Species = str(input('Organism name: '))

organism_ncbi = " “Ixodes“[Organism] OR  “Symphylella“[Organism] OR  “Glomeris“[Organism] OR  “Cypridininae“[Organism] OR  “Sarsinebalia“[Organism] OR  “Litopenaeus“[Organism] OR  “Celuca“[Organism] OR  “Lepeophtheirus“[Organism] OR  “Daphnia“[Organism] OR  “Speleonectes“[Organism] OR  “Acerentomon“[Organism] OR  “Sminthurus“[Organism] OR  “Tetrodontophora“[Organism] OR  “Anurida“[Organism] OR  “Polqonognathellus“[Organism] OR  “Folsomia“[Organism] OR  “Campodea“[Organism] OR  “Occasjapyx“[Organism] OR  “Meinertellus“[Organism] OR  “Machilis“[Organism] OR  “Tricholepidion“[Organism] OR  “Thermobia“[Organism] OR  “Atelura“[Organism] OR  “Calopteryx“[Organism] OR  “Cordulegaster“[Organism] OR  “Epiophlebia“[Organism] OR  “Baetis“[Organism] OR  “Isonychia“[Organism] OR  “Eurylophella“[Organism] OR  “Ephemera“[Organism] OR  “Zorotypus“[Organism] OR  “Forﬁcula“[Organism] OR  “Apachyus“[Organism] OR  “Leuctra“[Organism] OR  “Perla“[Organism] OR  “Cosmioperla“[Organism] OR  “Gryllotalpa“[Organism] OR  “Ceuthophilus“[Organism] OR  “Tetrix“[Organism] OR  “Prosarthria“[Organism] OR  “Stenobothrus“[Organism] OR  “Tanzaniophasma“[Organism] OR  “Galloisiana“[Organism] OR  “Grylloblatta“[Organism] OR  “Happoembia“[Organism] OR  “Aposthonia“[Organism] OR  “Timema“[Organism] OR  “Peruphasma“[Organism] OR  “Aretaon“[Organism] OR  “Metallyticus“[Organism] OR  “Empusa“[Organism] OR  “Mantis“[Organism] OR  “Blaberus“[Organism] OR  “Periplaneta“[Organism] OR  “Cryptocercus“[Organism] OR  “Mastotermes“[Organism] OR  “Prorhinotermes“[Organism] OR  “Zootermopsis“[Organism] OR  “Gynaikothrips“[Organism] OR  “Frankliniella“[Organism] OR  “Thrips“[Organism] OR  “Trialeurodes“[Organism] OR  “Bemisia“[Organism] OR  “Acanthocasuarina“[Organism] OR  “Planococcus“[Organism] OR  “Essigella“[Organism] OR  “Acyrthosiphon“[Organism] OR  “Aphis“[Organism] OR  “Acanthosoma“[Organism] OR  “Notostira“[Organism] OR  “Ranatra“[Organism] OR  “Velia“[Organism] OR  “Xenophysella“[Organism] OR  “Nilaparvata“[Organism] OR  “Cercopis“[Organism] OR  “Okanagana“[Organism] OR  “Ectopsocus“[Organism] OR  “Liposcelis“[Organism] OR  “Menopon“[Organism] OR  “Pediculus“[Organism] OR  “Tenthredo“[Organism] OR  “Orussus“[Organism] OR  “Cotesia“[Organism] OR  “Leptopilina“[Organism] OR  “Nasonia“[Organism] OR  “Chrysis“[Organism] OR  “Acromyrmex“[Organism] OR  “Harpegnathos“[Organism] OR  “Exoneura“[Organism] OR  “Apis“[Organism] OR  “Bombus“[Organism] OR  “Inocellia“[Organism] OR  “Xanthostigma“[Organism] OR  “Corydalus“[Organism] OR  “Sialis“[Organism] OR  “Conwentzia“[Organism] OR  “Osmylus“[Organism] OR  “Pseudomallada“[Organism] OR  “Euroleon“[Organism] OR  “Mengenilla“[Organism] OR  “Stylops“[Organism] OR  “Aleochara“[Organism] OR  “Dendroctonus“[Organism] OR  “Meloe“[Organism] OR  “Tribolium“[Organism] OR  “Lepicerus“[Organism] OR  “Priacma“[Organism] OR  “Gyrinus“[Organism] OR  “Carabus“[Organism] OR  “Rhyacophila“[Organism] OR  “Platycentropus“[Organism] OR  “Hydroptila“[Organism] OR  “Philopotamus“[Organism] OR  “Annulipalpia chim“[Organism] OR  “Micropterix“[Organism] OR  “Dyseriocrania“[Organism] OR  “Triodia“[Organism] OR  “Nemophora“[Organism] OR  “Yponomeuta“[Organism] OR  “zyyaena“[Organism] OR  “Polyommatus“[Organism] OR  “Parides“[Organism] OR  “Bombyx“[Organism] OR  “Manduca“[Organism] OR  “Ceratophyllus“[Organism] OR  “Archaeopsylla“[Organism] OR  “Ctenocephalides“[Organism] OR  “Boreus“[Organism] OR  “Nannochorista“[Organism] OR  “Bittacus“[Organism] OR  “Panorpa“[Organism] OR  “Anoppheles“[Organism] OR  “Aedes“[Organism] OR  “Phlebotomus“[Organism] OR  “Trichocera“[Organism] OR  “Tipula“[Organism] OR  “Bibio“[Organism] OR  “Bombylius“[Organism] OR  “Drosophila“[Organism] OR  “Lipara“[Organism] OR  “Rhagoletis“[Organism] OR  “Glossina“[Organism] OR  “Sarcophaga“[Organism] OR  “Triarthria“[Organism] OR  “Myrmica“[Organism] "
organism_uniprot = " Ixodes OR  Symphylella OR  Glomeris OR  Cypridininae OR  Sarsinebalia OR  Litopenaeus OR  Celuca OR  Lepeophtheirus OR  Daphnia OR  Speleonectes OR  Acerentomon OR  Sminthurus OR  Tetrodontophora OR  Anurida OR  Polqonognathellus OR  Folsomia OR  Campodea OR  Occasjapyx OR  Meinertellus OR  Machilis OR  Tricholepidion OR  Thermobia OR  Atelura OR  Calopteryx OR  Cordulegaster OR  Epiophlebia OR  Baetis OR  Isonychia OR  Eurylophella OR  Ephemera OR  Zorotypus OR  Forﬁcula OR  Apachyus OR  Leuctra OR  Perla OR  Cosmioperla OR  Gryllotalpa OR  Ceuthophilus OR  Tetrix OR  Prosarthria OR  Stenobothrus OR  Tanzaniophasma OR  Galloisiana OR  Grylloblatta OR  Happoembia OR  Aposthonia OR  Timema OR  Peruphasma OR  Aretaon OR  Metallyticus OR  Empusa OR  Mantis OR  Blaberus OR  Periplaneta OR  Cryptocercus OR  Mastotermes OR  Prorhinotermes OR  Zootermopsis OR  Gynaikothrips OR  Frankliniella OR  Thrips OR  Trialeurodes OR  Bemisia OR  Acanthocasuarina OR  Planococcus OR  Essigella OR  Acyrthosiphon OR  Aphis OR  Acanthosoma OR  Notostira OR  Ranatra OR  Velia OR  Xenophysella OR  Nilaparvata OR  Cercopis OR  Okanagana OR  Ectopsocus OR  Liposcelis OR  Menopon OR  Pediculus OR  Tenthredo OR  Orussus OR  Cotesia OR  Leptopilina OR  Nasonia OR  Chrysis OR  Acromyrmex OR  Harpegnathos OR  Exoneura OR  Apis OR  Bombus OR  Inocellia OR  Xanthostigma OR  Corydalus OR  Sialis OR  Conwentzia OR  Osmylus OR  Pseudomallada OR  Euroleon OR  Mengenilla OR  Stylops OR  Aleochara OR  Dendroctonus OR  Meloe OR  Tribolium OR  Lepicerus OR  Priacma OR  Gyrinus OR  Carabus OR  Rhyacophila OR  Platycentropus OR  Hydroptila OR  Philopotamus OR  Annulipalpia chim OR  Micropterix OR  Dyseriocrania OR  Triodia OR  Nemophora OR  Yponomeuta OR  zyyaena OR  Polyommatus OR  Parides OR  Bombyx OR  Manduca OR  Ceratophyllus OR  Archaeopsylla OR  Ctenocephalides OR  Boreus OR  Nannochorista OR  Bittacus OR  Panorpa OR  Anoppheles OR  Aedes OR  Phlebotomus OR  Trichocera OR  Tipula OR  Bibio OR  Bombylius OR  Drosophila OR  Lipara OR  Rhagoletis OR  Glossina OR  Sarcophaga OR  Triarthria OR  Myrmica "

# Class = str(input('Species family: '))
# create folder
import os
path = '/mnt/2a34c0e2-53cf-4ed8-ad5f-f8e31e63377f/1.All_for_Genetic_Project/Code_for_testing/' + AMP
os.mkdir(path)      #if the file is already existed, code will be error

#############################################################################################
# Search and Download sequences
    ## From ncbi
print ('##################################################')
print ('Getting data from ncbi...')
import Bio
from Bio import Entrez
Entrez.email = 'duyson1999@gmail.com'
search_term = AMP + '[All Fields]' + ' AND ' '(' + organism_ncbi + ')'
search = Entrez.esearch(db= 'protein', term = search_term, retmax=200, usehistory = 'y', idtype = 'acc')
search_result = Entrez.read(search)
search.close()
idlist = ','.join(search_result['IdList'])
record = Entrez.efetch(db = 'protein', id = idlist, rettype = 'fasta', retmode = 'text')
    ## From Uniprot
print('Getting data from Uniprot...')
from bioservices import UniProt
u = UniProt()
res = u.search(AMP + ' (' + organism_uniprot + ')', limit =200,  frmt = 'list')
uniprotdata = ''
for id in res.rstrip('\n').split('\n'):
    fastafile = u.retrieve(id, frmt='fasta')
    uniprotdata += fastafile
path2 = path + '/' + AMP + '.txt'
with open(path2, 'w') as f:
	f.write(record.read() + '\n' + uniprotdata) 
print ('--> Output file: ' + AMP + '.txt')
print ('##################################################')
############################################################################################
#Edit fasta file
print('Pre-processing dataset...')
a = open(path2).readlines()
b = ''
for line in a:
    if line.startswith('>'):
        line = line.replace('>', '\n>')
        b += line
    else: 
        line = line.strip()
        b += line 
b = b.lstrip('\n')
path3 = path + '/' + 'temperary_file.txt'
open(path3,'w').write(b)
x = open(path3).readlines()
    ## Filter querry sequences file
import re
length = len(x)
new_fasta = ''
for i in range (1, length, 2):
	if x[i-1].count(AMP) >0 or x[i-1].count(AMP.upper()) >0 or x[i-1].count(AMP.title()) > 0:
		if not re.search(x[i], new_fasta) and not re.search ('dorsal', x[i-1]):
			new_fasta += x[i-1]
			new_fasta += x[i]
path4 = path + '/' + 'newfile.txt'
open(path4, 'w').write(new_fasta)
print('---> Output file: newfile.txt')
print ('##################################################')

##############################################################################################
#Perform Blast: 
print('Running BLAST...')
path5 = path + '/result'
cmd = 'tblastn -db /mnt/2a34c0e2-53cf-4ed8-ad5f-f8e31e63377f/1.All_for_Genetic_Project/AMPs/Hermetia_illucens_WGS/genome_assemblies_genome_fasta/ncbi-genomes-2021-10-18/database/DB -query ' + path4 + ' -out ' + path5 + ' -word_size 3 -xdrop_gap 32 -xdrop_ungap 16 -evalue 0.0001 -outfmt "6 qseqid sseqid sstart send evalue bitscore length pident nident stitle"'
os.system(cmd)
print ('---> Output file: result')import
print ('##################################################')

##############################################################################################
### Retrieve reference gene location from ncbi
print ('Retrieve reference CDS location of AMP...')
Entrez.email = 'duyson1999@gmail.com'
search_term = AMP + ' AND ' + Species 
handle = Entrez.esearch(db="gene",term = search_term, retmax = 1000)
record = Entrez.read(handle)
handle.close()
#print (record)
handle_first_record = Entrez.efetch(db="gene",id=record["IdList"],rettype="tabular",retmode = 'text')
path6 = path + '/reference_result.txt'
open(path6, 'w').write(handle_first_record.read())
print ('---> Output file: reference_result.txt')
print ('##################################################')

################################################################################################
# Retrieve location from reference tabular file
def location_of_CDS(filenamee, column1, column2):
    df = pd.read_csv (filenamee, sep = '\t')
    aa = set(df.loc[:,'chromosome'])
    def location_list (dataframe):
        gene_result = []
        for i in range(len(dataframe)):
            a = dataframe.iloc[i][column1]
            b = dataframe.iloc[i][column2]
            list1 = [a,b]
            gene_result.append(list1)
        return gene_result
    gene_result_dict = {}
    for element in aa:
        name = 'chromosome_' + str(element)
        d = {name : location_list(df.loc[df['chromosome'] == element])}
        gene_result_dict.update(d)
    return gene_result_dict
################################################################################################
os.chdir(path)
if os.stat('result').st_size == 0:  # If there no hits, no need to do more
    print ('no hits from tBLASTn')
else:
    ### Edit output file (result from tblastn)
    ### add header, modify 'result' file, create 'result2.txt'
    print ("Edit 'result' file.....")
    import pandas as pd
    df2 = pd.read_csv('result', sep= '\t|\s+', header=None)
    df2.drop(df2.iloc[:,10:13], inplace=True, axis=1)
    df2.drop(df2.iloc[:,11:], inplace=True, axis=1)
    df2.iloc[:,-1] = df2.iloc[:,-1].map(lambda x: x.rstrip(','))
    HeaderList = ['qseqid', 'sseqid', 'sstart', 'send', 'evalue', 'bitscore', 'length', 'pident', 'nident', 'stitle', 'chromosome']
    df2.to_csv('result2.txt', sep='\t', header=HeaderList, index = False)
    print ('---> Output file: result2.txt')
    print ('##################################################')

    ################################################################################################
    # Create dict to contain location of CDS from result2.txt and reference_result,txt
    print ("Get location of CDS from 'result2.txt' and 'reference_result.txt'...")
    tBLASTn_dict = location_of_CDS('result2.txt', 'sstart', 'send')
    Reference_dict = location_of_CDS('reference_result.txt', 'start_position_on_the_genomic_accession', 'end_position_on_the_genomic_accession')
    reference_name_list = []
    for name in Reference_dict: 
        reference_name_list.append(name)

    ################################################################################################
    # Compare, list out new putative hits, old putative hits
    ## Compare 2 ranges
    def comparing (location_from_tblastn, location_from_reference):
        A = True
        a,b = [element1 for element1 in location_from_tblastn]              # location from tblastn result
        for list in location_from_reference:
            c,d = [y for y in list]                                         # Location from ncbi reference
            if c < a < d or c < b < d or a < c < d < b or b < c < d < a: 
                A = False
        return A
    ## Group hits into new_putative and old_putative
    new_CDS = {}
    old_CDS = {}
    for tblastn_chromosome in tBLASTn_dict:
        for reference_chromosome in Reference_dict:
            if reference_name_list.count(tblastn_chromosome) == 0:
                new_CDS[tblastn_chromosome] = tBLASTn_dict[tblastn_chromosome]
            else:
                if tblastn_chromosome == reference_chromosome:
                    for CDS_location in tBLASTn_dict[tblastn_chromosome]:
                        if comparing(CDS_location, Reference_dict[reference_chromosome]):
                            list_CDS = new_CDS.get(tblastn_chromosome, [])
                            list_CDS.append(CDS_location)
                            new_CDS[tblastn_chromosome] = list_CDS
                        if not comparing(CDS_location, Reference_dict[reference_chromosome]):
                            list_CDS = old_CDS.get(tblastn_chromosome, [])
                            list_CDS.append(CDS_location)
                            old_CDS[tblastn_chromosome] = list_CDS
    ## Re-search back 'result2.txt' file
    print ('Grouping new_putative_CDS and old_putative_CDS')
    df3 = pd.read_csv('result2.txt', sep = '\t')
    def Retrieve_Back_tBLASTn_result (CDS_dict, original_file, new_file_name):
        if len(CDS_dict) == 0:
            print ('no result in ' + new_file_name)
            open (new_file_name, 'a')
        else:
            for element1 in CDS_dict:
                for element2 in CDS_dict[element1]:
                    originaldf = pd.read_csv(original_file, sep = '\t')
                    newdf = originaldf.loc[originaldf['chromosome'] == element1.replace('chromosome_','')].loc[originaldf['sstart']== element2[0]].loc[originaldf['send']== element2[1]]
                    f = open(new_file_name, 'a')
                    newdf.to_csv(f, header = None, mode ='a', sep='\t')
            tempora_df = pd.read_csv(new_file_name, sep='\t', header = None).iloc[:,1:]
            tempora_df.to_csv(new_file_name, header=False , sep='\t')
            tempora_df.to_csv(new_file_name, header= list (df3.columns.values), sep='\t')
    # if repeat the code, there will be more headers !!!!!!        
    Retrieve_Back_tBLASTn_result(new_CDS, 'result2.txt','New_CDS.txt')
    Retrieve_Back_tBLASTn_result(old_CDS, 'result2.txt','Old_CDS.txt')
    print ("---> Output: 'New_CDS.txt' and 'Old_CDS.txt'")
