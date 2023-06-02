### Boas Pucker ###
### b.pucker@tu-bs.de ###

# some functions were copied from:
#KIPEs: https://doi.org/10.3390/plants9091103
#MYB_annotator: https://doi.org/10.1186/s12864-022-08452-5
#Pucker, B. & Iorizzo M. (2022): 10.1101/2022.02.16.480750

__cite__ = """Pucker, B., 2023. xxx"""

__version__ = "v0.1"

__usage__ = """
					Construction of species tree (""" + __version__ + """
					python3 phyl_constr.py
					--baits <MYB_SEQ_FILE>
					--in <MYB_CLASSIFICATION_FILE>
					--out <OUTPUT_DIR>
					
					optional:
					--mode <MODE_OF_ANALYSIS (cds|prot)>[cds]
					
					--blastp <PATH_TO_BLASTP>[blastp]
					--makeblastdb <PATH_TO_MAKEBLASTDB>[makeblastdb]
					--fasttree <PATH_TO_FastTree>[FastTree]
					--mafft <PATH_TO_MAFFT>[mafft]
					
					--minsim <MINIMAL_BLASTP_HIT_SIMILARITY>[80.0]
					--minlen <MINIMAL_BLASTP_HIT_LENGTH>[50]
					--repratio <MINIMAL_RATIO_OF_REPRESENTED_SPECIES>[1.0]
					--occupancy <MINIMAL_ALIGNMENT_OCCUPANCY>[0.1]
					
					bug reports and feature requests: b.pucker@tu-bs.de
					"""

import os, sys, subprocess
from operator import itemgetter
try:
	import hashlib
except ImportError:
	pass

# --- end of imports --- #

def load_sequences_clean( fasta_file, mode ):
	"""! @brief load candidate gene IDs from file """
	
	#illegal characters that will be removed or replaced: (, ),[,],|, /, @, space, tab, comma, simicolon
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		if "	" in header:
			header = header.split('	')[0]
		header = header.replace( "(", "-" ).replace( ")", "-" ).replace( "[", "-" ).replace( "]", "-" ).replace( "@", "-" ).replace( "/", "-" ).replace( "|", "-" ).replace( ",", "-" ).replace( ";", "-" )
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				if mode == "cds":
					if len( "".join( seq  ) ) % 3 == 0:	#load only proper CDS sequences
						sequences.update( { header: "".join( seq ) } )
				else:
					sequences.update( { header: "".join( seq ) } )
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
				if "	" in header:
					header = header.split('	')[0]
				header = header.replace( "(", "-" ).replace( ")", "-" ).replace( "[", "-" ).replace( "]", "-" ).replace( "@", "-" ).replace( "/", "-" ).replace( "|", "-" ).replace( ",", "-" ).replace( ";", "-" )
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		if mode == "cds":
			if len( "".join( seq  ) ) % 3 == 0:	#load only proper CDS sequences
				sequences.update( { header: "".join( seq ) } )
		else:
			sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_config_data( input_config_file, mode ):
	"""! @brief load data from config file """
	
	config_data = []
	with open( input_config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if len( parts ) > 2:
					if mode == "cds":
						config_data.append( { 'id': parts[0], 'spec': parts[1], 'cds_file': parts[2], 'cds': load_sequences_clean( parts[2], mode ) } )
					else:
						config_data.append( { 'id': parts[0], 'spec': parts[1], 'prot_file': parts[2], 'prot': load_sequences_clean( parts[2], mode ) } )
				else:
					sys.stdout.write( "ERROR: " + line + "\n" )
					sys.stdout.flush()
			line = f.readline()
	return config_data


def translate( seq, genetic_code ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
	seq = seq.upper()
	
	peptide = []
	for i in range( int( len( seq ) / 3.0 ) ):
		codon = seq[i*3:i*3+3]
		try:
			peptide.append( genetic_code[ codon ] )
		except:
			peptide.append( "*" )
	return "".join( peptide )


def transeq( cds ):
	"""! @brief translate given coding sequences """
	
	genetic_code = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'}
	peps = {}
	for key in list( cds.keys() ):
		if len( cds[ key ] ) % 3 != 0:
			sys.stdout.write( "ERROR: CDS length not multiple of 3 - " + key + "\t" + str( len( cds[ key ] ) ) + "\n" )
			sys.stdout.flush()
		peps.update( { key: translate( cds[ key ], genetic_code ) } )
	return peps


def load_best_BLAST_hits( blast_result_file, min_sim, min_len ):
	"""! @brief load best BLAST hit per query """
	
	blast_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > min_sim:
				if int( parts[3] ) > min_len:
					try:
						if float( parts[-1] ) > blast_hits[ parts[0] ]['score']:
							blast_hits[ parts[0] ] = { 'id': parts[1], 'score': float( parts[-1] ) }
					except KeyError:
						blast_hits.update( { parts[0]: { 'id': parts[1], 'score': float( parts[-1] ) } } )
			line = f.readline()
	return blast_hits


def construct_orthogroups( blast_results, ref_spec_ID, min_occupancy, number_of_species ):
	"""! @brief construct orthogroups for alignment and phylogenetic tree construction """
	
	specs = list( blast_results.keys() )	#list of species IDs (required to access BLAST results in dictionary)
	orthogroup_seed_seqs = list( blast_results[ ref_spec_ID ].keys() )	#list of sequences detected in reference species (seed for orthogroups)
	orthogroups = {}	#{ 'seqID1': { 'spec1': hit_in_spec1, 'spec2': hit_in_spec2, ... }, 'seqID2': { 'spec1': hit_in_spec1, 'spec2': hit_in_spec2, ... }, ... }
	for seq in orthogroup_seed_seqs:
		tmp_collection = {}
		for spec in specs:
			try:
				tmp_collection.update( { spec: blast_results[ spec ][ seq ] } )
			except KeyError:
				pass
		if len( tmp_collection.keys() ) / number_of_species >= min_occupancy:
			orthogroups.update( { seq: tmp_collection } )
	return orthogroups


def construct_corresponding_cds_aln_file( prot_aln, cds_collection, cds_aln_file ):
	"""! @brief construct corresponding CDS alignment file """
	
	with open( cds_aln_file, "w" ) as out:
		with open( prot_aln, "r" ) as f:
			header = f.readline().strip()[1:]
			line = f.readline()
			current_prot_seq = []
			while line:
				if line[0] == ">":
					spec, seqID = header.split('@')	#splitting protein alignment header to get species and sequence IDs
					cds = cds_collection[ spec ][ seqID ]	#get corresponding CDS
					cds_aln_seq = []	#construct a codon-aware aligned sequence
					i = 0
					current_prot_seq = "".join( current_prot_seq )
					for aa in current_prot_seq:
						if aa == "-":
							cds_aln_seq.append( "---" )	#add gaps corresponding to size of one codon
						else:
							cds_aln_seq.append( cds[ i*3:(i+1)*3 ] )	#get corresponding codon
							i += 1
					out.write( '>' + header + "\n" + "".join( cds_aln_seq ) + "\n" )
					header = line.strip()[1:]
					current_prot_seq = []
				else:
					current_prot_seq.append( line.strip() )
				line = f.readline()
	return cds_aln_file


def load_alignment( aln_file ):
	"""! @brief load alignment from input file """
	
	sequences, names = {}, []
	with open( aln_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update( { header: "".join( seq ) } )
				names.append( header )
				header = line.strip()[1:]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )
		names.append( header )
	return sequences, names


def aln_len_check( alignment ):
	"""! @brief check if all sequences in alignment have the same length """
	
	lengths = []
	for key in list( alignment.keys() ):
		lengths.append( len( alignment[ key ] ) )
	if len( list( set( lengths ) ) ) != 1:	#not all sequences have same length
		for key in list( alignment.keys() ):
			sys.stdout.write( key + " - len: " + str( len( alignment[ key ] ) ) + "\n" )
		sys.stdout.flush()
		sys.exit( "ERROR: sequences do not have the same length." )		


def alignment_trimming( aln_file, cln_aln_file, occupancy, sorting ):
	"""! @brief remove all alignment columns with insufficient occupancy """
	
	alignment, names = load_alignment( aln_file )
	
	if sorting:	#alphanumerical sorting of names if activated
		names = sorted( names )
	
	# --- if there is an alignment (expected case) --- #
	if len( names ) > 0:
		aln_len_check( alignment )	#perform alignment check (same sequence lengths?)
		
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( list(alignment.values())[0] ):
			counter = 0
			for key in names:
				if alignment[ key ][ idx ] != "-":
					counter += 1
			if counter / float( len( list(alignment.keys()) ) ) >= occupancy:
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in names:
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:
					new_seq.append( seq[ idx ] )
				new_seq  = "".join( new_seq )
				if new_seq.count('-') == len( new_seq ):	#exclude sequences that contain only gaps after trimming
					sys.stdout.write( "WARNING: only gaps remaining in sequence - " + key + " (sequence not included in output)\n" )
					sys.stdout.flush()
				else:
					out.write( ">" + key + '\n' + new_seq + '\n' )
	
	# --- just in case the alignment file is empty ---#
	else:
		sys.stdout.write( "WARNING: input file was empty (" + aln_file + ")\n" )
		sys.stdout.flush()
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def construct_concatenated_alignment_file( specs, clean_alignment_files, concatenated_alignment_file, mode ):
	"""! @brief construct a concatenated alignment file """
	
	concatenated_alignment = {}
	for spec in specs:
		concatenated_alignment.update( { spec: [] } )
	for aln_file in clean_alignment_files:
		alignment, names = load_alignment( aln_file )
		renamed_headers_alignment = {}
		for key in list( alignment.keys() ):
			renamed_headers_alignment.update( { key.split('@')[0]: alignment[ key ] } )
		for spec in specs:
			try:
				concatenated_alignment[ spec ].append( renamed_headers_alignment[ spec ] )
			except KeyError:
				if mode == "cds":
					concatenated_alignment[ spec ].append( "N"*len( list( renamed_headers_alignment.values() )[0] ) )	#add corresponding numbers of Ns if spec not represented for this gene
				else:
					concatenated_alignment[ spec ].append( "X"*len( list( renamed_headers_alignment.values() )[0] ) )	#add corresponding numbers of Xs if spec not represented for this gene
	with open( concatenated_alignment_file, "w" ) as out:
		for spec in specs:
			out.write( '>' + spec + "\n" + "".join( concatenated_alignment[ spec ] ) + "\n" )


def modify_names_in_tree( input_tree_file, output_tree_file, mapping_table ):
	"""! @brief replace short names in tree by original names """
		
	with open( input_tree_file, "r" ) as f:
		tree = f.read()
	
	for key in mapping_table.keys():
		tree = tree.replace( key, mapping_table[ key ] )
	
	with open( output_tree_file, "w" ) as out:
		out.write( tree )


#EXCLUDE sequences with multiple stop codons

def main( arguments ):
	"""! @brief run everything """
	
	input_config_file = arguments[ arguments.index('--in')+1 ]
	bait_sequence_file = arguments[ arguments.index('--baits')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--mode' in arguments:
		mode = arguments[ arguments.index('--mode')+1 ]
	else:
		mode = "cds"	#cds|pep
	
	if '--blastp' in arguments:
		blastp = arguments[ arguments.index('--blastp')+1 ]
	else:
		blastp = "blastp"	#path to BLASTp
	
	if '--makeblastdb' in arguments:
		makeblastdb = arguments[ arguments.index('--makeblastdb')+1 ]
	else:
		makeblastdb = "makeblastdb"	#path to makeblastdb
	
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
	
	if '--fasttree' in arguments:
		fasttree = arguments[ arguments.index('--fasttree')+1 ]
	else:
		fasttree = "FastTree"	#path to FastTree
	
	
	
	if '--minsim' in arguments:
		min_sim = float( arguments[ arguments.index('--minsim')+1 ] )
	else:
		min_sim = 80.0	#minimal similarity of a BLAST hit to be considered
	
	if '--minlen' in arguments:
		min_len = int( arguments[ arguments.index('--minlen')+1 ] )
	else:
		min_len = 50	#minimal length of a BLAST hit to be considered
	
	if '--repratio' in arguments:
		min_occupancy = float( arguments[ arguments.index('--repratio')+1 ] )
	else:
		min_occupancy = 1.0	#ratio of species that must be represented in an orthogroup
	
	if '--occupancy' in arguments:
		occupancy = float( arguments[ arguments.index('--occupancy')+1 ] )
	else:
		occupancy = 0.1	#minimal alignment occupancy
	
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	
	# --- load content of config file: ID, species name, CDS file (optional) --- #
	config_data = load_config_data( input_config_file, mode )	#[ { 'ID': "id1", 'spec': "species1", 'cds': "filex.fasta" } ]
	
	# --- construct corresponding PEP file --- #
	if mode == "cds":
		for idx, each in enumerate( config_data ):
			cds_file = output_folder + each['id'] + ".cds.fasta"	#clean CDS file
			if not os.path.isfile( cds_file ):
				with open( cds_file, "w" ) as out:
					for key in list( each["cds"].keys() ):
						out.write( '>' + key + '\n' + each["cds"][ key ] + "\n" )
			
			peps = transeq( each["cds"] )
			prot_file = output_folder + each['id'] + ".pep.fasta"	#clean PEP file
			if not os.path.isfile( prot_file ):
				with open( prot_file, "w" ) as out:
					for key in list( peps.keys() ):
						out.write( '>' + key + '\n' + peps[ key ] + "\n" )
			config_data[ idx ].update( { 'prot': peps, 'prot_file': prot_file } )	#adding the corresponding protein sequences
	
	# --- use BUSCO baits to identify single copy genes in all species (take best BLAST hit) --- #
	blast_results = {}
	for idx, each in enumerate( config_data ):
		blast_result_file = output_folder + each['id'] + ".bastp.txt"
		if not os.path.isfile( blast_result_file ):
			blastdb = output_folder + each['id'] + "_blastdb"
			p = subprocess.Popen( args= makeblastdb + " -in " + each['prot_file'] + " -out " + blastdb + " -dbtype prot", shell=True )
			p.communicate()
			
			p = subprocess.Popen( args= blastp + " -query " + bait_sequence_file + " -db " + blastdb + " -evalue 0.00001 -outfmt 6 -out " + blast_result_file, shell=True )
			p.communicate()
		best_hits = load_best_BLAST_hits( blast_result_file, min_sim, min_len )
		blast_results.update( { each['id']: best_hits } )
	
	
	# --- filter BLAST hits: present in all species, similarity in all species within certain range, length in all species within certain range --- #
	orthogroups = construct_orthogroups( blast_results, config_data[0]['id'], min_occupancy, float( len( config_data ) ) )
	
	# --- construct alignment for each gene (prot alignment) --- #
	alignment_folder = output_folder + "alignments/"
	if not os.path.exists( alignment_folder ):
		os.makedirs( alignment_folder )
	
	alignment_files = []
	for key in list( orthogroups.keys() ):
		seq_file = alignment_folder + key + ".fasta"
		with open( seq_file, "w" ) as out:
			for each in config_data:
				try:
					ID = orthogroups[ key ][ each['id'] ]['id']
					out.write( '>' + each['id'] + "@" + ID + "\n" + each['prot'][ ID ] + "\n" )	#write sequences of orthogroups into combined file
				except KeyError:	#species might not be represented in the alignment of this gene
					pass
		alignment_file = seq_file + ".prot.aln"
		tmp_file = seq_file + ".tmp"
		if not os.path.exists( alignment_file ):
			p = subprocess.Popen( args= mafft + " " + seq_file + " > " + alignment_file + " 2> " + tmp_file, shell=True )
			p.communicate()
		alignment_files.append( alignment_file )
	
	# --- optional: resubstitute prot alignment by CDS if available --- #
	if mode == "cds":
		cds_collection = {}
		for each in config_data:
			cds_collection.update( { each['id']: each['cds'] } )
		cds_alignment_files = []
		for prot_aln in alignment_files:
			cds_aln_file = prot_aln.replace( ".prot.aln", ".cds.aln" )
			cds_alignment_files.append( construct_corresponding_cds_aln_file( prot_aln, cds_collection, cds_aln_file ) )
	
	# --- perform alignment trimming --- #
	clean_alignment_files = []
	if mode == "cds":	#clean all CDS alignment files
		for aln_file in cds_alignment_files:
			cln_aln_file = aln_file + ".cln"
			alignment_trimming( aln_file, cln_aln_file, occupancy, sorting=False )
			clean_alignment_files.append( cln_aln_file )
	else:	#clean protein alignment files
		for aln_file in cds_alignment_files:
			cln_aln_file = aln_file + ".cln"
			alignment_trimming( aln_file, cln_aln_file, occupancy, sorting=False )
			clean_alignment_files.append( cln_aln_file )
	
	# --- concatenate alignments (fill in Ns for gaps) --- #
	concatenated_alignment_file = output_folder + "concatenated_alignment.fasta.aln.cln"
	if not os.path.isfile( concatenated_alignment_file ):
		specs = []
		for each in config_data:
			specs.append( each['id'] )
		construct_concatenated_alignment_file( specs, clean_alignment_files, concatenated_alignment_file, mode )
	
	
	# --- construct phylogenetic tree --- #
	tree_file = concatenated_alignment_file + ".tre"
	tmp_file = concatenated_alignment_file + ".tmp"
	if mode == "pep":
		cmd = " ".join( [ fasttree, "-wag -nosupport <", concatenated_alignment_file, ">", tree_file, "2>", tmp_file ] )	#protein tree
	else:
		cmd = " ".join( [ fasttree, "-gtr -nosupport <", concatenated_alignment_file, ">", tree_file, "2>", tmp_file ] )	#nucleotide tree
	p = subprocess.Popen( args=cmd, shell=True )
	p.communicate()
	
	
	# --- replace IDs with species names --- #
	output_tree_file = concatenated_alignment_file + ".FINAL_TREE.tre"
	mapping_table = {}
	for each in config_data:
		mapping_table.update( { each['id']: each['spec'] } )
	modify_names_in_tree( tree_file, output_tree_file, mapping_table )


if '--in' in sys.argv and '--baits' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
