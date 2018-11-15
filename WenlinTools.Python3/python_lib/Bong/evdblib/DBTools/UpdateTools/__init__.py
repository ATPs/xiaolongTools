'''
This package contains modules related to the batch update of 
SCOP or PDB databases.

For update related to SCOP, SCOPUpdate 
and for update related to PDBUpdate provide convenient tools.

This package also contains tools independent from update types.
'''
import os, sys, shutil, glob, shlex

from evdblib.DBTools import Settings
from evdblib.DBTools import DomainInformation
from evdblib.DBTools.UpdateTools.QualityControl import check_profile_integrity
from evdblib.Utils.Parsers import PDB
from evdblib.Utils.Parsers.Range import PDBRange
from evdblib.Utils import build_profile_filename, parse_profile_filename, parse_sequence_filename
from evdblib.Utils.Parsers.PDBMappingTools import PDBMap
from evdblib.Runners.StructureTools import DaliLiteDATRunner

from evdblib.Utils.SimpleThreadingTools import ThreadPool
from evdblib.SunGridEngineTools import ComputeNodes
cm = None #ComputeNodes()
pool = None #ThreadPool( len(cm.get_hostname_list()) )

verbose = 0

class SearchDatabasePreparationError( Exception ) :
	pass

def prepare_input_sequence_for_profile_building( dominfo ) :
	'''
	build a sequence files for profile and save the information into 
	dominfo.
	
	New items "profile_sequence_file" and "profile_sequence_range" will
	be added into dominfo dictionary.

	Note that the profile sequence file and range will be same as input file if the
	db type is sequence. For structure, profile sequence file will be biologically
	relavent region defined by DBREF in PDB.
	'''

	data_type = Settings.get( "data_type" ) 
	reference_sequence_suffix = Settings.get( "reference_sequence_suffix" )
	profile_sequence_suffix = Settings.get( "profile_sequence_suffix" )
	
	profile_sequence_file = os.path.join( dominfo['domain_path'], dominfo['uniqueid'] ) + profile_sequence_suffix
	#reference_sequence_file = os.path.join( dominfo['domain_path'], dominfo['uniqueid'] ) + profile_sequence_suffix

	if data_type == 'sequence' :
		if dominfo['original_input_path']  :
			try :
				shutil.copyfile( dominfo['original_input_path'], profile_sequence_file )
				dominfo[ 'profile_sequence_file' ] = profile_sequence_file
				dominfo[ 'profile_sequence_range' ] = dominfo[ 'range' ]

			except IOError :
				dominfo[ 'profile_sequence_file' ] = ''
				dominfo[ 'profile_sequence_range' ] = ''
				print("WARNING: Profile sequence file cannot be written.", profile_sequence_file, file=sys.stderr)
		else :
			dominfo[ 'profile_sequence_file' ] = ''#dominfo[ 'original_input_path' ]
			dominfo[ 'profile_sequence_range' ] = ''#dominfo[ 'range' ]


	elif data_type == 'structure' and dominfo['original_input_path'] and dominfo['domain_path'] :
		pdb = PDB.parse( dominfo[ 'original_input_path' ] )
		pdbrange = PDBRange()
		pdbrange.parse( dominfo['range'] )

		chainrange = PDBRange()
		chainrange.parse( ','.join([ cid+':' for cid in pdbrange.get_unique_chain_ids() ]) )

		#Full SEQRES sequence 
		profile_sequence = pdb.extract_sequence( chainrange, biological=False, standard_residue_name=False, atomrecord=True, backbone=False )[0]
		#convert the pdbrange into the sequence range matching the given
		#set of the residue indications.
		profile_sequence_range = pdb.pdbrange2sequencerange( pdbrange, biological=False, standard_residue_name=True, atomrecord=True, backbone=False )
		###########################

		#print "*"*10, "dominfo"
		#print dominfo
		#print profile_sequence_range
		#for i,contig in enumerate(profile_sequence_range) :
			#print 'contig:%d'%i, contig.get_start(), contig.get_end()

		header = '>%s' % ( dominfo['uniqueid']  )

		try: 
			if os.path.exists( profile_sequence_file ) :
				raise IOError

			fp = open( profile_sequence_file, 'w' )
			print(header, file=fp)
			print(profile_sequence, file=fp)
			fp.close()

		except IOError :
			dominfo[ 'profile_sequence_file' ] = ''
			dominfo[ 'profile_sequence_range' ] = ''
			print("WARNING: Profile sequence file cannot be written.", profile_sequence_file, file=sys.stderr)

		else :
			dominfo[ 'profile_sequence_file' ] = profile_sequence_file
			dominfo[ 'profile_sequence_range' ] = str(profile_sequence_range)
	else :
		print("Error!", file=sys.stderr)
		print(dominfo, file=sys.stderr)
		raise Exception("Cannot build profile sequence file.")

from evdblib.Runners.StructureTools import DaliLiteDATRunnerError		
def normalize_and_map_structures_for_searches( dominfo ) :
	'''
	Build normalized input structure files for structure searches.

	1. Build PDB file containing only CA ATOM records 
	   (without Unknown residues, unknown or hetero residues or 
	    alternative positions)
	   All Non-standard amino acids (cannot be mapped to standard AA)
	   will be converted for Alanine.
	   #this file is mainly for FAST and TMalign.

	2. Build PDB file with full atomes of residues having all backbone 
	   atoms. This type of PDB files will be used for DaliLite.
	   DaliLite needs full backbone atoms to accept as a valid amino acid.

	3. build map between sequence and structures.

	Note:
	This function should be run only on protein structure.
	'''

	#setting pdb
	pdb_file = dominfo['original_input_path']
	pdb = PDB.parse( pdb_file )
	
	#preparing pdbrange
	pdbrange = PDBRange()
	pdbrange.parse( dominfo['range'] )

	#preparing save directory
	domain_path = dominfo['domain_path'] #we want to save to the domain_path

	#major mapping genera
	pdbmap = PDBMap( pdb, pdbrange )
	pdbmap.save_mapping_file( save_dir=domain_path, basename=dominfo['uniqueid'] )

	#DaliLite DAT
	backbonefn = os.path.join( domain_path, dominfo['uniqueid']+'.bb' )
	try :	
		datrunner = DaliLiteDATRunner( save_dir=domain_path, pdbfn=backbonefn )
		datrunner.run()
	except DaliLiteDATRunnerError :
		print("DaliLite DAT building failed!")
		


def prepare_profile_search_database( method, domain_informations=None, iteration=None, prefix='', dir=None, use_between=None, selected_iterations=None, strict=False ) :
	'''
	Prepares local condensed search database file 
	for convienent and fast search.

	Methods should be one of the search method defined
	in database configuration.

	Currently "COMPASS", "HHsearch"
	can be used as a method keyword.

	Optional domain_informations list can be given for
	building local db for subset of the current content in the database

	This function returns a filename of a composite local db filename
	and the number of records in the database.
	The number of records in the db is helpful especially for running 
	HHsearch to print out all the available alignments.

	Note that the use_between option True makes
	the local db building procedure use in between selected iterations
	if the specified "iteraion" does not exist.
	e.g. for iteration=3, selected_iterations = [1,3,5,8]
	the domain converged at iteration 2.
	Then the iteration 2 will be used if the "use_between" option is true. 
	
	By Default, the following options will use values specified in config file;
		use_between
		selected_iterations
		domain_informations
	.
	'''
	##############################
	#preparation for options
	##############################
	if domain_informations == None :
		domain_informations = DomainInformation.get_all_records()
		
	if strict and domain_informations == None :
		raise TypeError( "Domain information fetch failed." )

	if iteration == None :
		raise TypeError( "Integer value of iteration should be given." )

	try :
		iteration = int( iteration )
	except ValueError :
		raise ValueError( "Iteration should be integer or should be convetable to an integer value." )

	if dir == None :
		dir = Settings.get( "local_db_space" )
	
	local_db_root = os.path.join( dir, prefix ) 
		
	if not os.path.exists( local_db_root ) :
		os.makedirs( local_db_root )

	if use_between == None :
		use_between_string = Settings.get( "use_between_selected_iterations" )

		#print >>sys.stderr, 'use_between_string', use_between_string

		if use_between_string in [ 'true', 'True' ] :
			use_between = True
		elif use_between_string in ['false', 'False'] :
			use_between = False
		else :
			raise ValueError( "database configuration of use_between_selected_iterations is wrong: %s"%use_between_string )

	if selected_iterations == None :
		selected_iterations = [ int(i) for i in Settings.get( "selected_iterations" ).split() ] 


	db_filename = os.path.join( local_db_root, '.'.join( [prefix, method, str(iteration)] ) )
	##############################
	#actual db building code
	##############################
	if method == 'HHsearch' :
		number_of_records = _prepare_hhsearch_search_db( 
			db_filename, 
			domain_informations, iteration, 
			use_between, selected_iterations  )
	elif method == 'COMPASS' :
		number_of_records = _prepare_compass_search_db( 
			db_filename, 
			domain_informations, iteration, 
			use_between, selected_iterations  )
	else :
		raise TypeError( "Profile search database method should be HHsearch or COMPASS." )

	return db_filename, number_of_records

def integer2four_letter_id( i ) :
        '''convert unique number id into unique alphabet id'''
        n = int(i)

        if n/(26**4) or n <= 0:
		raise IndexError( i, "out of the range that can be converted" )

        n = n-1 #adjusting starting from 1
        # 1 -> aaaa, 456976 (or 26*26*26*26) -> zzzz

        a3 = n/(26**3)
        r = n%(26**3)
        a2 = r/(26**2)
        r = r%(26**2)
        a1 = r/26
        r = r%26
        a0 = r

        acount = 'abcdefghijklmnopqrstuvwxyz'
        return acount[a3]+acount[a2]+acount[a1]+acount[a0]


def generate_temporary_domain_ids( domain_informations ) :
	'''
	For Structural comparison programs it is quite necessary
	to have the normalized domain names.

	domain_informations is a list of dominfo, a dictionary of
	domain related infos.
	'''
	uniqueid2tempid = {}
	tempid2uniqueid = {}
	
	for i, dominfo in enumerate( domain_informations ) :
		uniqueid = dominfo['uniqueid']
		tempid = integer2four_letter_id( i )
		
		uniqueid2tempid[ uniqueid  ] = tempid
		tempid2uniqueid[ tempid ] = uniqueid

		#save it into dominfo directly!!!
		dominfo['tempid'] = tempid

	return uniqueid2tempid, tempid2uniqueid


def prepare_structure_search_database( method, domain_informations=None, prefix='', dir=None, strict=False, compute_node_dir=None ) :
	'''
	Prepares local condensed search database directory 
	for convienent and fast search.

	Methods should be one of the search method defined
	in database configuration.

	Currently "DaliLite", "FAST", "TMalign"
	can be used as a method keyword.

	Optional domain_informations list can be given for
	building local db for subset of the current content in the database
	'''

	##############################
	#preparation for options
	##############################
	if domain_informations == None :
		domain_informations = DomainInformation.get_all_records()
	else :
		pass
		
	if strict and domain_informations == None :
		raise TypeError( "Domain information fetch failed." )

	#specific directory to have all data files into one directory.
	if dir == None :
		dir = Settings.get( "local_db_space" )
	
	local_db_root = os.path.join( dir, prefix ) 
		
	if not os.path.exists( local_db_root ) :
		os.makedirs( local_db_root )

	#unlike to profile search case,
	#the main db path is directory name
	#not the filename.
	#since all of the methods are essentially pairwise
	db_dir = os.path.join( local_db_root, '.'.join([prefix, method]) )
	if os.path.exists( db_dir ) :
		shutil.rmtree( db_dir )
		
	if not os.path.exists( db_dir ) :
		os.makedirs( db_dir )

	##############################
	#actual db building code
	##############################
	if method == 'DaliLite' :
		number_of_records = _prepare_dalilite_search_db( 
			db_dir, domain_informations )
	elif method == 'FAST' :
		number_of_records = _prepare_fast_search_db( 
			db_dir, domain_informations )
	elif method == 'TMalign':
		number_of_records = _prepare_tmalign_search_db( 
			db_dir, domain_informations )
	else :
		raise TypeError( "Structure search database method should be DaliLite, FAST or TMalign.", method )
	
	################################
	#copying LOCAL_SCRATCH if it is set in the db.config.
	if compute_node_dir == None :
		compute_node_db_dir = Settings.get( "compute_node_db_space" )

	if compute_node_db_dir :
		compute_node_db_root = os.path.join( compute_node_db_dir, prefix )

	db_dir = _prepare_compute_node_db( db_dir, compute_node_db_root )


	return db_dir, number_of_records

def _prepare_computing_pools() :
	'''
	prepare cm and pool variable in the module
	'''
	global cm 
	cm = ComputeNodes()
	global pool 
	pool = ThreadPool( len(cm.get_hostname_list()) )


def _prepare_compute_node_db( local_db_dir, remote_db_dir ) :
	'''
	remotely prepare db files.
	to speed up access to files in compute nodes.
	'''
	if not remote_db_dir :
		return local_db_dir

	if not pool :
		_prepare_computing_pools()

	#remove previous junk
	cleanup = "rm -rf %s" % remote_db_dir
	cm.run_remote_command_thread_saving( shlex.split(cleanup), pool )

	#build up the directory tree
	parent_dir = ''
	for dir in remote_db_dir.split('/') :
		
		if dir == '' :
			continue
		else :
			parent_dir += '/' + dir
			cm.run_remote_command_thread_saving( ['mkdir', parent_dir], pool )
	else :
		#except the last dir.
		cm.run_remote_command_thread_saving( ['rmdir', remote_db_dir], pool )

	#copying the files
	cm.run_remote_command_thread_saving( ['cp', '-r', local_db_dir, remote_db_dir], pool )

	return remote_db_dir


#Internal functions for actual Local DB building procedures tailored to 
#each method.
def _prepare_dalilite_search_db( local_db_dir, domain_informations ) :
	'''
	copy DAT files into the local_db_dir directory.
	'''

	#dat = DaliLiteDAT() 
	for dominfo in domain_informations :
		uniqueid = dominfo[ 'uniqueid' ] 
		dompath = dominfo[ 'domain_path' ] 
		dat_path = os.path.join( dompath, uniqueid + '.dat' )

		if not os.path.exists( dat_path ) :
			raise IOError( dat_path, "Dali DAT file does not exist!" )

		new_dat_path = os.path.join( local_db_dir, uniqueid + '.dat' )
		if os.path.exists( new_dat_path ) :
			raise IOError( new_dat_path, "Dali DAT file already exists!" )

		shutil.copy( dat_path, new_dat_path )
	
		#blocked out due to the dalilite id mapping is done
		#in the DaliLite runner
		"""
		if not dominfo['tempid'] :
			#raise SearchDatabasePreparationError( "Each domain should be assigned to a temporary id!" )

		new_dat_path = os.path.join( local_db_dir, dominfo['tempid']+'A.dat' )

		#write the converted file
		dat.convert_identifier( input=dat_path, input_identifier=dominfo['tempid']+'A', 
					output=new_dat_path, output_identifier=dominfo['tempid']+'A' )
		"""
	return len(dominfo)


def _prepare_fast_search_db( local_db_dir, domain_informations ) :
	'''
	copy CA pdb files into the local_db_dir directory.
	'''

	for dominfo in domain_informations :
		uniqueid = dominfo[ 'uniqueid' ] 
		dompath = dominfo[ 'domain_path' ] 
		capdb_path = os.path.join( dompath, uniqueid + '.ca' )


		if not os.path.exists( capdb_path ) :
			raise IOError( capdb_path, "CA only PDB file does not exist!" )
	
		new_capdb_path = os.path.join( local_db_dir, parse_sequence_filename(capdb_path)[1]+".pdb" )

		if os.path.exists( new_capdb_path ) :
			raise IOError( new_capdb_path, "CA only PDB file already exists!" )

		#write the converted file
		shutil.copy( capdb_path, new_capdb_path )
	return len(dominfo)


def _prepare_tmalign_search_db( local_db_dir, domain_informations ) :
	'''
	copy CA pdb files into the local_db_dir directory.
	'''

	for dominfo in domain_informations :
		uniqueid = dominfo[ 'uniqueid' ] 
		dompath = dominfo[ 'domain_path' ] 
		capdb_path = os.path.join( dompath, uniqueid + '.ca' )


		if not os.path.exists( capdb_path ) :
			raise IOError( capdb_path, "CA only PDB file does not exist!" )
	
		new_capdb_path = os.path.join( local_db_dir, parse_sequence_filename(capdb_path)[1]+".pdb" )

		if os.path.exists( new_capdb_path ) :
			raise IOError( new_capdb_path, "CA only PDB file already exists!" )

		#write the converted file
		shutil.copy( capdb_path, new_capdb_path )
	return len(dominfo)


def _prepare_hhsearch_search_db( 
	db_filename,
	domain_informations, 
	iteration, 
	use_between, 
	selected_iterations ) :
	'''
	prepare HHsearch DB.
	and returns the number of records saved in the database file.
	'''
	db_fp = open( db_filename, 'w' )
	ext = Settings.get( 'hhm_suffix' )

	#getting previous iteration
	#for selecting iteration bigger than before.
	previous_iteration = selected_iterations[ max(selected_iterations.index( iteration ) - 1, 0) ] 

	number_of_records = 0
	for dominfo in domain_informations :

		#read domain path
		if 'domain_path' in dominfo :
			domain_path = dominfo[ 'domain_path' ]
		else :
			raise ValueError( 'domain_path does not exists', dominfo )

		if not domain_path :
			if verbose :
				print("WARNING: Dominfo does not have domain_path...")
				print(dominfo)
			continue

		domid = dominfo['uniqueid']
		hhsearch_file = build_profile_filename( domain_path, domid+'.prof', iteration, ext )

		if not os.path.exists( hhsearch_file ) and use_between :
		#in case the hhsearch file of the iteration
		#does not exists
		#and the use_between flag is On...
			#means the profile is generation is good!
			#and the value is max iteratoin!!
			last_available_iteration = check_profile_integrity( dominfo )
			if not last_available_iteration :
				print("WARNING: Profile is bad!", domain_path, domid)
				continue
			
			if last_available_iteration > previous_iteration :
				hhsearch_file = build_profile_filename( domain_path, domid+'.prof', last_available_iteration, ext )
			else :
				if verbose :
					print("WARNING: No between iteration available!", iteration, last_available_iteration)
				continue

		elif not os.path.exists( hhsearch_file ) and not use_between :
			#when hhsearch file is not availble and use  between flag is off.
			last_available_iteration = check_profile_integrity( dominfo )
			hhsearch_file = build_profile_filename( domain_path, domid+'.prof', last_available_iteration, ext )
			

		#final check!
		if not os.path.exists( hhsearch_file ) :
			#error!
			print("Error: HHsearch HMM file should be available but not found!", hhsearch_file, file=sys.stderr)
			raise SearchDatabasePreparationError( "HHsearch HHM file is not availble!", hhsearch_file )

		fp = open( hhsearch_file )
		content = fp.read()
		fp.close()

		db_fp.write( content )
		number_of_records += 1

	db_fp.close()
	return number_of_records

def _prepare_compass_search_db( 
	db_filename,
	domain_informations, 
	iteration, 
	use_between, 
	selected_iterations ) :

	'''
	prepare compass DB.
	'''

	db_fp = open( db_filename, 'w' )
	ext = Settings.get( 'compass_suffix' )

	db_size_fp = open( db_filename + ".len", 'w' ) #need to be built
	compass_db_size = 0

	previous_iteration = selected_iterations[ max( selected_iterations.index(iteration)-1, 0 ) ] 

	number_of_records = 0
	for dominfo in domain_informations :

		#read domain path
		domain_path = dominfo[ 'domain_path' ]
		if not domain_path :
			if verbose :
				print("WARNING: Dominfo does not have domain_path...")
				print(dominfo)
			continue

		domid = dominfo['uniqueid']
		compass_file = build_profile_filename( domain_path, domid+'.prof', iteration, ext )

		if not os.path.exists( compass_file ) and use_between :
		#in case the hhsearch file of the iteration
		#does not exists
		#and the use_between flag is On...
			#find the last iteration
			last_available_iteration = check_profile_integrity( dominfo )
			if not last_available_iteration :
				print("WARNING: Profile is bad!", domain_path, domid)
				continue
			
			if last_available_iteration > previous_iteration :
				compass_file = build_profile_filename( domain_path, domid+'.prof', last_available_iteration, ext )
			else :
				if verbose :
					print("No between iteration available!", iteration, last_available_iteration)
				continue

		elif not os.path.exists( compass_file ) and not use_between :
			#when hhsearch file is not availble and use  between flag is off.
			last_available_iteration = check_profile_integrity( dominfo )
			compass_file = build_profile_filename( domain_path, domid+'.prof', last_available_iteration, ext )

		#final check!
		if not os.path.exists( compass_file ) :
			#error!
			print("WARNING: COMPASS file should be available but not found!", compass_file, file=sys.stderr)
			raise SearchDatabasePreparationError( "COMPASS numerical profile file is not availble!", compass_file )
			continue

		fp = open( compass_file )
		content = fp.read()
		fp.close()

		db_fp.write( content )
		number_of_records += 1

		compass_size_file =  compass_file+".len"
		try :
			fp = open( compass_size_file )
			compass_db_size += int( fp.read().strip() )
			fp.close()
		except :
			print("WARNING: Cannot read compass profile size file.", compass_size_file)

	db_fp.close()
	print(compass_db_size, file=db_size_fp)
	db_size_fp.close()

	return number_of_records


