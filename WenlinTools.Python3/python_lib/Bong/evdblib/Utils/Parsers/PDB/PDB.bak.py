'''
PDB package main module.

Code example:
pdb = PDB.parse( fn='12as.pdb' ) #generating PDB instance
pdb_null = PDB.parse() #will generate a PDB class instance with No data in it.

'''

from evdblib.Utils.Parsers.PDB import *
from .PDBInfo import PDBInfo
from .Models import Models

class PDB :
	'''
	PDB main class.
	
	This class contains instances of the Header and Model classes.
	The two classes are the main workhorse for understanding and
	parsing the PDB file.

	Note that PDB module and PDB class is different.
	If you want to directly use without using PDB module's parse command
	you should refer PDB by PDB.PDB.
	'''
	def __init__( self, fn='', fp=None ) :
		#initializing variables
		self.fn = fn
		self.pdbinfo = None
		self.models = None

		#variable for saving selection 
		#to extract sequences/coordinates
		self.selections = []

		##########
		#use fn when the filename is given.
		if fn :
			self.fp = open( fn )
		else :
			self.fp = fp

		self.input_stream_start_point = self.fp.tell()

		########################
		#start parsing if 
		#an input stream is defined!
		if self.fp :
			self.parse()
		########################
	
	def __getitem__( self, chainid ) :
		'''
		returns chain from the default model.
		'''
		return self.get_model()[chainid ]
		

	def get_fp( self ) :
		'''
		returns file pointer of the current input stream.

		This function returns None 
		when there is no input stream currently set.

		Designed for internal usage.
		'''
		return self.fp


	def rewind_fp( self ) :
		'''
		rewinds the file pointer of PDB input.
		and returns the pointer

		Designed for internal usage.
		'''
		self.fp.seek( self.input_stream_start_point )
		return self.fp


	def parse( self ) :
		'''
		parses the PDB file or input stream defined 
		at the class initialization time.

		Note that this function will raise NullError when
		no file or no file pointer is given at the initialize time.
		'''
		if not self.fp :
			raise NullError( "No input stream is given." )

		#self or PDB class instance pointer should be given
		#to the header class so that the header class can
		#add information to PDB class.
		#header class will add various information to the PDB class
		try :
			self.pdbinfo = PDBInfo( pdb=self )
		except IndexError :
			self.pdbinfo = PDBInfo()
			
		self.rewind_fp()

		#model class will add various information to the PDB class to.
		self.models = Models( pdb=self ) #PDB class should be given.


	def get_exp_type( self ) :
		'''
		returns experimental type of the structure.
		'''
		if self.pdbinfo :
			return self.pdbinfo.get_experiment()
		else :
			return ''


	def get_resolution(self ) :
		'''
		returns resolution of the structure for experiments
		that have resolution defined.

		If no resolution is defined or found by the title section,
		it returns -1.
		'''
		if self.pdbinfo :
			return self.pdbinfo.get_resolution()
		else :
			return -1.0


	def get_pdb_id( self ):
		'''
		returns PDB ID parsed from the header.
		'''
		if self.pdbinfo :
			return self.pdbinfo.pdb_id
		else :
			return ''

	def get_model( self, model_id='default' ) :
		'''
		returns the model.

		By default, the "default" MODEL will be returned. 
		The defualt model is the representative model 
		if specified in REMARK 210 section,
		otherwise the model with lowest model id.
		'''
		return self.models.get( model_id )

	def get_protein_chain_ids( self ) :
		''' 
		returns Chain IDs for protein, or polypeptide chains only.
		This function is convenient for choosing only protein chains.
		
		Note:
		Current code determines the protein chains by checking SEQRES records.
		It might be rewritten to check based actual ATOM records to be more
		robust.
		'''
		chain_ids = []
		for chain_id, chaininfo in self.pdbinfo.sequenceinfos.items() :
			if chaininfo.is_protein() :
				chain_ids.append( chain_id )

				if verbose : print(chain_id, "is a protein.")
			else :
				if verbose : print(chain_id, "is not a protein.")

		return chain_ids

	def get_residues( self, chain_id='', model_id=-1, pdbrange=None, residue_type='biological' ) :
		'''
		returns selected residues by the chain 
		(with model if model_id is specified)
		or defined by the pdbrange object.

		This function is the main sequence extractor from the chain.
		'''

		#not implemented yet!
		if pdbrange :
			for irange in pdbrange :
				irange.start() #this should return 
		

	def get_atom_sequence_modres_corrected( self, model_id=-1, chain_id='A' ) :
		'''
		generate sequence for a given chain_id.
		'''
		chaininfo = self.remarks.chaininfos[ chain_id ]
		chain = self.get_model(model_id).get_chain(chain_id)
		residue_ids = [residue.id for residue in self.get_residue_selection( model_id=model_id, chain_id=chain_id )]
		seq_fp = cStringIO.StringIO()
		for residue_id in residue_ids :
			residue = chain.residues[residue_id]
			if not residue.is_belong_to_chain() :
				continue
				
			resname = residue.get_name()
			if resname in amino_acids :
				seq_fp.write( amino_acids[resname] )
			else :
				std_resname = chaininfo.get_standard_resname( residue_id )
				if std_resname :
					if std_resname in amino_acids :
						seq_fp.write( amino_acids[std_resname]  )
					else :
						print("WARNING: MODRES record have non-standard residue name!", resname, residue_id, self.get_pdb_id(), self.pdb_fn, file=sys.stderr)
						seq_fp.write( 'X'  )
				else :
					print("WARNING: MODRES record does not have this modification", self.pdb_fn, residue_id, resname, file=sys.stderr)
					seq_fp.write( 'X' )
		return seq_fp.getvalue()
		
	def print_chain_summary( self, model_id=-1 ) :
		'''
		prints chain_id, # of missing residues, # of engineered residues, experimental method, resolution, average B-factor
		
		mainly developed for debugging.
		'''
		model = self.get_model( model_id )
		print("ChainID	length	# of missing residues Name")
		for chain_id in sorted(self.pdbinfo.get_chain_ids()) :
			chain = model.get_chain( chain_id )
			chaininfo = chain.seqinfo
			print(chain_id,"\t", len(chain.extract_sequence(atomrecord=True)),"/",len(chain.get_seqres_sequence()), "\t", len(chaininfo.get_missing_residue_ids( model_id )), "\t",chaininfo.get_name())
			print(chain.get_bfactor(), self.get_resolution(), self.get_exp_type())

	def get_uniq_longest_protein_chain_ids( self ) :
		'''
		returns list of chain ids from unique molecules.
		The selected chains have longest chains or coverages of
		the molecules. And if the molecules have same number of atoms,
		then the chains having lower B-factor 

		This function is useful for getting chains for template library.
		'''
		#protein chain ids to filter protein ids.
		protein_chain_ids = set(self.get_protein_chain_ids())

		longests = []
		mol_ids, id_lists = self.pdbinfo.moleculeinfo.get_molecule_id_and_chain_groups()
		for id_list in id_lists :
			id_list2 = protein_chain_ids.intersection( id_list )
			if not id_list2 :
				continue

			chain_id = self.get_longest_protein_chain_id( id_list2 ) 
			longests.append( chain_id )
	
		return longests

	def get_longest_protein_chain_id( self, chain_ids ) :
		'''
		returns longest chain id.
		If tie occurs, the tie will be broken by
		selecting lowest B-factor.
		'''
		id_list = []
		for chain_id in chain_ids :
			chain = self.get_model().get_chain( chain_id )
			id_list.append( (-len(chain.extract_sequence(atomrecord=True)), chain.get_bfactor(), chain_id) )
		if id_list :
			return min(id_list)[-1]
		else :
			return None

	def get_lowest_bfactor_chain_id( self, chain_ids ) :
		'''
		returns lowest bfactor chain_id.
		If error occurs returns None.
		'''
		bfactor_list = []
		for chain_id in chain_ids :
			chain=self.get_model().get_chain( chain_id )
			bfactor_list.append( (chain.get_bfactor(), chain_id) )
			
		if bfactor_list :
			return min(bfactor_list)[-1]
		else :
			return None
			

	def get_uniq_lowest_bfactor_protein_chain_ids( self ) :
		'''
		returns list of chain ids from unique molecules.
		The selected chains have lowest b-factors.

		This function is useful for getting chains for template library.
		'''
		#protein chain ids to filter protein ids.
		protein_chain_ids = set(self.get_protein_chain_ids())

		lowests = []
		mol_ids, id_lists = self.pdbinfo.moleculeinfo.get_molecule_id_and_chain_groups()
		for id_list in id_lists :
			id_list2 = protein_chain_ids.intersection( id_list ) 
			if not id_list2 :
				continue

			chain_id = self.get_longest_protein_chain_id( id_list2 ) 
			lowests.append( chain_id )
	
		return lowests



	def select_residues_by_scop_range( self, scop_range, model_id=0 ) :
		'''
		This function select range of scop
		'''
		ranges = scop_range.split(',')
		for r in ranges :
			chain_id = None
			start_resnum = None
			end_resnum = None
			start_insertion_code = ' '
			end_insertion_code = ' '
			chain_and_numbers = r.split(':')
			chain_id = chain_and_numbers[0]
			if chain_and_numbers[1] :
				number_array = chain_and_numbers[1].split('-')
				numbers = []
				icodes = [] #need to incorporate insertion code part. :(
				while number_array :
					i = number_array.pop(0)
					if i :
						pass
					else : #this part is written for differentiating negative sign (-) 
					       #from the range -. Not that readable. :(
						i = number_array.pop(0)
						i = '-'+i
					#dealing with insertion code
					if i[-1].isdigit() :  #no insertion code
						numbers.append( int(i) )
						icodes.append( ' ' )
					else : 	#with insertion code following digits
						numbers.append( int(i[:-1]) )
						icodes.append( i[-1] )
				start_resnum, end_resnum = numbers
				start_insertion_code, end_insertion_code = icodes
					
			if not self.select_residues( model_id=model_id, chain_id=chain_id, start_resnum=start_resnum, end_resnum=end_resnum, start_insertion_code=start_insertion_code, end_insertion_code=end_insertion_code) :
				return False 
		return True
	
			
	def check_selected_residue( self, residue ) :
		for selection in self.residue_selections :
			if residue in selection :
				return True
		return False

	def select_residues( self, model_id=0, chain_id='A', start_resnum=None, end_resnum=None, start_insertion_code = ' ', end_insertion_code = ' ' ) :
		'''
		saves selection given by the 
		'''
		selection = self.get_residue_selection( model_id, chain_id, start_resnum, end_resnum, start_insertion_code, end_insertion_code )
		if selection :
			self.residue_selections.append( selection )
		return True

	def get_residue_selection_by_residue_id( self, start_id, end_id, model_id=0 ) :
		selection = []
		for chain_id in self.get_model( model_id ).chains.keys() :
			chain = self.get_model( model_id ).get_chain(chain_id)
			for residue_ids in chain.get_within_chain_residue_ids() :
				if start_id <= residue_id  <= end_id :
					selection.append( residue_id )
		return selection
			
			
	def get_residue_selection( self, model_id, chain_id, start_resnum=None, end_resnum=None, start_insertion_code= ' ', end_insertion_code=' ' ) :
		selection = []
		chain = self.get_model( model_id ).get_chain( chain_id )
		if not chain :
			print("Warning: No Chain", chain_id, "exists.", file=sys.stderr)
			return False
		for residue_id in chain.get_within_chain_residue_ids() :
			residue = chain.get_residue( residue_id )
			#resnum = residue.get_number()
			
			start_residue_id = None
			end_residue_id = None
			if start_resnum != None :
				start_residue_id = build_residue_id( chainid, start_resnum, start_insertion_code )
			if end_resnum != None :
				end_residue_id = build_residue_id( chainid, end_resnum, end_insertion_code )
				
			if start_residue_id == None or residue.id >= start_resnum :
				if end_residue_id == None or residue.id <= end_resnum :
					if not self.check_selected_residue( residue ) :
						selection.append( residue )
					else :
						print('WARNING: Duplicated Selection', residue_id, file=sys.stderr)
		return selection
	
	def reset_selection( self ) :
		self.residue_selections = []

	def selection2sequence( self ) :
		str_fp = cStringIO.StringIO()
		if verbose :
			print(self.residue_selections)
		for selection in self.residue_selections :
			for residue in selection :
				#print residue.id, residue.get_name(), residue.get_name_letter()
				str_fp.write( residue.get_name_letter() )
		return str_fp.getvalue()
	
	def selection2fragment_range_string( self ) :
		fragment_lengths = []
		if verbose :
			print("@ selection2fragment_range_string")
			
		for selection in self.residue_selections :
			str_fp = cStringIO.StringIO()
			for residue in selection :
				str_fp.write( residue.get_name_letter() )
			fragment_lengths.append( len(str_fp.getvalue()) )
		if not self.residue_selections :
			print("Error! No residues selected!", file=sys.stderr)
			return ""
		if verbose :
			print('fragment_lengths:',fragment_lengths)
		line = '0-' + str(fragment_lengths[0])
		current_offset = fragment_lengths[0]
		for length in fragment_lengths[1:] :
			line += ' ' + str(current_offset) + '-'+str(length+current_offset)
			current_offset += length
		return line

