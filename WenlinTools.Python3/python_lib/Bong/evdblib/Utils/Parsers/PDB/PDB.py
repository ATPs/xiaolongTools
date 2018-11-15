'''
PDB package main module.

Code example:
pdb = PDB.parse( fn='12as.pdb' ) #generating PDB instance
pdb_null = PDB.parse() #will generate a PDB class instance with No data in it.

'''
import sys

from evdblib.Utils.Parsers.PDB import *
from .PDBInfo import PDBInfo
from .Models import Models
from .Atom import Atom
from ..Range import SequenceRange, SequenceContiguousRange, PDBRange, PDBContiguousRange

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
			self.fp.close()
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

	

	
	#Deprecated
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

		Note that this function relies on the PDB Header.
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
			if chain :
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
			if chain :
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

	def pdbrange2sequencerange( self, pdbrange, *args, **kwargs ) :
		'''
		Convert PDBRange like SCOP domain definition into
		SequenceRange defining blocks of contiguous regions
		with sequential index values (1 based).

		The args should be None, it is just a place holder required 
		by the Python grammer.

		The kwargs are same as the extract_sequence function.

		Note that the PDBRange can shrink depending on the kwargs.
		For example, the PDBRange including His TAG at the N-terminus,
		and the kwarg "biological=True" is given, which indicates the 
		N-terminal HIS TAG removal, the sequence range will not include
		the TAG region.

		More generally, the kwargs define the sequence that should be conceptually 
		considered as a original protein sequence, 
		resulting sequence range is a domain defintion based on 
		the original protein sequence.
		'''
		residues = []
		uniq_chain_ids = []
		#saves start and end of chain residue indices
		#I.e. for X:,Y: case
		#{'X':(0,100), 'Y':(100,200)}
		#assuming that 100 residues in both chains X and Y.
		start_and_end_of_chains = {}

		#prepare non-redundent list of residues comprising 
		#multiple chains
		for contig in pdbrange :
			chain = self[contig.chain_id]
			if contig.chain_id in uniq_chain_ids :
				pass
			else :
				uniq_chain_ids.append( contig.chain_id )
				start_index = len(residues)
				#extract full seqres residues
				extracted_residues = chain.extract_residues(atomrecord=False, 
							biological=False, backbone=False)

				residues.extend( extracted_residues )
				end_index = len(residues)-1
				if end_index <= start_index :
					raise ValueError( "Extracting chain residues has a problem!" )
				start_and_end_of_chains[contig.chain_id] = (start_index, end_index)

		#convert list of residues into list of trues and falses
		#to build blocks of residues according to the given 
		#parameters in **kwargs
		binary_decision_values = self.filter_residues( residues, **kwargs )
		selected_residue_blocks = self.blockify( binary_decision_values )

		#converting indices for each blcok
		pdb_residue_blocks = []
		for contig in pdbrange :
			start_rid = contig.get_start()
			end_rid = contig.get_end()

			if start_rid == None :
				start_rid = residues[start_and_end_of_chains[contig.chain_id][0]].id
			if end_rid == None :
				end_rid = residues[start_and_end_of_chains[contig.chain_id][1]].id

			#converting indices from the PDBRange
			start_index = self.residue_id2index( start_rid, residues )
			if start_index == -1 :
				raise RangeConversionError( "Start Residue ID not found! %s" % str(start_rid) )
			end_index = self.residue_id2index( end_rid, residues )
			if end_index == -1 :
				raise RangeConversionError( "End Residue ID not found! %s" %str(end_rid) )
			if end_index < start_index :
				raise RangeConversionError( "End Residue Index is bigger than Start Residue!" )

			#blocks are assumed to be 0 based slices as in python
			pdb_residue_blocks.append( (start_index, end_index+1) ) 

		#overlap between PDBrange blocks and residue filtering blocks
		overlapping_blocks = self.get_overlap_between_blocks( selected_residue_blocks, pdb_residue_blocks )

		#debug
		#print selected_residue_blocks
		#print pdb_residue_blocks
		#print overlapping_blocks

		filtered_residues = []
		for residue, bin in zip( residues, binary_decision_values ) :
			if bin :
				filtered_residues.append( residue )

		sequencerange = SequenceRange()
		for block in overlapping_blocks :
			#getting Residue IDs from original full length seqres based list.
			start_rid = residues[ block[0] ].id
			end_rid = residues[ block[1]-1 ].id #adjust slice index

			start_index = self.residue_id2index( start_rid, filtered_residues )
			if start_index == -1 :
				raise RangeConversionError( "Start Residue ID not found! %s" % str(start_rid) )
			end_index = self.residue_id2index( end_rid, filtered_residues )
			if end_index == -1 :
				raise RangeConversionError( "End Residue ID not found! %s" %str(end_rid) )
			if end_index < start_index :
				raise RangeConversionError( "End Residue Index is bigger than Start Residue!" )

			#1 based blocks like in SCOPrange
			sequencerange.add_contiguous_range( SequenceContiguousRange( start_index+1, end_index+1) )

		return sequencerange


	def filter_residues( self, residues, *args, **kwargs ) :
		'''
		Convert Residue list into a list of Trues and Falses
		according to the give kwargs (key word arguments).

		If kwargs does not contain certain keywords, the following
		default values will be used.

		Default values
			'atomrecord': True
			'biological': True
			'backbone': False

		Note that this function has same filtering as extract_residues.
		'''
		bin_list = []

		#setting up the default values
		if 'atomrecord' in kwargs :
			atomrecord = kwargs['atomrecord']
		else :
			atomrecord = True

		if 'biological' in kwargs :
			biological = kwargs['biological']
		else :
			biological = True

		if 'backbone' in kwargs :
			backbone = kwargs[ 'backbone' ]
		else :
			backbone = False

		for residue in residues :
			#atomrecord
			if atomrecord :
				if residue.is_missing_residue(): 
					bin_list.append( False )
					continue

			if biological :
				if not residue.is_biological_residue() :
					bin_list.append( False )
					continue

			if backbone :
				if not residue.has_full_backbone_atoms() :
					bin_list.append( False )
					continue
		
			bin_list.append( True )

		return bin_list

	def blockify( self, binary_decision_list ) :
		'''
		Convert binary decision_list (Trues and Falses) into contiguous block index arrays

		For example, 

		01234567890
		FTTTTFFTTFF
		
		will be converted as the following
		[(1,5), (7,9)]
		.
		'''

		blocks = []
		in_block = False
		for i, bin in enumerate(binary_decision_list) :
			if bin and not in_block :
				block_start = i
				in_block = True
			elif bin and in_block :
				continue
			elif not bin and in_block :
				blocks.append( (block_start, i) )
				in_block = False
			elif not bin and not in_block :
				continue
			else :
				print("This message should not be seen!", file=sys.stderr)
				sys.exit(-1)
		else :
			if in_block :
				blocks.append( (block_start, i+1) )
				
		return blocks
		
			
	def get_overlap( self, start1, end1, start2, end2 ) :
		'''
		Finds overlapping indices between indice pairs of block.

		Written for the internal use to facilitate the pdbrange2sequencerange.
		'''

		start = max(start1, start2)
		end = min( end1, end2 )

		if end <= start :
			return None
		else :
			return (start, end)


	def get_overlap_between_blocks( self, blocks1, blocks2 ) :
		'''
		Finds overlapping indices between two sets of blocks.
	
		Written for the internal usage.
		'''
		blocks = []
		for block1 in blocks1 :
			for block2 in blocks2 :
				overlap = self.get_overlap( block1[0], block1[1], block2[0], block2[1] )
				if overlap :
					blocks.append( overlap )

		return blocks
		
			
	def residue_id2index( self, residue_id, residue_list ) :
		'''
		converts the residue id into the index value in the resdiue list.
		If no matching residue_id is found, -1 will be returned.
		
		This is helper function for pdbrange2sequencerange.
		'''


		for i, residue in enumerate(residue_list) :
			if residue.id == residue_id :
				return i
		else :
			return -1


	def extract_sequence( self, pdbrange, *args, **kwargs ) :
		'''
		returns a sequence falls within the range boundary.
		and an alignment between the reference SEQRES sequence and
		the extracted sequence.

		This function works as shallow wrapper for Chain.extract_sequence function.
		PDBRange class will be passed to it
		and Keywords are same as Chain.extract_sequence.
		'''

		fp = cStringIO.StringIO()
		rfp = cStringIO.StringIO()
		sfp = cStringIO.StringIO()

		for contig in pdbrange :
			chain = self[contig.chain_id]
			#print contig
			#print contig.get_start(), chain.residue_id2residue_index( contig.get_start() )
			#print contig.get_end(), chain.residue_id2residue_index( contig.get_end() )
			fp.write( chain.extract_sequence( start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), **kwargs ) )
		
			ref, seq = chain.generate_sequence_mapping( start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), **kwargs )
			
			rfp.write( ref )
			sfp.write( seq )
			

		return fp.getvalue(), rfp.getvalue(), sfp.getvalue()


	def extract_seqres_residues( self, pdbrange, model_id='default' ) :
		'''
		returns a list of residues falls within the range boundary.
		Note that the list supposed to contain all residues in
		SEQRES even if they are missing residues IF the PDB header
		information is correct.
		'''
		residues = []
		model = self.get_model( model_id )

		for contig in pdbrange :
			chain = model[contig.chain_id]
			residues.extend( 
				chain.extract_residues( 
					atomrecord = False,
					biological = False,
					start_residue_id=contig.get_start(), 
					end_residue_id=contig.get_end() 
				) 
			)

		return residues

	def extract_residues( self, pdbrange, model_id='default', **kwargs ) :
		'''
		returns a list of residues falls within the range boundary.

		This function is supposed to return list of residues
		that is subset of the SEQRES list.
		
		Optional arguments are same as Chain.extract_residues.
		briefly speaking this function can select biological, or 
		residues with actual coordinates, and etc.
		'''
		residues = []
		model = self.get_model( model_id )

		for contig in pdbrange :
			chain = model[contig.chain_id]
			residues.extend( 
				chain.extract_residues( 
					start_residue_id=contig.get_start(), 
					end_residue_id=contig.get_end() ,
					**kwargs
				) 
			)

		return residues




	def extract_ca_atom_records( self, pdbrange, use_one_chainid=True, *args, **kwargs ) :
		'''
		returns CA ATOM records falls within the the range boundary.
		Also reference seqres sequence and ca atom based sequence are returned.

		The sequences are returned for the mapping purpose.

		Note that this function is a shallow wrapper of 
		Chain.extract_ca_atom_records function in Chain module.
		keyword arguments are same as the Chain.extract_ca_atom_records function.
		'''

		afp = cStringIO.StringIO()
		rfp = cStringIO.StringIO()
		mfp = cStringIO.StringIO()
		start_atom_number=1
	
		if use_one_chainid == True :
			one_chainid = pdbrange[0].chain_id
		elif use_one_chainid :
			one_chainid = use_one_chainid[0]
	
		for contig in pdbrange :
			chain = self[contig.chain_id]

			if use_one_chainid :
				ca_records, ref, mseq = chain.extract_ca_atom_records(start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), start_number=start_atom_number, chainid=one_chainid, **kwargs )
			else :
				ca_records, ref, mseq = chain.extract_ca_atom_records(start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), start_number=start_atom_number, **kwargs )

			start_atom_number += ca_records.count('\n')

			afp.write( ca_records )
			rfp.write( ref )
			mfp.write( mseq )

		return afp.getvalue(), rfp.getvalue(), mfp.getvalue()


	def extract_atom_records( self, pdbrange, use_one_chainid=True, *args, **kwargs ) :
		'''
		returns ATOM records falls within the the range boundary.
		Also reference seqres sequence and atom based sequence are returned.
		The sequences are returned for the mapping purpose.

		Note taht this function is a shallow wrapper of 
		Chain.extract_atom_records function in Chain module.
		Keyword arguments are same as the Chain.extract_ca_atom_records function.

		if use_one_chain option is true,
		also use_one_chain option can be a ChainID.
		then chain ID will be the modified to the value.
		'''

		afp = cStringIO.StringIO()
		rfp = cStringIO.StringIO()
		mfp = cStringIO.StringIO()
		start_atom_number=1
		start_residue_number=1

		if use_one_chainid == True :
			one_chainid = pdbrange[0].chain_id
		elif use_one_chainid :
			one_chainid = use_one_chainid[0]
		
		for contig in pdbrange :
			chain = self[contig.chain_id]
			if use_one_chainid :
				records, ref, mseq = chain.extract_atom_records(start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), start_number=start_atom_number, start_residue_number=start_residue_number, chainid=one_chainid, **kwargs )
			else :
				records, ref, mseq = chain.extract_atom_records(start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), start_number=start_atom_number, start_residue_number=start_residue_number,**kwargs )
				

			start_atom_number += records.count('\n')
			if records :
				try :
					start_residue_number += int( records.rsplit('\n',2)[-2][22:26] )
				except IndexError :
					pass
				except TypeError :
					print("WARNING: problem during extract atom records.")
					start_residue_number += records.count('\n')

			afp.write( records )
			rfp.write( ref )
			mfp.write( mseq )

		return afp.getvalue(), rfp.getvalue(), mfp.getvalue()

	def extract_atom_records_for_dalidat( self, pdbrange, use_one_chainid='A', *args, **kwargs ) :
		'''
		returns ATOM records falls within the the range boundary.
		Also reference seqres sequence and atom based sequence are returned.
		The sequences are returned for the mapping purpose.

		Note taht this function is a shallow wrapper of 
		Chain.extract_atom_records function in Chain module.
		Keyword arguments are same as the Chain.extract_ca_atom_records function.

		if use_one_chain option is true,
		also use_one_chain option can be a ChainID.
		then chain ID will be the modified to the value.
		'''

		afp = cStringIO.StringIO()
		rfp = cStringIO.StringIO()
		mfp = cStringIO.StringIO()
		start_atom_number=1
		start_residue_number=1

		if use_one_chainid == True :
			one_chainid = pdbrange[0].chain_id
		elif use_one_chainid :
			one_chainid = use_one_chainid[0]
		
		for contig in pdbrange :
			chain = self[contig.chain_id]
			if use_one_chainid :
				records, ref, mseq = chain.extract_atom_records(start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), start_number=start_atom_number, start_residue_number=start_residue_number, chainid=one_chainid, icode=' ', altloc=' ', **kwargs )
			else :
				records, ref, mseq = chain.extract_atom_records(start_residue_id=contig.get_start(), end_residue_id=contig.get_end(), start_number=start_atom_number, start_residue_number=start_residue_number, icode=' ', altloc=' ', **kwargs )
				

			start_atom_number += records.count('\n')
			if records :
				try :
					start_residue_number += int( records.rsplit('\n',2)[-2][22:26] )
				except IndexError :
					pass
				except TypeError :
					print("WARNING: problem during extract atom records.")
					start_residue_number += records.count('\n')

			afp.write( records )
			rfp.write( ref )
			mfp.write( mseq )

		else :
			atom = Atom()
			atom.parse( atom_line=records.rsplit('\n',2 )[-2] )
			afp.write( atom.get_ter_record(number=atom.get_number()+1) )
			print("END", file=afp)

		return afp.getvalue(), rfp.getvalue(), mfp.getvalue()




class RangeConversionError( Exception ) :
	pass
