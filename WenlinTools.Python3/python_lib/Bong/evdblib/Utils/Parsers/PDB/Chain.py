import sys, io

from evdblib.Utils import is_sorted
from evdblib.Utils.Parsers.PDB import resnames2sequence, verbose, SEQRESMappingError
from evdblib.Utils.Parsers.PDB import build_residue_id, build_one_bigger_residue_id, amino_acids, amino_acids_rev

from .Residue import Residue

#verbose=1


class Chain :
	'''
	Chain class contains list of residues that are in the same chain.
	'''
	def __init__( self, chainid, model=None )  :
		self.model = model #pointer to the model that contains this chain
		self.chainid = chainid #chain id
		self.molid = None   #molecule ID in Header
		self.molname = None #molecule name in Header
		self.residue_ids = []
		self.residues = {}
		self.seqinfo = None #pointer to SequenceInfo class instance
		
		self.terminated = False #variable finding the info for the chain

		#flag for correct mapping
		self.no_sequence_mapping = False
		if model :
			self.no_sequence_mapping = model.no_sequence_mapping

		##################################
		#SEQRES mapping will be done later.
		#Just keeping the record of mapping.
		##################################
		self.seqres_atom_mapping_flag = 0 
		# 2: mapping done sucessfully covered from alignment procedure
		# 1: mapping done successfully!!, 
		# 0: mapping not done, 
		#-1: mapping is vague, due to missing MISSING remarks!!
		#-2: mapping has error even after the alignment procedure

		###########################
		#MODRES mapping will be done very late.
		#When the residue.get_standard_amino_acid_name() is called,
		#MODRES mapping for the residue will be done implicitly.
		###########################


		#setting information about the chain
		if model and model.pdb and model.pdb.pdbinfo :
			pdbinfo = model.pdb.pdbinfo
			self.molid = pdbinfo.get_molecule_id( self.chainid )
			self.molname = pdbinfo.get_molecule_name( self.chainid )
			self.seqinfo = pdbinfo.get_sequence_info( self.chainid )

			#need to be called after the building is completed!!
			#self.map_seqres_record()
			#self.seqinfo.set_biological_sequence( self )
			#currently those thow important functions will be called after the building
			#of model is done in Model class!


		#to reduce loading time
		#the mapping of missing sequence
		#to other sequences will be done
		#when the mapping is necessary.

		#######################
		#deprecated variables
		#######################
		#self.residues_within_chain = set()

		#all Residue IDs (including missing residue ids) 
		#are supposed to be here mapped
		#self.all_residue_ids = [] 

		#SEQRES index is mapping to ATOM index, 
		#Note that missing are mapped to the last ordered ATOM index
		#self.seqres_index2atom_index = []  #note that this is 1 based array!!!!

	def set_biological_sequence(self) :
		'''
		proxy function for the set_biological_sequence in seqinfo object.
		Due to error handling this is added here.
		'''
		if self.seqinfo :
			self.seqinfo.set_biological_sequence( self )


	def is_seqres_atom_mapping_done_correctly( self ) :
		''' returns one of the following values
		2: mapping done sucessfully covered from alignment procedure
		1: mapping done successfully!!
		0: mapping has not been done yet.
		-1: mapping is done but vaguely, due to missing MISSING remarks!!
		-2: mapping has error even after the alignment procedure
		'''
		return self.seqres_atom_mapping_flag
		

	#It might need to think about hetero elements
	def write_record( self, fp ) :
		for rid in self.residue_ids :
			residue = self.get_residue( rid )
			residue.write_record( fp )

	
	def get_seqres_resnames( self, start=None, end=None ) :
		'''
		returns list of 3 letter residue names in SEQRES record from start to end.
		
		Note that the start and end indices should be given python style.
		'''
		if self.seqinfo :
			return self.seqinfo.seqres_sequence[start:end]
		else :
			return []

	def get_seqres_sequence( self, start=None, end=None ) :
		'''
		returns one letter amino acid sequence derived from seqres record.
		'''
		seqres_resnames = self.get_seqres_resnames( start=start, end=end )
		#print seqres_resnames
		if seqres_resnames :
			return resnames2sequence( seqres_resnames )
		else :
			return ''

	def extract_residues( self, atomrecord=True, biological=True, standard_residue_name=True, backbone=False,
				start_residue_id=None, end_residue_id=None ) :
		'''
		returns list of Residues.
		
		If atomrecord is true, the list will be based on ATOM records, or else
		the list will be basing on SEQRES record.
		If biological is true, the list will be biologically relevant,
		most likely to be the seqeunce without vector contaminants.
		If standard residue name is true, the list content does not change since
		this option is developed for the extract_sequence function.
		If backbone is true, the residues do not have all the backbone residues will be
		removed from the list.
		
		Also the boundaries of the list can be determined by the start residue id
		and the end residue id.
		
		This is the best list generation tool for general usages.
		'''
		start, end = self.residue_id2residue_index(start_residue_id), self.residue_id2residue_index(end_residue_id)
		#adjust end index for the slicing.
		if end != None :
			end = end + 1

		residues = []
		for residue_id in self.residue_ids[start:end] :
			residue = self.get_residue(residue_id)
			if not residue.is_belong_to_chain() :
				continue

			if atomrecord and residue.is_missing_residue():
				continue
			
			if biological and not residue.is_biological_residue() :
				continue
		
			if backbone and not residue.has_full_backbone_atoms() :
				continue

			residues.append( residue )

		return residues
		

	def extract_sequence( self, atomrecord=True, biological=True, standard_residue_name=True, backbone=False,
				start_residue_id=None, end_residue_id=None ) :
		'''
		returns sequence of ATOM records or SEQRES record.
		
		If atomrecord is true, the sequence will be based on ATOM records, or else
		the sequence will be basing on SEQRES record.
		If biological is true, the sequence will be biologically relevant,
		most likely to be the seqeunce without vector contaminants.
		If standard residue name is true, the sequence will be unmodified standard residues
		found in the sequence database.
		If backbone is true, the residues do not have all the backbone residues will be
		removed from the sequence.
		
		Also the boundaries of the sequence can be determined by the start residue id
		and the end residue id.
		
		This is the best sequence generation tool for general usages.
		'''
		start, end = self.residue_id2residue_index(start_residue_id), self.residue_id2residue_index(end_residue_id)
		#adjust end index for the slicing.
		if end != None :
			end = end + 1

		fp = io.StringIO()
		for residue_id in self.residue_ids[start:end] :
			residue = self.get_residue(residue_id)
			if not residue.is_belong_to_chain() :
				continue

			if atomrecord and residue.is_missing_residue():
				continue
			
			if biological and not residue.is_biological_residue() :
				continue
		
			if backbone and not residue.has_full_backbone_atoms() :
				continue

			if standard_residue_name :
				fp.write( residue.get_standard_amino_acid_name_letter() )
			else :
				fp.write( residue.get_name_letter() )

		return fp.getvalue()
		

	def extract_ca_atom_records( self, biological=True, standard_residue_name=True, backbone=False, start_residue_id=None, end_residue_id=None, start_number=None, **kwargs ) :
		'''
		returns ATOM records, reference seqres based sequence, and aligned sequence of CA ATOM records to the SEQRES records.
		
		If biological is true, the records will be biologically relevant,
		most likely to be the ATOM records without vector contaminants.
		If standard residue name is true, the record will be unmodified standard residues
		found in the sequence database instead of the original line.
		If backbone is true, the residues do not have all the backbone residues will be
		removed from the records.
		
		Also the boundaries of the records can be determined by the start residue id
		and the end residue id.
		'''

		atomrecord = True
		start, end = self.residue_id2residue_index(start_residue_id), self.residue_id2residue_index(end_residue_id)
		if end != None : #adjust end index for slicing.
			end = end + 1

		atom_serial_number=1
		residue_serial_number=1

		if start_number != None :
			atom_serial_number = start_number
			residue_serial_number = start_number

		rfp = io.StringIO() #reference sequence fp
		fp = io.StringIO() #sequence fp
		afp = io.StringIO() #atom record fp

		for residue_id in self.residue_ids[start:end] :
			residue = self.get_residue(residue_id)
			if not residue.is_belong_to_chain() :
				continue

			if atomrecord and residue.is_missing_residue():
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue
			
			if biological and not residue.is_biological_residue() :
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue
		
			if backbone and not residue.has_full_backbone_atoms() :
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue

			if standard_residue_name :
				rfp.write( residue.get_standard_amino_acid_name_letter() )
				if residue.has_ca():
					fp.write( residue.get_standard_amino_acid_name_letter() )
					afp.write( residue.get_ca_record( use_standard_aa=True, number=atom_serial_number, resnum=residue_serial_number, **kwargs) )
					atom_serial_number += 1
					residue_serial_number += 1
				else :
					fp.write( '=' )
			else :	
				rfp.write( residue.get_name_letter() )
				if residue.has_ca() :
					fp.write( residue.get_standard_amino_acid_name_letter() )
					afp.write( residue.get_ca_record( use_standard_aa=False, number=atom_serial_number, resnum=residue_serial_number, **kwargs) )
					atom_serial_number += 1
					residue_serial_number += 1
				else :
					fp.write( '=' )

		return afp.getvalue(), rfp.getvalue(), fp.getvalue()
		

	def extract_atom_records( self, biological=True, standard_residue_name=True, backbone=False, start_residue_id=None, end_residue_id=None, start_number=None, start_residue_number=None, **kwargs ) :
		'''
		returns ATOM records, reference seqres based sequence, and aligned sequence of ATOM records to the SEQRES records.
		
		If biological is true, the records will be biologically relevant,
		most likely to be the ATOM records without vector contaminants.
		If standard residue name is true, the record will be unmodified standard residues
		found in the sequence database instead of the original line.
		If backbone is true, the residues do not have all the backbone residues will be
		removed from the records.
		
		Also the boundaries of the records can be determined by the start residue id
		and the end residue id.
		'''

		atomrecord = True
		start, end = self.residue_id2residue_index(start_residue_id), self.residue_id2residue_index(end_residue_id)
		if end != None : #adjust end index for slicing.
			end = end + 1

		atom_serial_number=1
		residue_serial_number=1

		if start_number != None :
			atom_serial_number = start_number
		if start_residue_number != None :
			residue_serial_number = start_residue_number

		rfp = io.StringIO() #reference sequence fp
		fp = io.StringIO() #sequence fp
		afp = io.StringIO() #atom record fp

		for residue_id in self.residue_ids[start:end] :
			residue = self.get_residue(residue_id)
			if not residue.is_belong_to_chain() :
				continue

			if atomrecord and residue.is_missing_residue():
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue
			
			if biological and not residue.is_biological_residue() :
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue
		
			if backbone and not residue.has_full_backbone_atoms() :
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue

			if standard_residue_name :
				rfp.write( residue.get_standard_amino_acid_name_letter() )
				if residue.has_ca():
					fp.write( residue.get_standard_amino_acid_name_letter() )
					atom_records = residue.get_atom_records( use_standard_aa=True, number=atom_serial_number, resnum=residue_serial_number, **kwargs )
					afp.write( atom_records )
					atom_serial_number += atom_records.count('\n')
					residue_serial_number += 1
				else :
					fp.write( '=' )
			else :	
				rfp.write( residue.get_name_letter() )
				if residue.has_ca() :
					fp.write( residue.get_standard_amino_acid_name_letter() )
					atom_records = residue.get_atom_records( use_standard_aa=False, number=atom_serial_number, resnum=residue_serial_number, **kwargs )
					afp.write( atom_records )
					atom_serial_number += atom_records.count('\n')
					residue_serial_number += 1
				else :
					fp.write( '=' )

		return afp.getvalue(), rfp.getvalue(), fp.getvalue()
		
		
	def get_biological_sequence( self ) :
		'''
		returns a sequence based on SEQRES yet biologically relevant region only.
		'''
		fp = io.StringIO()
		for residue_id in self.residue_ids :
			residue = self.get_residue( residue_id )
			if residue.is_belong_to_chain() and residue.is_biological_residue() :
				fp.write( residue.get_standard_amino_acid_name_letter() )

		return fp.getvalue()
				
	
	def generate_sequence_mapping( self, biological=True, atomrecord=True, backbone=False, 
				standard_residue_name=True, start_residue_id=None, end_residue_id=None ) :
		'''
		returns two sequences.
		The first is the reference is SEQRES sequence.
		The reference is SEQRES sequence since the SEQRES sequence is longest possible.
		And the second is sequence of ATOM records or SEQRES record.
		
		If atomrecord is true, the sequence will be based on ATOM records, or else
		the sequence will be basing on SEQRES record.
		If biological is true, the sequence will be biologically relevant,
		most likely to be the seqeunce without vector contaminants.
		If standard residue name is true, the sequence will be unmodified standard residues
		found in the sequence database.
		If backbone is true, the residues do not have all the backbone residues will be
		removed from the sequence.
		
		Also the boundaries of the sequence can be determined by the start residue id
		and the end residue id.
		'''
		start, end = self.residue_id2residue_index(start_residue_id), self.residue_id2residue_index(end_residue_id)
		if end != None : #adjust end index for slicing.
			end = end + 1

		rfp = io.StringIO()
		fp = io.StringIO()
		for residue_id in self.residue_ids[start:end] :
			residue = self.get_residue(residue_id)
			if not residue.is_belong_to_chain() :
				continue

			if atomrecord and residue.is_missing_residue():
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue
			
			if biological and not residue.is_biological_residue() :
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue
		
			if backbone and not residue.has_full_backbone_atoms() :
				if standard_residue_name :
					rfp.write( residue.get_standard_amino_acid_name_letter() )
					fp.write( '=' )
				else :
					rfp.write( residue.get_name_letter() )
					fp.write( '=' )
				continue

			if standard_residue_name :
				rfp.write( residue.get_standard_amino_acid_name_letter() )
				fp.write( residue.get_standard_amino_acid_name_letter() )
			else :	
				rfp.write( residue.get_name_letter() )
				fp.write( residue.get_name_letter() )

		return rfp.getvalue(), fp.getvalue()
		



	def get_sequence( self, residue_ids=None ) :
		'''
		returns one letter amino acid sequence derived from given residue sequence.
	
		Note that this function is designed to be used for internal use.
		extracting sequences needs to be used more specific functions.

		#get_atom_sequence()
		#get_seqres_seqeuence()

		#get_sequence_mapping()
		'''
		if residue_ids == None :
			residue_ids = self.get_within_chain_residue_ids()

		seqinfo = self.seqinfo
		fp = io.StringIO()
		#start to process the given residue_ids
		for residue_id in residue_ids :
			residue = self.get_residue( residue_id )

			if residue and residue.is_belong_to_chain() :
				fp.write( residue.get_name_letter() )

			elif seqinfo.is_missing_residue( residue_id ) :
				resname = seqinfo.missing_residue_id2resname( residue_id )
				fp.write( amino_acids.get(resname, 'X') )
			else :
				if verbose :
					print("Residue ID is not belong to the chain, nore in missing residue id.", residue_id, file=sys.stderr)

		sequence = fp.getvalue()
		fp.close()
		return sequence
	
	def is_same_sequences( self, seq1, seq2 ) :
		''' 
		returns True if the two sequences are approximately same.

		The approximation here means 
		if the two sequences have same length but matching standard amino acids 
		to the non-standard (appeared to be X) amino acids.
		'''
		if not seq1 and not seq2 :
			return True
		
		if len(seq1) != len(seq2) :
			return False
		for i, (aa1, aa2) in enumerate(zip( seq1, seq2 )) :
			aa1 = aa1.upper()
			aa2 = aa2.upper()
			if aa1 == aa2 or aa1 == 'X' or aa2 == 'X' :
				continue
			else :
				return False
		return True

	def check_good_id_order( self, residue_ids, index, residue_id ) :
		''' 
		returns True if the residue_id is larger than the index-1 and
		smaller than the index position of the residue_ids.

		Designed for my internal usage in mapping SEQRES to ATOM records.
		'''
		if not residue_ids :
			print("WARNING: No elements found in residue_ids!", file=sys.stderr)
			sys.exit()
		if 0 < index < len(residue_ids) :
			if residue_ids[index-1] < residue_id < residue_ids[index] :
				return True
			else :
				return False
		elif index == 0 :
			if residue_id < residue_ids[index]:
				return True
			else :
				return False
		elif index >= len(residue_ids) :
			if residue_ids[-1] < residue_id :
				return True
			else :
				return False

	"""
	Potentially solves all the problems.
	
	def map_seqres_record_without_missing_records( self ) :
		'''sets up the chain info missing records info
		and makes some fake ids for the non-existing 
		'''
		chaininfo = self.get_chaininfo()
	"""
			
	def map_seqres_record( self ) :
		''' 
		maps without changing the order of residues in the chain.
		This is nessary because of the some cases of reversed order 
		in the residue numberings.
		
		Alse the mapping will be saved into the chain instance.
		'''
		if verbose :
			print("mapping of seqres record is starting...")
			print("chainid:", self.chainid)
			print("molname:", self.molname)
	
		#######################################
		#Preliminary checking for the mapping
		#######################################

		#if already mapped
		#do not need to do it again!
		if self.is_missing_residue_id_mapped() :
			#return self.get_all_residue_ids()
			if verbose :
				print("WARNING: missing residue mapping is already done!")
			return

		#getting primary sequence information
		seqinfo = self.seqinfo

		#if no primary sequence information is available
		#stop it!
		if not seqinfo :
			if verbose :
				print("WARNING: SequenceInfo object is not ready!")
			return

		#If no missing residues exist,
		#we do not need to map.
		#And set the mapping flag success!
		if not seqinfo.has_missing_residues() :
			if verbose :
				print("WARNING: Primary Structure indicates that there are no missing residues!")
			self.seqres_atom_mapping_flag = 1
			return
		#######################################

		########################################
		#Preparing to map SEQRES to ATOM records
		########################################

		#get list of residues in the order appearing in the PDB file
		residue_ids = self.get_within_chain_residue_ids()

		#make copy of residues for mapping
		#to use it as a test case!
		mapping_ids = residue_ids[:]
		
		####################################################
		#using the following stuff to fill the missing spots
		#seqinfo.missing_resnames
		#seqinfo.missing_residues
		#seqinfo.missing_residue_models
		####################################################
		#build missing residue id list 
		model_num = self.model.model_id
		missing_ids = []

		#########################
		#check if the missing_id
		#########################

		#####################
		#this block  changes the missing residue id
		#if they are overlapping...
		#####################
		for rid in seqinfo.get_missing_residue_ids( model_num ) :
			if rid in mapping_ids :
				if verbose :
					print("WARNING: Redundency in missing IDs found.", rid, file=sys.stderr)

				#add non clashing id?
				'''
				old_rid = rid
				while rid in mapping_ids :
					rid = self._decrease_icode( rid )
				
				seqinfo.change_missing_residue_id( old_rid, rid )
				missing_ids.append( rid )
				'''
				#or remove it?
				seqinfo.remove_missing_residue_id( rid, model_num )
				pass

			else :
				missing_ids.append( rid )


		#for missing_residue_id, model_number_list in zip(seqinfo.missing_residues, seqinfo.missing_residue_models) :
			#if model_num in model_number_list :
				#missing_ids.append( missing_residue_id )

		if verbose :
			print("seqinfo.missing_resdiues:", seqinfo.missing_residues)
			print("seqinfo.missing_resdiue_models:", seqinfo.missing_residue_models)
			print("missing_ids", missing_ids) 

		#If no missing IDs.
		#if not missing_ids :
			#print >>sys.stderr, "WARNING! The record indicates that missing residues exist, but building list of missing residue has failed!", self.molname, self.chainid, self.model.pdb.get_pdb_id()
			#return 

		if not mapping_ids :
			if verbose :
				print("WARNING: No mapping_ids are available.")
			return

		seqres_sequence = self.get_seqres_sequence()
		#Preparation ends.
		########################################

		########################################
		#simpler yet more error prone method is to 
		#add all residue ids into one list 
		#and sort them all.
		########################################
		sorted_mapping_ids = mapping_ids[:]
		sorted_missing_ids = missing_ids

		if is_sorted( sorted_mapping_ids ) and is_sorted(sorted_missing_ids) :
			sorted_mapping_ids.extend( sorted_missing_ids )
			
			sorted_mapping_ids.sort()
					
			if len(seqres_sequence) == len(sorted_mapping_ids) :
				if verbose :
					print("Mapping length is good", file=sys.stderr)
				pass #good

			#elif len(seqres_sequence) > len( sorted_mapping_ids) :
			else :
				#Some Insertions are missed!
				#try to recover by putting the missing missing infos in.

				#for i, (a,b,c) in enumerate( zip( seqres_sequence, self.get_sequence(residue_ids=sorted_mapping_ids), sorted_mapping_ids) ) :
					#print i, a, b 
				#print seqres_sequence[ len(sorted_mapping_ids): ]

				"""
				#old code to map roughly..
				#assuming the easy mapping.
				#but seems to be not good!
				self._recover_missing_missing_residues( sorted_mapping_ids )

				#reset the arrays
				sorted_mapping_ids = mapping_ids[:]

				missing_ids = seqinfo.get_missing_residue_ids(model_num)
				sorted_missing_ids = missing_ids[:]

				if is_sorted(sorted_mapping_ids) and is_sorted(sorted_missing_ids) :
					if verbose :
						print >>sys.stderr, "Second trial of simple ordering method"

					sorted_mapping_ids.extend( sorted_missing_ids )
					sorted_mapping_ids.sort()
				"""
				#new code run global alignment
				#to bypass the problem in 1e4f 
				#this part is tried after one-by-one

				if verbose :
					print("Lengths of mapping ids + missing and seqres record is not same...", file=sys.stderr)
					print("Global alignment procedure should be used!", file=sys.stderr)

				mapping_ids = self.map_seqres_record_by_global_alignment() 

			'''
			else :
				#if verbose :
				self.seqinfo.print_seqres_information()
				print >>sys.stderr, "WARNING: SEQRES sequence is shorter than ATOM + MISSING."
				print >>sys.stderr, "len(seqres_sequence)", len(seqres_sequence)
				print >>sys.stderr, "len(sorted_mapping_ids):", len(sorted_mapping_ids)
				print >>sys.stderr, seqres_sequence
				print >>sys.stderr, self.get_sequence( residue_ids = sorted_mapping_ids )

				raise SEQRESMappingError( "SEQRES sequence is shorter than ATOM + MISSING.", len(seqres_sequence), len(sorted_mapping_ids) )
			'''
			
		#Test if the simple method works or not..
		if self.is_same_sequences( seqres_sequence, 
					  self.get_sequence( residue_ids = sorted_mapping_ids ) ) :
			#simple ordering works!!
			if verbose :
				print("simple mapping succeeded.")

			mapping_ids = sorted_mapping_ids
		else :
			if verbose :
				print("simple mapping failed")
				print("seqres    :", self.get_seqres_sequence())
				print("simple_map:", self.get_sequence(residue_ids=sorted_mapping_ids ))
				for i, (a, b) in enumerate( zip( self.get_seqres_sequence(), self.get_sequence(residue_ids=sorted_mapping_ids ) )) :
					print(i, a, b, sorted_mapping_ids[i])
					if a != b :
						print("Not same!")
						break

			#######################################
			#Simple ordering is not applicable nor failed..
			#due to chains with the residue numbering is not monotonically increasing... 
			#Now the more robust yet complex adding one by one is tried.
			########################################
			self.map_seqres_record_by_adding_one_by_one( mapping_ids, missing_ids )
		
		#################################################
		#Final check part!!!
		#if seqres sequence is different from the mapping result,
		#and if somehow the SEQRES record does not exists in the 
		#PDB file..
		#I will perfom the global alignment with affine gap penalty written for this mapping purpose.
		#This algorithm uses abnormally high mismatch value (so almost always insert gaps) and fairly low gap penalities.
		#Then based on the mapping alignment, I will generate residue ids and other stuff in the chaininfo class instance.
		#################################################
		if self.is_same_sequences( self.get_seqres_sequence(), 
					  self.get_sequence( residue_ids = mapping_ids ) ) :
			self.seqres_atom_mapping_flag = 1

			if verbose :
				print("Comparing the two sequences is successful at one_by_one method")
				print(self.get_seqres_sequence())
				print(self.get_sequence( residue_ids = mapping_ids ))
		else:  
			#chaininfo.make_missing_residue_information is used to build 
			#information into the ChainInfo class instance in the self.map_seqres_record_by_global_alignment..
			if verbose :
				print("Comparing the two sequences is failed at one_by_one method")
				print(self.get_seqres_sequence())
				print(self.get_sequence( residue_ids = mapping_ids ))
				print("Global alignment method will be tried.")

			mapping_ids = self.map_seqres_record_by_global_alignment() 

		#mapping is done!!
		#################################################
		
		#################################################
		#Adjusting internal structure of the class
		#################################################
		if mapping_ids :
			#change new residue id list
			#that have missing residues

			#rebuild the information..
			#due to fixes for wrong PDBs.
			#the information should be rebuilt.
			model_num = self.model.model_id
			missing_ids = seqinfo.get_missing_residue_ids( model_num )
			missing_residues = {}
			for rid in missing_ids :
				residue = Residue( id=rid, chain=self, name=seqinfo.missing_residue_id2resname(rid), number=rid[1] )
				missing_residues[rid] = residue
		
			#############################
			#Update residue ids list
			old_residue_ids = self.residue_ids
			self.residue_ids = mapping_ids[:] #making copy
			#since just replacing residue_ids to mapping_ids
			#will remove non-polymer residues like ligands, etc
			#add them seperately.
			[self.residue_ids.append(rid) for rid in old_residue_ids if rid not in mapping_ids] 

			#checking if the missing residues 
			#are overlapping with the current set of 
			#residues (without the missing residues).
			intersection = set(missing_residues).intersection( set(self.residues) )
			if not intersection :
				self.residues.update( missing_residues )

			#to check further to correct cases 
			#of the non protein atoms have same residue ID as in amino acid. :(
			#ex) 2jw6
			elif all( [(not self.residues[residue_id].is_belong_to_chain()) for residue_id in intersection ] ) :
				self.residues.update( missing_residues )
				
			else :
				#if the intersection exists, there must be some problem!!
				raise SEQRESMappingError( "Missing residues already in the residues.", 
					str(intersection) )
		else :
			#if no element in mapping_ids,
			#it is a big problem!
			raise SEQRESMappingError( "Final mapping failed!" )
			
		#Old code
		#self.set_all_residue_ids( mapping_ids )
		##################################
		#mapping SEQRES index to ATOM index
		#for easy conversion between the two numbers
		##################################
		#self.seqres_index2atom_index = []
		#rcount = 0
		#for residue_id in mapping_ids :
			#if self.is_residue_within_chain( residue_id ) :
				#rcount += 1
				#self.seqres_index2atom_index.append( rcount )
			#else :
				#self.seqres_index2atom_index.append( rcount )
				
		#return mapping_ids	

	def _increase_icode( self, residue_id ) :
		'''
		returns a new residue id with increased icode in the alphabetical order.
		'''
		max_ascii = 90 #Z

		chainid, resnum, icode = residue_id
		if icode == ' ' :
			return build_residue_id( chainid, resnum, 'A' )

		elif icode.isupper() :
			if ord(icode)+1 > max_ascii :
				return build_residue_id(chainid, resnum+1, ' ')
			return build_residue_id( chainid, resnum, chr(ord(icode)+1) )
		else :
			raise SEQRESMappingError( "This error should not be seen!" )

	def _decrease_icode( self, residue_id ) :
		'''
		returns a new residue id with decreased icode in the alphabetical order.
		'''
		min_ascii = 65 #A

		chainid, resnum, icode = residue_id
		if icode == ' ' :
			return build_residue_id( chainid, resnum-1, 'Z' )

		elif icode == 'A' :
			return build_residue_id( chainid, resnum, ' ' )

		elif icode.isupper() :
			return build_residue_id( chainid, resnum, chr(ord(icode)-1) )
		else :
			raise SEQRESMappingError( "This error should not be seen!" )



	def _populate_middle_residue_ids( self, before_id=None, after_id=None, nmiddle=0, missing_ids=[] ) :
		'''
		Builds middle residues ids bounded by the given ids.
		This function is intended to be used in _recover_misssing_missing_residues.
		'''

		#if nothing is given, nothing can be done!
		if before_id == None and after_id == None :
			return  []

		#if before and after is same I cannot build anything.
		if before_id == after_id :
			return []

		residue_ids = []
		if not before_id and after_id :
			new_rid = after_id

			for i in range( 1, nmiddle+1 ) :
				new_rid = self._decrease_icode( new_rid )
				while new_rid in self.residues or new_rid in missing_ids:
					new_rid = self._decrease_icode( new_rid )

				residue_ids.insert( 0, new_rid )

		elif before_id and not after_id :
			new_rid = before_id
			for i in range( 1, nmiddle+1 ) :
				new_rid = self._increase_icode( new_rid ) 
				while new_rid in self.residues or new_rid in missing_ids:
					new_rid = self._increase_icode( new_rid )
			
				residue_ids.append( new_rid )
		else :
			chain1, resnum1, icode1 = before_id
			chain2, resnum2, icode2 = after_id

			if chain1 == chain2 :
				pass
			else :
				raise SEQRESMappingError( "Chain IDs are different in before and after the missing_residues.", chain1, chain2 ) 

			if resnum2 - resnum1 -1 == nmiddle :
				for i in range( resnum1+1, resnum2 ) :
					new_rid = build_residue_id( chain1, i, ' ' )
					while new_rid in self.residues or new_rid in missing_ids :
						new_rid = self._increase_icode( new_rid )

					residue_ids.append( new_rid )
			else :
				new_rid = before_id
				for i in range( 1, nmiddle+1 ) :
					new_rid = self._increase_icode( new_rid ) 
					while new_rid in self.residues or new_rid in missing_ids:
						new_rid = self._increase_icode( new_rid )
				
					residue_ids.append( new_rid )

		if len(residue_ids) == nmiddle :
			return residue_ids
		else :
			raise SEQRESMappingError( "Populating residue ids for missing missing failed!", str(residue_ids) )
		

	def _determine_insert_point( self, sorted_mapping_ids, sorted_missing_ids, mapping_index, after=False ) :
		'''
		This function finds the index sorted_missing_ids 
		corresponding to the mapping_index 
		(or the closest in sorted_mapping_ids).

		This function is designed to find correct insertion point in
		the recovering missing missing residue id.
		'''
		model_num = self.model.model_id
		#convert the mapping index 
		#to actural residue_id
		if mapping_index < len( sorted_mapping_ids ) :
			residue_id = sorted_mapping_ids[ mapping_index ]
		else :
			mapping_index = len( sorted_mapping_ids )-1
			residue_id = sorted_mapping_ids[-1]


		#easy case where the residue_id is in
		#sorted_missing_ids.. meaning the 
		#residue id is a missing resdiue.
		#if residue_id in sorted_missing_ids :
		if self.seqinfo.exists_missing_residue( residue_id, model_num ) :
			if after :
				return self.seqinfo.find_missing_residue_index( residue_id, model_num ) + 1
			else :
				return self.seqinfo.find_missing_residue_index( residue_id, model_num )

		#things are not that easy.
		#the residue should be found
		#through the proxy sorted_mapping_ids.

		#start from back toward front
		while mapping_index >= 0 :
			if self.seqinfo.exists_missing_residue( residue_id, model_num ) :
			#if residue_id in sorted_missing_ids :
				#if after :
					#return sorted_missing_ids.index( residue_id )+1
				#else :
					#return sorted_missing_ids.index( residue_id )
				#maybe all cases should be like this..
				#return sorted_missing_ids.index( residue_id ) + 1
				return self.seqinfo.find_missing_residue_index( residue_id, model_num ) + 1

			mapping_index -= 1
			residue_id = sorted_mapping_ids[mapping_index]
		else :
			return 0
			
				
	def _recover_missing_missing_residues( self, mapping_ids, seqres_aln, atomres_aln, missing_sequence, insert=True ) :
		'''
		This function will try to recover missing information 
		about missing residues.

		If insert is True,
		recovery by adding lost missing residue info.
		If insert is False,
		recovery by removing wrong missing residue info.
		'''
		if verbose :
			print("WARNING: Trying to recover missing missing residue information.", file=sys.stderr)

		model_num = self.model.model_id
		seqres_sequence = self.get_seqres_sequence()
		mapped_sequence = self.get_sequence( residue_ids=mapping_ids )

		#initial conditioning.
		i,j = -1,-1 #seqres_sequence & mapped_sequence indices
		x = -1 #missing sequence index

		n,m = len(seqres_sequence), len(mapped_sequence)
		z = len( missing_sequence )

		diff = n-m-z
		missing_block = 0
		before_missing_id = None

		for h, (a, b) in enumerate(zip(seqres_aln, atomres_aln)) :
			if a.isalpha() : 
				i += 1

			c = '' #missing residue

			#b is either atomres or missing
			if b.isalpha() : 
				j += 1 #atomres
			else :
				x += 1 #missing 
				#a, b are seqres and atom respectively.
				if x < z :
					c = missing_sequence[x]#missing residue letter
				else :
					c = ' '
			
			if a == b :
				#residue atom and seqres is matching
				continue

			elif a == c :
				#the match is already in missing residue set.
				continue
			else :
					
				if insert :
					insert_point = x
					fillin_seq = a
					missing_ids = self.seqinfo.get_missing_residue_ids( model_num )
					fillin_residue_ids = self._populate_middle_residue_ids( after_id=mapping_ids[max(j,0)], nmiddle=1, missing_ids=missing_ids )

					if verbose :
						print("inserting missing residue information", insert_point, file=sys.stderr)
						print(fillin_residue_ids, fillin_seq, file=sys.stderr)
						print("a, b, c", [a, b, c], file=sys.stderr)
						print("i, j, x", [i, j, x], file=sys.stderr)
						print("aln index", h, file=sys.stderr)
						print(seqres_aln, file=sys.stderr)
						print(atomres_aln, file=sys.stderr)


					self.seqinfo.insert_missing_residue_information( fillin_residue_ids, fillin_seq, model_num, insert_point )

					if verbose :
						self.seqinfo.print_seqres_information()

					diff -= 1
					if diff < 0 :
						if verbose :
							print('seqres_aln ', seqres_aln, file=sys.stderr)
							print('atomres_aln', atomres_aln, file=sys.stderr)
							self.seqinfo.print_seqres_information()
							
						raise SEQRESMappingError( "Recovering missing missing residue failed due to incompatible sequences." )

					missing_sequence = resnames2sequence( self.seqinfo.get_missing_resnames( self.model.get_model_id() ) )
					z = len( missing_sequence )
				else :
					missing_ids = self.seqinfo.get_missing_residue_ids( model_num )

					if verbose :
						print("removing missing residue_information", file=sys.stderr)
						print("a, b, c", [a, b, c], file=sys.stderr)
						print("i, j, x", [i, j, x], file=sys.stderr)
						print("aln index", h, file=sys.stderr)
						print(seqres_aln, file=sys.stderr)
						print(atomres_aln, file=sys.stderr)
					
					self.seqinfo.remove_missing_residue_id( missing_ids[x], model_num )
		
					if verbose :
						self.seqinfo.print_seqres_information()

					diff += 1
					if diff > 0 :
						raise SEQRESMappingError( "Recoving excessive missing residue failed due to incompatible sequences." )
					missing_sequence = resnames2sequence( self.seqinfo.get_missing_resnames( self.model.get_model_id() ) )
					z = len( missing_sequence )
						
						

		#####################
		#final data consistency checking.
		if diff != 0 :
			raise SEQRESMappingError( "After the recovering the diff should be 0.", diff)

		if i == n-1 and j == m-1 and x == z-1 :
			pass
		else :
			raise SEQRESMappingError( "After the recovering all i,j,x values should be same as n,m,z.", str( [i,j,x])+str([n,m,z]) )
		#####################

		if verbose :
			print("Recovering finsiehd!", file=sys.stderr)
		#everything is perfectly recovered!!
		return True

	#############################3
	# This function is deprecated!
	# The idea of adding missing
	# info is good but the method
	# detemine how to add is too crude.
	##############################
	def _recover_missing_missing_residues_old( self, sorted_mapping_ids ) :
		'''
		This function will try to recover missing information 
		about missing residues.

		Note that the function was originally developed for
		recovery missing residue information from sorted mapping ids.
		But this function does not requite the "sorted" mapping ids.
		'''
		print("WARNING: Trying to recover missing missing residue information.", file=sys.stderr)
		seqres_sequence = self.get_seqres_sequence()
		mapped_sequence = self.get_sequence( residue_ids=sorted_mapping_ids )

		#initial conditioning.
		i,j = 0,0
		n,m = len(seqres_sequence), len(mapped_sequence)
		diff = n-m
		missing_block = 0
		before_missing_id = None

		while i < n and j < m :
			if i < n :
				a = seqres_sequence[i]
			else :
				a = ' '
		
			if j < m :
				b = mapped_sequence[j]
			else :
				SEQRESMappingError( "Inserting failed, j should be lower than m always.", str([i,j,m,n])  )

			if a == b :
				if missing_block :
					#this means missing block is ended.
					after_missing_index = i
					#populates residue ids for fill in missing.
					model_num = self.model.model_id
					sorted_missing_ids = self.seqinfo.get_missing_residue_ids( model_num )
				
					if j >= m :
						fillin_residue_ids = self._populate_middle_residue_ids( before_id=before_missing_id, nmiddle=missing_block )
					else :
						fillin_residue_ids = self._populate_middle_residue_ids( before_id=before_missing_id, after_id=sorted_mapping_ids[j], nmiddle=missing_block )
					fillin_seq = seqres_sequence[start_missing_index:after_missing_index]
					insert_point = self._determine_insert_point( sorted_mapping_ids, sorted_missing_ids, j, after=False )


					if verbose :
						print("inserting missing residue information", file=sys.stderr)
						print(before_missing_id, fillin_residue_ids, sorted_mapping_ids[j], missing_block, file=sys.stderr)
						print(fillin_residue_ids, fillin_seq, file=sys.stderr)
					self.seqinfo.insert_missing_residue_information( fillin_residue_ids, fillin_seq, self.model.model_id, insert_point )
					missing_block = 0
						
				if j >= m :
					before_missing_id = sorted_mapping_ids[-1]
				else :
					before_missing_id = sorted_mapping_ids[j]
				
				i, j = i+1, j+1
			else :
				diff -= 1
				if missing_block :
					missing_block += 1
				else :
					missing_block = 1
					start_missing_index = i

				if diff < 0 :
					raise SEQRESMappingError( "Recovering missing missing residue failed due to incompatible sequences." )
				
				i = i+1

		if verbose :
			print("Recovering finsiehd!", file=sys.stderr)
		#everything is perfectly recovered!!
		return True

	def map_seqres_record_by_adding_one_by_one( self, mapping_ids, missing_ids ) :
		'''
		Generally good heuristic solution to map sequence to structure.
		It is developed for fast mapping yet to handle those midly erroneous cases of 
		mapping sequences to structures.
		'''
		#id_count = 0 #total number of id's mapped
		id_index = 0 #index of the id so far checked
		
		while missing_ids :
			###############
			# First par of missing id mapping
			# get possible insertion points
			###############
			#use max_id and min_id for reducing the ambiguity.. 
			#but maybe it is not right idea.. (PDBID:1IAO showed that sometimes they can be really messed up!!)
			#recalculate max and min for each round since missing ids might change them
			#max_id = max( mapping_ids )
			#min_id = min( mapping_ids )
			mid = missing_ids.pop(0)
			possible_insertion_points = []
			if id_index < len(mapping_ids) :
				if mid < mapping_ids[id_index] :
					possible_insertion_points.append(id_index)
				else :
					if id_index > 0 and mapping_ids[id_index] < mapping_ids[id_index-1] :
						if mid > mapping_ids[id_index] or mid < mapping_ids[index-1] :
							possible_insertion_points.append( id_index )
			for i, (rid1, rid2) in enumerate( zip( mapping_ids[id_index:-1], mapping_ids[id_index+1:] ) ) :
				i = id_index + i
				if rid1 < mid < rid2 :
					possible_insertion_points.append( i+1 )
					continue
				elif rid1 > rid2 : #when the continuity is violated!
					#if rid1 == max_id and mid > rid1 :
					if mid > rid1 :
						possible_insertion_points.append( i+1 )
						continue
					#elif rid2 == min_id and mid < rid2 :
					elif mid < rid2 :
						possible_insertion_points.append( i+1 )
						continue
			else :
				possible_insertion_points.append( len(mapping_ids) )
			if verbose :
				print('mid:', mid, "id_index:", id_index)
				print("possible_insertion_points:", possible_insertion_points)

			###############
			# Second part of missing id mapping
			# Analyse the possible insertion points
			###############
			#easy case for the problem!!
			if len( possible_insertion_points ) == 1 :
				insertion_point = possible_insertion_points[0]
				mapping_ids.insert( insertion_point, mid )
				#check for the insertion validity by checking the sequences
				if self.is_same_sequences(  #this is function relaxed than brutal comparison
					self.get_seqres_sequence(  start=id_index, end=insertion_point+1 ),
					self.get_sequence( residue_ids =  mapping_ids[id_index: insertion_point+1] ) ) :
					#this is expected!!
					pass
				else :
					if verbose :
						print('WARNING:, Mapping at current insertion %d is not right!' % insertion_point, file=sys.stderr)
						print('WARNING1', self.get_seqres_sequence( start=id_index, end=insertion_point+1 ), file=sys.stderr)
						print('WARNING2', self.get_sequence( residue_ids =  mapping_ids[id_index: insertion_point+1] ), file=sys.stderr)
				#need to adjust 1 bigger than the insertion point to move forward
				id_index = insertion_point + 1
				continue
					
			else :
				#not easy to sort out the cases!!!
				fipoints = [idx for idx in possible_insertion_points if self.is_same_sequences(
								self.get_seqres_sequence( start=id_index, end=idx ),
								self.get_sequence( residue_ids= mapping_ids[id_index: idx] ) ) 
							and 
							self.is_same_sequences( 
								self.get_seqres_sequence( start=idx, end=idx+1 ),
								self.get_sequence( residue_ids=[ mid ] ) )]
				if verbose :
					print("##############################################")
					print("multiple possible insertion points were found")
                                	print(self.get_seqres_sequence( start=possible_insertion_points[-1], end=possible_insertion_points[-1]+1 ))
                                	print(mid, self.get_sequence( residue_ids=[ mid ] ))
					for i in possible_insertion_points :
						print("##### Insertion position:", i, "start:", 0, "end:", id_index, "#####", file=sys.stderr)
						print(self.get_seqres_sequence( end=id_index ), file=sys.stderr)
						print(self.get_sequence( residue_ids= mapping_ids[:id_index] ), file=sys.stderr)
						print("##### Insertion position:", i, "start:", id_index, "end:", i, "#####", file=sys.stderr)
						print(self.get_seqres_sequence( start=id_index, end=i ), file=sys.stderr) 
						print(self.get_sequence( residue_ids= mapping_ids[id_index:i] ), file=sys.stderr)
						print('SEQRES[%d]:'%i, self.get_seqres_sequence( start=i, end=i+1 ), file=sys.stderr)
						print('Candidate :', self.get_sequence( residue_ids=[ mid ] ), file=sys.stderr)
						
					print("possible_insertion_points:", possible_insertion_points)
					print("After first check...")
					print("fipoints:", fipoints)
				
				if len(fipoints) == 1 :
					insertion_point = fipoints[0]
					#do not need to further check already checked during filtering process
				elif len(fipoints) > 1 :
					#the logest insertion point selected 
					#if ambiguoty exists
					#because the longest one is least ambiguous!!!
					
					#but I also found problem in this thought!!
					#longest one might not be right especially the first first is randomly matching with
					#the missing residue identity!!!
					
					if verbose :
						print("fipoints:", fipoints, file=sys.stderr)
						for i in fipoints :
							print("##### Insertion position:", i, "start:", 0, "end:", id_index, "#####", file=sys.stderr)
							print("", self.get_seqres_sequence( end=id_index ), file=sys.stderr)
							print("", self.get_sequence( residue_ids= mapping_ids[:id_index] ), file=sys.stderr)
							print("##### Insertion position:", i, "start:", id_index, "end:", i, "#####", file=sys.stderr)
							print("", self.get_seqres_sequence( start=id_index, end=i ), file=sys.stderr) 
							print("", self.get_sequence( residue_ids= mapping_ids[id_index:i] ), file=sys.stderr)
						
					#too loose 
					fipoints2 = [i for i in fipoints if self.check_good_id_order( mapping_ids, i, mid )]
					if verbose :
						print("AFTER second filtering\nfipoints2:", fipoints2, file=sys.stderr)
						for i in fipoints2 :
							print("##### Insertion position:", i, "start:", 0, "end:", id_index, "#####", file=sys.stderr)
							print("", self.get_seqres_sequence( end=id_index ), file=sys.stderr)
							print("", self.get_sequence( residue_ids= mapping_ids[:id_index] ), file=sys.stderr)
							print("##### Insertion position:", i, "start:", id_index, "end:", i, "#####", file=sys.stderr)
							print("", self.get_seqres_sequence( start=id_index, end=i ), file=sys.stderr) 
							print("", self.get_sequence( residue_ids= mapping_ids[id_index:i] ), file=sys.stderr)
						
					#possibly the biggest index in the second good residue ordering filter might be good
					if fipoints2 :
						insertion_point = fipoints2[-1] 
						if verbose :
							print("insertion point determined!", insertion_point, file=sys.stderr)
					else :
						if verbose :
						#in case the filterring is too much!!
							print("WARNING! fipoints are too ambiguous..", fipoints, file=sys.stderr)
						insertion_point = fipoints[-1] 
				else : #No potential insertion point detected
					if verbose :
						print("WARNING! No possible insertion point was detected!", fipoints, possible_insertion_points, mid,missing_ids, file=sys.stderr)
						print('mid:',  mid, 'missing_ids:', missing_ids, file=sys.stderr)
					#printint out error message
					#for i in possible_insertion_points :
						#print >>sys.stderr, "##### Insertion position:", i, "start:", 0, "end:", id_index, "#####"
						#print >>sys.stderr, "Error1", self.get_seqres_sequence( end=id_index )
						#print >>sys.stderr, "Error2", self.get_sequence( residue_ids= mapping_ids[:id_index] )
						#print >>sys.stderr, "##### Insertion position:", i, "start:", id_index, "end:", i, "#####"
						#print >>sys.stderr, "Error1", self.get_seqres_sequence( start=id_index, end=i )
						#print >>sys.stderr, "Error2", self.get_sequence( residue_ids= mapping_ids[id_index:i] )
						#print >>sys.stderr, chaininfo.missing_resnames
						#print >>sys.stderr, chaininfo.missing_residues
						#print >>sys.stderr, chaininfo.missing_residue_models
					#since the mapping through MISSING record did not work,
					#just break the loop.
					#alignment through generic alignment procedure will be done

					#reverting is not necessary.
					#since this part is taken into a separate function.
					#revert all intermediate mapping...
					#mapping_ids = residue_ids[:]
					if verbose: print("Could not finish mapping without generic alignment procedure!", file=sys.stderr)
					if verbose: print("Need to use generic alignment procedure!", file=sys.stderr)
					self.seqres_atom_mapping_flag = -2 #first mark for problem
					break
						
					
				if verbose :
					print("inserting", insertion_point, mid)
				mapping_ids.insert( insertion_point, mid ) 
				id_index = insertion_point + 1
				if verbose :
					print(mapping_ids[insertion_point], file=sys.stderr)
				continue

		
	def map_seqres_record_by_global_alignment( self, loose_mapping=True, fix_redundent_residue_id=True ) :
		'''
		This function use function "align_sequence_with_affine_gap" which is an implementation of
		Goto's global alignment dynamic programming algorithm with one special modification.
		We use identity matrix with high penalty of mismatch and low gap penalty. 
		So my implementation almost always do not misalign SEQRES sequences and ATOM sequences but just
		putting gaps. (That is right in many cases).

		Loose mapping flag is defined to bypass MISSING record problems.
		So far PDBID|1USM found.
		'''
		######################################
		# SEQRES sequence and ATOM residue sequence is not same!
		# We need to recover the 
		########################################
		seqinfo = self.seqinfo
		#sequence extracted from SEQRES
		seqres_seq = self.get_seqres_sequence()
		#this is the main residue ID list to be modified!!
		residue_ids = self.get_within_chain_residue_ids()
		mapping_ids = residue_ids[:] 
		atomres_seq = self.get_sequence( residue_ids = residue_ids )
		if verbose :
			print('mapping by global alignment starting!', file=sys.stderr)
			print(seqres_seq, file=sys.stderr)
			print(atomres_seq, file=sys.stderr)

		#simple checking if the two sequences are same..
		if seqres_seq == atomres_seq :
			#Everything is good!!
			return residue_ids
			
		#another simple checking 
		#if the seqres sequence contains atomres seq
		if atomres_seq in seqres_seq :
			index = seqres_seq.index( atomres_seq )
			seqres_aln = seqres_seq
			atomres_aln= '-'*index + atomres_seq + '-'*( len(seqres_seq)-len(atomres_seq)-index )
		#if simple mapping does not wor,
		#go for the more complex alignment mapping!!
		else :
			#/evolvable_database_source_code/evdblib/Utils/Algorithms/SequenceAlignment.py
			from evdblib.Utils.Algorithms.SequenceAlignment import align_sequences_with_affine_gap_position_specific_gap_opening_penalty
			seqres_aln, atomres_aln, alnscore = align_sequences_with_affine_gap_position_specific_gap_opening_penalty( self.get_seqres_sequence(), self.get_sequence( mapping_ids ), chain_continuous=self.get_chain_continuous() )
			
		#if gap is present in seqres_aln
		#it is a problem!!
		if '-' in seqres_aln  :
			if verbose :
				print("Error: gap is present in seqres alignment line!", file=sys.stderr) 
				print("Error1", seqres_aln, file=sys.stderr)
				print("Error2", atomres_aln, file=sys.stderr)
			raise SEQRESMappingError( "Error: ATOM + MISSING has more residues than SEQRES." )
				
	
		if len(seqres_aln) == len(atomres_aln) and '-' in atomres_aln :
			#extract list of insertioned in seqres_aln
			#by comparing to atomres_aln
			seqres_insertions = ''.join ( [ p[0] for p in [p for p in zip( seqres_aln, atomres_aln ) if p[1]=='-'] ] )
			missing_sequence = resnames2sequence( seqinfo.get_missing_resnames( self.model.get_model_id() ) )
			insertion_points = [ p[0] for p in [p for p in enumerate(atomres_aln) if p[1] == '-'] ]

			if verbose :
				print(insertion_points)
				print(" seqres_aln:", seqres_aln, file=sys.stderr)
				print("atomres_aln:", atomres_aln, file=sys.stderr)

			#no missing residue definition presnt!!
			if not missing_sequence :
				#to save into chaininfo
				missing_residue_ids = [] 
				largest_residue_id = max( mapping_ids )
				while insertion_points :
					insertion_point = insertion_points.pop(0) #FIFO
					largest_residue_id = build_one_bigger_residue_id( largest_residue_id )
		
					mapping_ids.insert( insertion_point, largest_residue_id )
					missing_residue_ids.append( largest_residue_id )

				#this function build information into the ChainInfo instace
				#so that later we can use the mapping info correctly.
				making_missing_residue_result = seqinfo.make_missing_residue_information( missing_residue_ids, seqres_insertions, self.model.model_id )
				if verbose :
					print(seqres_insertions)
					print(seqinfo.missing_residues)
					print(seqinfo.missing_resnames)

			#Check if the mapping is correct or not
			elif seqres_insertions == missing_sequence :
				missing_residue_ids = seqinfo.get_missing_residue_ids( self.model.model_id )[:] #copy of the 
				for insertion_point, residue_id in zip( insertion_points, missing_residue_ids) :
					#checking if the redundency between missing vs atom records.
					changed = 0
					if residue_id in mapping_ids :
						changed = 1
						old_residue_id = residue_id
						while residue_id in mapping_ids or residue_id in missing_residue_ids :
							residue_id = self._decrease_icode( residue_id )

					if changed :
						#redundency found!
						if fix_redundent_residue_id :
							seqinfo.change_missing_residue_id( old_residue_id, residue_id )
						else :
							raise SEQRESMappingError( "A Redundent Residue ID found!", old_residue_id )
						
					mapping_ids.insert( insertion_point, residue_id )

			#loose mapping is applied only if
			#the loose mapping flag is true
			#and the length of the mapping sequence is same.
			elif loose_mapping and len(seqres_insertions) ==  len(missing_sequence) :
				print("WARNING: Loose mapping (just checking the length of insertion) is applied to bypass the problems like 1usm incorrect missing residue name.", file=sys.stderr)
				if verbose :
					print("seqres_insertions:", seqres_insertions)
					print("missing_sequence :", missing_sequence)

				#checking eligibility of loose mapping bypass
				count_differences = 0
				max_differences = 20
				for i, (a, b) in enumerate(zip( seqres_insertions, missing_sequence )) :
					if a == b :
						continue

					count_differences += 1
					
					if count_differences > max_differences :
						print("Error failed to recover: More than %d residue is wrong at missing residue section" % max_differeneces, file=sys.stderr)
						break

					#need to fix the possibly wrong record in missing_resname in seqinfo.
					seqinfo.missing_resnames[i] = amino_acids_rev[a]


				missing_residue_ids = seqinfo.get_missing_residue_ids( self.model.model_id)[:]
				while insertion_points :
					insertion_point = insertion_points.pop(0)
					residue_id = missing_residue_ids.pop(0)
					mapping_ids.insert( insertion_point, residue_id )

			#If the missing residue information is missing..
			elif len(seqres_insertions) > len(missing_sequence) :
				if verbose :
					print("MISSING residue information is missing!")
					print('seqres_aln ', seqres_aln, file=sys.stderr)
					print('atomres_aln', atomres_aln, file=sys.stderr)
					print("seqres_insertions:", seqres_insertions, file=sys.stderr)
					print("missing_sequence :", missing_sequence, file=sys.stderr)
					for i ,(a, b) in enumerate( zip(seqres_aln, atomres_aln)) :
						print(i, a, b, file=sys.stderr) 

				#recovering by add lost missing residue info.
				self._recover_missing_missing_residues(  mapping_ids, seqres_aln, atomres_aln, missing_sequence )

				missing_residue_ids = seqinfo.get_missing_residue_ids( self.model.model_id)[:]
				while insertion_points :
					insertion_point = insertion_points.pop(0)
					residue_id = missing_residue_ids.pop(0)
					mapping_ids.insert( insertion_point, residue_id )


			elif len(seqres_insertions) < len(missing_sequence) :
				if verbose :
					print("MISSING Reisude information has surplus.", file=sys.stderr) 
					print('seqres_aln ', seqres_aln, file=sys.stderr)
					print('atomres_aln', atomres_aln, file=sys.stderr)
					print('seqres  insertions:', seqres_insertions, file=sys.stderr)
					print('missing_sequence  :', missing_sequence, file=sys.stderr)
					seqinfo.print_seqres_information()

				#recovering by removing excessive missing residues
				self._recover_missing_missing_residues( mapping_ids, seqres_aln, atomres_aln, missing_sequence, insert=False )

				missing_residue_ids = seqinfo.get_missing_residue_ids( self.model.model_id)[:]
				while insertion_points :
					insertion_point = insertion_points.pop(0)
					residue_id = missing_residue_ids.pop(0)
					mapping_ids.insert( insertion_point, residue_id )


				#raise SEQRESMappingError( "SEQRES mapping failed. Missing sequence insertion is longer than SEQRES insertions." )

			###########################
			#final checking
			###########################
			if self.is_same_sequences( self.get_seqres_sequence(), self.get_sequence( residue_ids=mapping_ids) ) :
				self.seqres_atom_mapping_flag = 2
				if verbose :
					print("Global alignment is done successfuly!", file=sys.stderr) 

				return mapping_ids
			else :
				if verbose :
					print("Error: Final checking in mapping failed!", file=sys.stderr)
					print("seqres_sequence:", self.get_seqres_sequence(), file=sys.stderr)
					print("mapped_sequence:", self.get_sequence( residue_ids=mapping_ids ), file=sys.stderr)
					for i, (a, b, rid) in enumerate( zip( self.get_seqres_sequence(), self.get_sequence(residue_ids=mapping_ids), mapping_ids) ) :
						print(i, rid, a, b, file=sys.stderr)
				self.seqres_atom_mapping_flag = -2

		else :
			if verbose :
				print('Error: alignment procedure failed to generate good alignment', file=sys.stderr)
				print('Error1', len(seqres_aln), seqres_aln, file=sys.stderr)
				print('Error2', len(atomres_aln), atomres_aln, file=sys.stderr)
			raise SEQRESMappingError( "Alignment Prodecure Failed!" )
			
		if verbose :
			print("Error: This message should not be seen.", file=sys.stderr)
			print("Error1 seqres_aln ", seqres_aln, file=sys.stderr)
			print("Error2 atomres_aln", atomres_aln, file=sys.stderr)
			print("Error3 missing_sequence", missing_sequence, file=sys.stderr)
			print("Error4 seqres_insertions", seqres_insertions, file=sys.stderr)
		raise SEQRESMappingError( "Final Checking in global alignment mapping failed." )

			
	def is_missing_residue_id_mapped( self ) :
		'''
		returns true if missing residue ids are mapped, otherwise returns false.
		'''
		if self.seqres_atom_mapping_flag :
			return True
		else :
			return False

	def get_residue( self, residue_id ) :
		'''
		returns Residue for the given residue_id.
		
		If residue_id does not have corresponding Residue,
		this function will return None.
		'''
		if residue_id in self.residues :
			return self.residues[residue_id]

	def __getitem__( self, *args ) :
		chainid = self.chainid
		resnum = None
		icode = ' '

		if len(args) == 1 :
			resnum = args[0]
			
		elif len(args) == 2 :
			resnum = args[0]
			icode = args[1]
			if icode :
				icode = ' '

		elif len(args) == 3 :
			chainid, resnum, icode = args
			if icode :
				icode = ' '

		residue_id = build_residue_id( chainid, resnum, icode )
			
		return self.get_residue( residue_id )

	###############################
	#Good.
	#does not need modification!!
	###############################
	#one concern is that the residues
	#might need to be checked for the
	#missing residue or not.
	#since the new seqres mapping will
	#add missing residues into the residue list!
	#But the missing residue information
	#might be handled by residue.is_belong_to_chain.
	#Let me think about it.
	###############################
	def get_within_chain_residue_ids( self ) :
		'''
		returns list of residue ids. 
		This function is the standard way of getting 
		a list of residues in the protein chain!
		'''
		residue_ids = []
		for id in self.residue_ids :
			residue = self.get_residue( id )
			if residue.is_belong_to_chain() :
				residue_ids.append( id )
		return residue_ids

	def get_within_chain_residues( self ) :
		residues = []
		for id in self.residue_ids :
			residue = self.get_residue( id )
			if residue.is_belong_to_chain() :
				residues.append( residue )
		return residues

	#?Is this necessary?
	def is_residue_within_chain( self, residue_id ) :
		''' returns True for the residue_id belong to the chain.
		Otherwise, returns False. 
		Note that "belong to the chain" means the residue ID appears before
		chain termination (TER record).
		'''
		if  residue_id in self.residues :
			residue = self.get_residue(residue_id)
			return residue.is_belong_to_chain()
		else :
			return False
	'''
	#does not make sense.. 
	def get_residue_within_chain( self, residue_id ) :
		if residue_id in self.residues :
			residue.
			return self.residues[ residue_id ]
	'''

	def add_residue( self, residue_id, atom ) :
		'''
		adds a new residue into the Chain class.
		'''

		atom_line = atom.atom_line

		residue = None
		if residue_id in self.residues :
			residue = self.residues[residue_id]
		else :
			residue = Residue( id=residue_id, name=atom.get_residue_name() )
			self.residues[residue_id] = residue
			self.residue_ids.append( residue_id )

		if atom_line[:6] == 'TER   ' :
			self.terminated = True
		
		if self.terminated == False :
			#changed!
			#now residues need to ask the question
			#self.residues_within_chain.add(residue_id) 
			#registering self (chain) to the residue

			#this should be changed!
			if residue.chain :
				pass
			else :
				residue.chain = self
			residue.add_atom(atom)

		return residue
	
	def get_bfactor( self ) :
		'''
		Average B-factor of the residues belonged to the chain
		'''
		bfactors = [ residue.get_bfactor() for residue in self.residues.values() ]
		
		#count_none = bfactors.count(None)
		#for i in xrange(count_none) :
			#bfactors.remove( None )
		bfactors = [s for s in bfactors if s!=None]
		
		if bfactors :
			return sum(bfactors)/len(bfactors)
		else :
			return None


	def get_simple_ca_based_sequence( self ) :
		''' 
		returns sigle residue letter string for all residues in the chain.

		Caution!! This function does not check if the chain is ended or not!
		Trailing amino acids (e.g. ligand or enzyme substrate, product, etc will be
		considered as a valid amino acid).
		'''
		fp = io.StringIO()
		
		for rid in self.residue_ids :
			residue = self.get_residue( rid )
			resname = residue.get_name()
	
			if residue.has_ca() :
				if resname in amino_acids :
					fp.write( amino_acids[resname] )
				elif resname == 'MSE' :
					fp.write( 'M' )
				else :
					fp.write( 'X' )
		return fp.getvalue()



	def get_chain_continuous( self ) :
		'''
		returns list of residues with coordinates
		And check if the residue is connect to the next residue
		'''
		residues = self.get_within_chain_residues()
		
		continuous = []
		for r1, r2 in zip(residues[:-1], residues[1:]) :
			if r1.is_continuous( r2 ) :
				if verbose :
					print(r1.id, r1.name, "continuous", r1.ca_distance( r2 ), file=sys.stderr)
				continuous.append( True )
			else :
				if verbose :
					print(r1.id, r1.name, "discontinuous", r1.ca_distance( r2 ), file=sys.stderr)
				continuous.append( False )

		return continuous
			

	def residue_id2residue_index( self, residue_id ) :
		'''
		returns residue index for the residue id.

		If residue ID is not found or erroneously formatted,
		then this function returns None.
		'''

		if residue_id in self.residue_ids :
			return self.residue_ids.index( residue_id )
		
