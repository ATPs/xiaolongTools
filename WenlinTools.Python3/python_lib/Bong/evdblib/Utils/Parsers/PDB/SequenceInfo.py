import sys
from evdblib.Utils.Parsers.PDB import *
from evdblib.Utils import is_sorted, find_all_indices

#verbose = 1

class SequenceInfo :
	'''
	This class will manage information about the primary sequences
	of a protein molecule.
	'''
	def __init__( self ) :
		self.chain_id = "" #linked chain
		self.mol_id = "" #mol_id in COMPND records
		self.name = "" #to print as a name of a sequence
		self.SEQRES = ""
		self.MODRES = "" 
		self.type = "" #type protein, rna, dna
		self.chain = None #link to Chain class instance
		self.seqres_sequence = [] #seqres residues (3 letter as in SEQRES) in list 

		##############################
		#information related to MODRES
		##############################
		#modres residues mapping to the 3 letter code as appeared in ATOM records
		self.residue_id2modified_resname = {} 
		#modres residues mapping to the STANDARD (unmodified) 3 letter code
		self.residue_id2standard_resname = {} 

		########################################
		#information related to MISSING residues
		########################################
		#the three lists have related items at the same index
		self.missing_resnames = []
		self.missing_residues = []
		self.missing_residue_models = []

		self.biological_range = [] #range given by residue_id from DBREF
		self.full_residue_ids = [] #contain the result of missing residue id mapping


	def set_biological_sequence( self, chain ) :
		'''
		builds information about biological sequence
		into residues in the chain.
		'''
		#if biological range is not
		#parsed do nothing!
		if not self.biological_range :
			return

		#reset the biological residue information.
		for residue in chain.residues.values() :
			residue.biological_residue = False

		for biological_range in self.biological_range :
			self._set_biological_sequence_pair( chain, *biological_range )

	def _guess_residue_index( self, all_residue_ids,  within_chain_ids, residue_id ) :
		'''
		fix for wrong index boundary information in DBREF.
		This function tries to guess the missing range 
		if the residues are in increasing ordering.
		'''
		if is_sorted( within_chain_ids ) :
			for i, rid in enumerate( within_chain_ids ) :
				if rid >= residue_id :
					i = max( 0, i-1 )
					return all_residue_ids.index( within_chain_ids[i] )
			else :
				return all_residue_ids.index( within_chain_ids[-1] )
		else :
			return None
				

	def _set_biological_sequence_pair( self, chain, start_residue_id, end_residue_id  ) :
		'''
		builds information about biological sequence
		into residues in the chain.
		
		This function adds info for a range.
		'''
		#read start and end of the DBREF range.
		#and convert DBREF residue ids into 
		#indices of the residues.
		#start_residue_id, end_residue_id = self.biological_range
		#fallback option 

		within_chain_residue_ids = chain.get_within_chain_residue_ids()

		if start_residue_id in within_chain_residue_ids :
			start_index = chain.residue_ids.index(start_residue_id)
		else :
			start_index = self._guess_residue_index( chain.residue_ids, 
				within_chain_residue_ids, start_residue_id )
			
			if start_index == None :
				print("DEREF range", start_residue_id, file=sys.stderr)
				if start_residue_id < within_chain_residue_ids[0] :
					print("WARNING: DBREF start position recovery succeeded.", file=sys.stderr)
				else :
					print("WARNING: DBREF start position recovery failed in mapping biological sequence from DBREF.", file=sys.stderr)
				start_index = 0


		if end_residue_id in within_chain_residue_ids :
			end_index = chain.residue_ids.index(end_residue_id)
		else :
			end_index = self._guess_residue_index( chain.residue_ids, 
				within_chain_residue_ids, end_residue_id )

			if end_index == None :
				if end_residue_id > within_chain_residue_ids[-1] :
					print("WARNING: DBREF end position recovery succeeded.", file=sys.stderr)
				else :
					print("WARNING: DBREF end position recovery failed in mapping biological sequence from DBREF.", file=sys.stderr)
				end_index = chain.residue_ids.index( within_chain_residue_ids[-1] )


		#set biological residue information
		for residue_id in chain.residue_ids[start_index:end_index+1] :
			residue = chain.get_residue( residue_id )
			if residue.is_belong_to_chain() :
				residue.biological_residue = True

		return True
				

	
	def make_missing_residue_information( self, residue_ids, res_seq, model_num, override=0 ) :
		''' 
		builds Missing Residue infomration when the record is missing!!
		if oeveride value is True, pre-exisiting missing residue information will be
		replaced. if not overide is False, the existing missing residue information will not be
		replaced.
		residue_ids are the list of residue IDs.
		res_seq is the sequence of missing residues in one letter.
		model_num is the the number of model.
		'''
		if self.missing_resnames :
			if not override :
				print("WARNING: Missing residue information update cannot be overriden. Set override 1.", file=sys.stderr)
				return False
		if not len(res_seq) and len(residue_ids) :
			print("WARNING: Residue sequence and number of residue ids do not match!!, The record will not be built!", file=sys.stderr)
			return False
			
		self.missing_resnames = [ amino_acids_rev[aa] for aa in res_seq ]
		self.missing_residues = residue_ids[:]
		self.missing_residue_models = [ (model_num,) ] * len(res_seq )
		return True

	def change_missing_residue_id( self, old_missing_residue_id, new_missing_residue_id ) :
		'''
		This function is designed to change some wrong missing residue ID.
		I found out this some cases (e.g. 1GPU) the missing residue ids were same as the 
		already existing residue IDs in polymer chain.

		This function provide ad-hoc remedy for those cases just simply put non-colliding
		residue id.
		'''
		target_indices = find_all_indices( self.missing_residues, old_missing_residue_id )
		for i in target_indices :
			self.missing_residues[i] = new_missing_residue_id

	def exists_missing_residue( self, residue_id, model_num ) :
		indices = find_all_indices( self.missing_residues, residue_id )
		for i in indices :
			if model_num in self.missing_residue_models[i] :
				return True
		else :
			return False

	def find_missing_residue_index( self, residue_id, model_num ) :
		'''
		returns the missing residue index.
		'''
		indices = find_all_indices( self.missing_residues, residue_id )
		ret = []
		for i in indices :
			if model_num in self.missing_residue_models[i] :
				ret.append( i )
		if len(ret) == 1 :
			return ret[0]
		elif len(ret) == 0 :
			raise ValueError( "No element is found." )
		else :
			raise ValueError( "More than one residue_id is found.", str([ret,model_num]) )

	def remove_missing_residue_id( self, residue_id, model_num ) :
		'''
		Removes missing residue ID.

		This function is designed for fixing overlapping records 
		in Missing Residue remarks.

		returns False if no residue index is found.
		'''
		try :
			index = self.find_missing_residue_index( residue_id, model_num ) 
		except ValueError:
			return False

		if verbose :
			print("removing residue id:", residue_id, model_num, "@ index", index, file=sys.stderr)

		del self.missing_residues[index]
		del self.missing_resnames[index]
		del self.missing_residue_models[index]

	def insert_missing_residue_information( self, residue_ids, res_seq, model_num, insert_index ) :
		''' 
		Adds Missing Residue infomration when the record is missing!!

		residue_ids are the list of residue IDs.
		res_seq is the sequence of missing residues in one letter.
		model_num is the the number of model.
		'''
		if not len(res_seq) and len(residue_ids) :
			print("WARNING: Residue sequence and number of residue ids do not match!!, The record will not be built!", file=sys.stderr)
			return False
			
		aa_seq = list( res_seq )

		if verbose :
			print("In the insertion function:", residue_ids, aa_seq, model_num, insert_index, file=sys.stderr)

		while residue_ids :
			residue_id = residue_ids.pop()
			aa = aa_seq.pop()

			if verbose and insert_index < len(self.missing_residues):
				i=insert_index
				print("before:", self.missing_residues[i-1], self.missing_resnames[i-1], self.missing_residue_models[i-1], file=sys.stderr)
				print(residue_id, aa, "inserted at position", insert_index, file=sys.stderr) 
				print("after:", self.missing_residues[i], self.missing_resnames[i], self.missing_residue_models[i], file=sys.stderr)
				print(file=sys.stderr)

			self.missing_resnames.insert(insert_index, amino_acids_rev[aa] )
			self.missing_residues.insert(insert_index, residue_id )
			self.missing_residue_models.insert( insert_index, (model_num,) )
			
		if verbose :
			print("###after insert_missing_residue_information", file=sys.stderr)
			for i, (a,b,c) in enumerate( zip( self.missing_residues, 
				self.missing_resnames, self.missing_residue_models ) ) :
				print(i,":", a, b, c, file=sys.stderr)

	def print_seqres_information( self, fp=sys.stderr ):
		for i, (a,b,c) in enumerate( zip( self.missing_residues, 
			self.missing_resnames, self.missing_residue_models ) ) :
			print(i,":", a, b, c, file=fp)
	
			

		return True
	
	def has_missing_residues( self ) :
		'''
		returns True if missing residues are defined in the Remarks section.
		Otherwise returns False.
		'''
		if self.missing_resnames :
			return True
		else :
			return False

	def set_chain( self, chain ) :
		self.chain = chain

	def print_full_residue_ids( self, chain ) :
		'''
		prints out full residue ids by chainid, residue number, and insertion code.
		But written hastely for Qian.. :(
		Only treated as an simple toy example, not to be seriously used.
		'''
		count = 0
		for residue_id in self.full_residue_ids  :
			residue = chain.get_residue( residue_id )
			letter = "x"
			if residue :
				if residue.is_belong_to_chain() :
					pass
				else :
					continue
				
				#mapping to the modified residues
				letter = residue.get_name_letter()
				if letter == "X" :
					if residue_id in self.residue_id2standard_resname :
						resname = self.residue_id2standard_resname[ residue_id ] 
						if resname in amino_acids :
							letter = amino_acids[resname]
				#print for the residues with coordinates
				count += 1
				print("%s\t%s\t%s"%residue_id + "\t" + letter + '\t' + str(count))
			else :
				try :
					index = self.missing_residues.index( residue_id ) 
				except :
					index = -1
				if index == -1 :
					print("Error: residue_id is not found!", residue_id, file=sys.stderr)
				else :
					resname = self.missing_resnames[index]
					if resname in amino_acids :
						letter = amino_acids[resname]
					else :
						letter = "X"
				
				#print for the residues without coordinates
				print("%s\t%s\t%s"%residue_id + "\t" + letter)

	def get_missing_residue_ids( self, model_num ) :
		residue_ids = []
		for residue_id, model_nums in zip( self.missing_residues, self.missing_residue_models ) :
			if model_num in model_nums :
				residue_ids.append( residue_id ) 
			elif -1 in model_nums :
				residue_ids.append( residue_id )

		return residue_ids

	def get_missing_resnames( self, model_num ) :
		resnames = []
		for resname, model_nums in zip( self.missing_resnames, self.missing_residue_models ) :
			if model_num in model_nums :
				resnames.append( resname )
			elif -1 in model_nums :
				resnames.append( resname )
		return resnames

	def missing_residue_id2resname( self, residue_id ) :
                try :
                        index = self.missing_residues.index( residue_id )
                except :
                        print("WARNING: Missing residue not found", residue_id, file=sys.stderr)
                        return None
         
                return self.missing_resnames[index]
        
        def is_missing_residue( self, residue_id ) :
                return residue_id in self.missing_residues

	def get_standard_resname( self, residue_id ) :
		if residue_id in self.residue_id2standard_resname :
			return self.residue_id2standard_resname[ residue_id ]
		else :
			print("WARNING: Modified residue not found", residue_id, file=sys.stderr)
			return None

	def get_name( self ) :
		return self.name

	def set_chain_id( self, chain_id ) :
		self.chain_id = chain_id

	def set_name( self, chain_name ) :
		self.name = chain_name
		
	def print_seqres_sequence( self ) :
		for s in self.seqres_sequence :
			print(s)

	####################
	#Better function to check if the given Chain is protein or not.
	####################
	def is_protein( self ) :
		'''
		returns True if the sequence is polypeptide chain.
		Otherwise, returns False.
		'''
		count = 0 
		for s in self.seqres_sequence :
			if verbose :
				print(s)
			if s in amino_acids :
				count += 1
		if count : #if more than 1 amin acid in seqres sequence.. 
			return True
		else :
			return False

	def parse_SEQRES( self, seqres_list, chainid="" ) :
		'''
		parse given SEQRES record that matches with the given chain ID.
		If the chainid is not given (default), the chain_id in the class instance will be used.
		'''
		if not chainid :
			chainid = self.chain_id
		
		number_of_residues = 0
		seqres_sequence = []
		for line in seqres_list :
			
			if line[11] == chainid :
				number_of_residues = int( line[13:17].strip() )
				for i in range( 19, 70, 4 ) :
					s = line[i:i+3]
					if s.isspace() :
						continue
					seqres_sequence.append( s )
		if len(seqres_sequence) != number_of_residues :
			print("WARNING: parsing of SEQRES record failed", number_of_residues, seqres_sequence, file=sys.stderr)
			return False
		
		self.seqres_sequence = seqres_sequence
		return True

	#this function is important to get the Engineered residues.. 
	def parse_SEQADV( self, seqadv_list ) :
		pass

	def parse_DBREF( self, dbref_list ) :
		if not dbref_list :
			print("WARNING: no DBREF found!", self.name, file=sys.stderr)
			return False

		dbref_lines = []
		dbref1_lines = []
		dbref2_lines = []
		for line in dbref_list :
			if line[12] == self.chain_id :
				if line.startswith( "DBREF1" ) :
					dbref1_lines.append( line )
				elif line.startswith( "DBREF2" ) :
					dbref2_lines.append( line )
				else :
					dbref_lines.append( line )

		#check the reading of DBREF record.
		if dbref_lines or dbref1_lines or dbref2_lines  :
			pass
		else :
			print("WARNING: DBREF with the chain id is not found!", self.chain_id, file=sys.stderr)
			return False


		if len(dbref1_lines) != len(dbref2_lines) :
			raise PDBInfoRecordFormatError( "DBREF1 and DBREF2 lines should be matched", "dbref1", dbref1_lines, "dbref2", dbref2_lines )
			
	
		dbref_infos = []

		#parse DBREF lines
		#there might be vectors..
		#need to remove those non-biological definitions.
		for dbref_line in dbref_lines :
			dbref_info = self.parse_DBREF_line( dbref_line )
			
			#IF db name is PDB
			#and the reference points to itself.
			#means that the region defined is vector contaminents.
			#if  dbname=='PDB' and pdbid == dbaccession :
			if dbref_info[3] == 'PDB' and dbref_info[0] == dbref_info[4] :
				continue

			dbref_infos.append( dbref_info )

		#parse DBREF1 and DBREF2 cases.
		#they supposed to be true.
		for dbref1, dbref2 in zip( dbref1_lines, dbref2_lines ) :
			dbref_info = self.parse_DBREF12_line( dbref1, dbref2 )
			dbref_infos.append( dbref_info )
		
		if not dbref_infos :
			print("WARNING: DBREF with the chain id is self PDB!", self.chain_id, file=sys.stderr)
			return False

		elif len(dbref_infos) > 1 :
			dbaccessions = set([dbref_info[4] in dbref_infos])
			if len(dbaccessions) == 1:
				for dbref_info in dbref_infos :
					self.biological_range.append( [dbref_info[1], dbref_info[2]] )
			else :
				raise PDBInfoRecordParsingError( "DBREF is ambiguous.", dbref_infos )

		else :
			#add start and end index.
			self.biological_range.append( [dbref_infos[0][1], dbref_infos[0][2]] )

		return True


	def parse_DBREF_line( self, line ) :
		'''
		Parse Single line of DBREF that is relavent to the 
		primary sequence of the protein.
		'''
		pdbid = line[7:11]
		start_id = build_residue_id( line[12], int(line[14:18]), line[18] )
		end_id = build_residue_id( line[12], int(line[20:24]), line[24] )
			
		dbname = line[26:32].strip()
		dbaccession = line[33:41].strip()
		dbid = line[42:54].strip()
		#self.biological_range.append( start_id )
		#self.biological_range.append( end_id )
		return (pdbid, start_id, end_id, dbname, dbaccession, dbid)

	def parse_DBREF12_line( self, dbref1, dbref2 ) :
		'''
		Parse DBREF1 and DBREF2 lines
		that is relevant to the protein.
		'''
		pdbid1 = dbref1[7:11]
		pdbid2 = dbref2[7:11]
		chain1 = dbref1[12]
		chain2 = dbref2[12]
		
		start_id = build_residue_id( dbref1[12], int(line[14:18]), line[18] )
		end_id = build_residue_id( dbref1[12], int(line[20:24]), line[24] )
		dbname = dbref1[26:32].strip()
		dbid = dbref1[47:67].strip()
		dbaccession = dbref2[18:40].strip()

		if pdbid1 == pdbid2 and chain1 == chain2 :
			pass
		else :
			raise PDBParsingError( "No DBREF1 and DBREF2 does not match!" )

		return (pdbid1, start_id, end_id, dbname, dbaccession, dbid)
		

	#MODRES contains non-standard amino acid mapping to the standard amino acids
	def parse_MODRES( self, modres_list ) :
		'''
		Parses MODRES record and saves the modified residues
		into the following variables;

		self.residue_id2modified_resname
		self.residue_id2standard_resname

		The keys for the two dictionary is ResidueID. 
		They looks like this: ('A', 120, ' ') 
			'A' is the Chain ID, 120 is the residue number, ' ' is the insertion code.

			This residue ID is built by the build_residue_id() defined 
			in the PDB package.
		'''
		if not modres_list :
			#since MODRES is optional it is OK without it. :)
			return True
		modified = {}
		standard = {}
		for line in modres_list :
			modified_resname = line[12:15]  #residue name in the ATOM record
			#the following line's definition came from the following example
			#^MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE                                 $
			chainid = line[16]
			if chainid == self.chain_id :
				pass
			else :
				continue
			resnum = line[18:22]
			insertion_code = line[22]
			standard_resname = line[24:27]
			#checking if the resnum is correct :)
			try:
				test = int( resnum.strip() )
			except :
				print("WARNING: Resnum is not a valid number!", resnum, file=sys.stderr)
			#building residue id for mapping
			#it is important to use this residue ID in the mapping
			#to uniquely identify the target residue
			residue_id = build_residue_id( chainid, resnum, insertion_code )
			if residue_id in modified :
				print("WARNING: Duplicated record in modres list!", residue_id, modified, file=sys.stderr)
			#in cases where modified residue names are
			#the same as standard residue names 
			#just skip the record line!!
			#save the residue id to resnames mapping
			modified[residue_id] = modified_resname
			standard[residue_id] = standard_resname
		#saving the final results
		self.residue_id2modified_resname = modified
		self.residue_id2standard_resname = standard
		return True

	def parse_MISSNG( self, missng_list ) :
		'''
		parses MISSING residues defined in REMARK 465.
		
		This function add the missing residue information
		into three lists in the class.

			self.missing_resnames
			self.missing_residues
			self.missing_residue_models
		.
		'''
		if not missng_list :
			return True
		start_mark_non_nmr = "REMARK 465   M RES C "#SSSEQI" #use shinked version to find the start..
		start_mark_nmr =     "REMARK 465     RES C "#SSSEQI" #1r5r has SSEQI instead of standard SSSEQI. :(
		parse_nmr_template = 0
		list_of_models = []
		start_line = 0
		for i, line in enumerate( missng_list ) :
			if line.startswith( start_mark_non_nmr ) :
				start_line = i+1
				break
			elif line.startswith( start_mark_nmr ) :
				start_line = i+1
				parse_nmr_template = 1
				break
		else :
			#print >>sys.stderr, "WARNING: No start point in parsing MISSING RESIDUES!", missng_list
			raise MissingResidueFormatError('No start point in parsing Missing Residues', str(missng_list))

		if parse_nmr_template and start_line >= 2 :
			assert start_line >= 2, 'start_line should be more than 2'
			line = missng_list[start_line-2] #changed for debugging??
			#normal NMR model line
			if line.startswith( "REMARK 465   MODELS" ): 
				l = line.split()
				temp_range = l[-1].split('-')
				if len( temp_range ) == 2:
					model_start, model_end = temp_range
				else :
					print("WARNING: The model range is not correct", self.name, line, file=sys.stderr)
					#raise MissingResidueFormatError( 'The model range is not correct.' )
					model_start, model_end = 0, 300
				
				try :
					model_start = int(model_start)
					model_end = int(model_end)
				except :
					model_start, model_end = 0, 300
						
				list_of_models = list(range(model_start, model_end+1))
			else :
				#assume that the PDB file contains only 1 model 
				#even if the PDB is about NMR structure
				#case 2k0j...
				list_of_models = (-1,) #treat like non-nmr..

		#list_of_models contains NMR model list..
		#for all of the missing residues.
			
		for line in missng_list[start_line:] :
			final_models = () 
			#parse nmr missing residues
			#list of models that have missing residues are predefined in list_of_models
			if parse_nmr_template == 1 and list_of_models :
				final_models = list_of_models
			elif parse_nmr_template == 1 and not list_of_models :
				print("WARNING: NMR tempalte should have list of models pre-defined!", file=sys.stderr)
				final_models = list(range(0,300))
			#parse non-nmr (X-ray) missing residues
			#models are defined in the line
			elif parse_nmr_template == 0 :
				model = line[12:14].strip()
				
				if model == '' :
					model = -1
				final_models = ( int(model), )
				
			else :
				print("WARNING: Parsing problem while MISSING residues! This messsage should not be seen!", file=sys.stderr)
				print(''.join(missng_list))
				
			chainid = line[19]
			resnum = line[22:26]
			resname = line[15:18]
			insertion_code = line[26]
			if chainid != self.chain_id :
				continue
			residue_id = build_residue_id( chainid, resnum, insertion_code )
			self.missing_resnames.append( resname )
			self.missing_residues.append( residue_id )
			self.missing_residue_models.append( final_models )

