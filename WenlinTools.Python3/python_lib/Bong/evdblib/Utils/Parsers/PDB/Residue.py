from evdblib.Utils.Parsers.PDB import *

class Residue :
	def __init__( self, id=None, chain=None, name='', number=0 ) :
		self.atoms = []
		self.id = id
		self.name = name
		self.number = number
		#amino acid
		self.three2one_common = amino_acids
		self.chain = chain

		self.standard_aa = ''
		if self.name in amino_acids :
			self.standard_aa = self.name
		self.biological_residue = True

	def get_chainid( self ) :
		return self.chain.chainid

	def get_resnum( self ) :
		return self.id[1]

	def __len__( self ) :
		'''
		len( residue ) -> number of atoms.
		'''
		return len(self.atoms)

	def print_summary( self ) :
		'''
		use for debugging.
		'''
		print(self.id, self.name, self.standard_aa, 'missing:',self.is_missing_residue(), 'biological:', self.is_biological_residue(), 'chain:', self.is_belong_to_chain())

	def is_missing_residue( self ) :
		'''
		returns True if the residue is missing residue
		that does not have any atoms with modelled coordinates.

		Otherwise, it returns False.
		'''
		for atom in self.atoms :
			if atom.has_coordinate() :
				return False
		else :
			return True


	def is_biological_residue( self ) :
		'''
		returns Trues if the residue belongs to the biological sequences
		defined in DBREF.
		Otherwise, returns False.
		'''
		if self.biological_residue :
			return True
		else :
			return False

	def has_ca( self ) :
		'''
		returns True if the residue has CA atom, otherwise it returns False.
		'''
		return len( [x for x in self.atoms if x.is_ca_atom()] )


	def has_full_backbone_atoms( self ) :
		'''
		return True if the residue all backbone atoms,
		N, CA, C, O.
		'''
		nitrogens = len( [x for x in self.atoms if x.is_nitrogen()] )
		ca_atoms = len( [x for x in self.atoms if x.is_ca_atom()] )
		carbons = len( [x for x in self.atoms if x.is_carbon()] )
		oxygens = len( [x for x in self.atoms if x.is_oxygen()] )

		if nitrogens and ca_atoms and carbons and oxygens :
			return True
		else :
			return False

	def get_ca_atoms( self ) :
		'''
		returns list of CA atoms.
		It returns list becasue there might be multiple alternative locations.
		'''
		ca_atoms = [atom for atom in self.atoms if atom.is_ca_atom()]
		if len( ca_atoms ) >= 1 :
			return ca_atoms
		else :
			return []


	def get_ca_atom( self ) :
		'''
		returns a CA atom.
		In case, there are multiple CA atoms exist, 
		the highest occupancy CA atom will be returns.

		TypeError will be raised when no CA atom is found.
		'''

		atoms = self.get_ca_atoms()

		if len(atoms) == 1 :
			return atoms[0]

		elif atoms :
			#returning the atom of biggest occupancy
			oa_list = [ (atom.get_occupancy(), atom) for atom in atoms ]
			oa_list.sort()
			return oa_list[-1][-1] 
		else :
			raise TypeError( "No CA atom found in the residue %s"% str(self.id) )
		
	def get_ca_record( self, use_standard_aa=True, alanine_fall_back=False, number=None, resnum=None, chainid=None, **kwargs ) :
		'''
		returns a CA ATOM record.

		If the residue is modified and designated as HETATM,
		you can use "use_standard_aa" flag to change HETATM into ATOM record.

		This standardization is sometimes useful, 
		because some structural aligners including
		FAST and TMalign cannot use HETATM record as alignable positions.

		Also the standard designation cannot be done due to various reasons,
		the HETATM can be replaced to ATOM with alanine. 
		However, this might be problematic, so the default is not using this
		fallback option.

		Importantly atom serial number (keyword: number), residue number (resnum), and Chain ID (chainid)
		can be specified. This is very useful feature for building or normalizing PDB file.

		Optionally, other keyward arguments for Atom.get_atom_record function can be specified.

		This function raises an TypeError exception, when use_standard_aa is on 
		but it could not find the standard amino acid.
		'''

		#first get the CA atom
		atom = self.get_ca_atom()

		#For the case of GOOD AMINO ACID
		if self.name in amino_acids :
			return atom.get_atom_record( number=number, resnum=resnum, chainid=chainid, **kwargs )

		elif use_standard_aa :
			resname = self.get_standard_amino_acid_name()

			if resname in amino_acids :
				return atom.get_atom_record( type='ATOM  ', number=number, resnum=resnum, chainid=chainid, resname=resname, **kwargs )
			
			elif alanine_fall_back :
				return atom.get_atom_record( type='ATOM  ', number=number, resnum=resnum, chainid=chainid, resname='ALA', **kwargs)

			raise TypeError( "No Standard amino acid was found." )
		else :
			return atom.get_atom_record( number=number, resnum=resnum, chainid=chainid, **kwargs )

	def get_atom_records( self, use_standard_aa=True, alanine_fall_back=False, number=None, resnum=None, chainid=None, **kwargs ) :
		'''
		returns ATOM records of Atoms belong to this resiude.

		If the residue is modified and designated as HETATM,
		you can use "use_standard_aa" flag to change HETATM into ATOM record.

		This standardization is sometimes useful, 
		because some structural aligners including
		FAST and TMalign cannot use HETATM record as alignable positions.

		Also the standard designation cannot be done due to various reasons,
		the HETATM can be replaced to ATOM with alanine. 
		However, this might be problematic, so the default is not using this
		fallback option.

		Importantly atom serial number (keyword: number), residue number (resnum), and Chain ID (chainid)
		can be specified. This is very useful feature for building or normalizing PDB file.

		Optionally, other keyward arguments for Atom.get_atom_record function can be specified.

		This function raises an TypeError exception, when use_standard_aa is on 
		but it could not find the standard amino acid.
		'''

		#first get all of atoms
		atoms = self.atoms
		atom_fp = cStringIO.StringIO()

		#For the case of GOOD AMINO ACID
		if self.name in amino_acids :
			for atom in atoms :
				atom_rec = atom.get_atom_record( number=number, resnum=resnum, chainid=chainid, **kwargs )
				if atom_rec :
					atom_fp.write( atom_rec )
					if number != None :
						number = number+1
			return atom_fp.getvalue()

		elif use_standard_aa :
			resname = self.get_standard_amino_acid_name()

			if resname in amino_acids :
				for atom in atoms :
					atom_rec = atom.get_atom_record( type='ATOM  ', number=number, resnum=resnum, chainid=chainid, resname=resname, **kwargs )
					if atom_rec :
						atom_fp.write( atom_rec )

						if number != None :
							number = number+1

				return atom_fp.getvalue()
			
			elif alanine_fall_back :
				for atom in atoms :
					atom_rec = atom.get_atom_record( type='ATOM  ', number=number, resnum=resnum, chainid=chainid, resname='ALA', **kwargs)
					if atom_rec :
						atom_fp.write( atom_rec )
						if number != None :
							number = number+1
				return atom_fp.getvalue()

			raise TypeError( "No Standard amino acid was found." )
		else :
			for atom in atoms :
				atom_rec = atom.get_atom_record( number=number, resnum=resnum, chainid=chainid, **kwargs ) 
				if atom_rec :
					atom_fp.write( atom_rec )
					if number != None :
						number = number+1
			return atom_fp.getvalue()
			
		
	def get_closest_distance( self, residue ) :
		'''
		returns smallest distance between two residues.
		'''
		distances = []
		for atom1 in self.atoms :
			if atom1.has_coordinate() :
				for atom2 in residue.atoms :
					if atom2.has_coordinate() :
						d = atom1.distance( atom2 )
						distances.append( d )
		return min( distances )
			
	def is_continuous( self, residue, min_ca=3.65, max_ca=3.95 ) :
		'''check the current CA atom and the given residue
		and returns True if the two resdiues has CA and the distance between them is around 3.8A.
		'''
		if self.has_ca() and residue.has_ca() :
			ca1 = self.get_ca_atom() 
			ca2 = residue.get_ca_atom()
		
			d = ca1.distance( ca2 ) 
			if min_ca < d < max_ca :
				return True
			else :
				False
		else :
			return False
		

	def ca_distance( self, residue ) :
		'''check the current CA atom and the given residue
		and returns True if the two resdiues has CA and the distance between them is around 3.8A.
		'''
		if self.has_ca() and residue.has_ca() :
			ca1 = self.get_ca_atom() 
			ca2 = residue.get_ca_atom()
		
			d = ca1.distance( ca2 ) 
			return d

	def is_belong_to_chain( self ) :
		''' 
		returns True if the residue is belong to the chain.
		Otherwise returns False.
		
		The checking of belonging to the chain is simply done by
		checking the appearance of the residue before the "TER" record in
		the PDB file.
		
		Note:
		This simple checking sometimes causes problems,
		when "HOH"s (waters) appear at the begining of the chain.
		I added an additional checking routine to fix those cases!
		'''
		if self.chain == None :
			return False 

		if self.chain != None and self.name == 'HOH' :
			return False 
		#potentially need to check for the HETATM 
		#since they are likely not to belong to the chain
		return True

		
	def get_name( self ) :
		'''
		returns residue name appeared in the ATOM record.
		'''
		return self.name

	def get_standard_amino_acid_name( self ) :
		'''
		returns standard amino acid residue name.
		'''
		if self.standard_aa :
			return self.standard_aa

		name = self.name
		if name in amino_acids :
			return name
		else :
			if self.chain :
				seqinfo = self.chain.seqinfo
				if seqinfo :
					stdname = seqinfo.get_standard_resname( self.id )
					if stdname :
						name = stdname
						self.standard_aa = stdname
		return name

					
	def get_standard_amino_acid_name_letter( self ) :
		'''
		returns standard amino acid residue name in one letter
		'''
		resname = self.get_standard_amino_acid_name()
		if resname in self.three2one_common :
			return self.three2one_common[resname]
		elif resname :
			return 'X'
		else :
			return ''

	def get_name_letter( self ) :
		'''
		returns one letter aa residue name appeared in the ATOM record.
		
		Note that the amino acid name is not the standard amino acid, 
		this function will return X.
		'''
		if self.name in self.three2one_common :
			return self.three2one_common[self.name]
		else :
			return 'X'

	def get_bfactor( self ) :
		'''
		Average B-factor of the atoms belonged to the residue
		'''
		bfactors = [ atom.get_bfactor() for atom in self.atoms ]
		
		
		#removes None from the bfactor list
		#count_none = bfactors.count(None)
		#for i in xrange(count_none) :
			#bfactors.remove( None )
		bfactors = [s for s in bfactors if s!=None]
		
		if bfactors :
			return sum(bfactors)/len(bfactors)
		else :
			None

	def is_part_of_protein( self ) :
		'''
		checking for the belonging of the protein.
		
		WARNING: this function is not safe to get the
		correct amino acid residue.

		Instead, is_belong_to_chain() should be used.
		'''
		for atom in self.atoms :
			if atom.is_ca_atom() :
				return True
		return False

	##The following two member functions use external definition of residue_id
	##currently residue_id means tuple ( chainid, int(resnum), insertion_code )
	## the build_residue_id function on top. :)
	def add_atom( self, atom ) :
		self.atoms.append( atom )
		if not self.name :
			self.name = atom.get_residue_name()
		if not self.id :
			self.id = atom.get_residue_id()
			self.number = self.id[1]
	def get_id_string( self ) :
		if self.id :
			return str(self.id[0]) + '\t' + str(self.id[1]) + '\t' + str( self.id[2] )
			
	def get_number( self ) :
			return self.number
	def get_insertion_code( self ) :
			return self.id[2]
	def write_record( self, fp ) :
		for atom in self.atoms :
			atom.write_record( fp )
