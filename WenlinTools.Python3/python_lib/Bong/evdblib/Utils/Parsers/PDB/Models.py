'''
Models module read in the PDB file
and manage multiple models in the PDB file
for NMR.
'''
import sys
from evdblib.Utils.Parsers.PDB import verbose
from evdblib.Utils import is_eof

from .Atom import Atom
from .Chain import Chain

class Models :
	'''
	Container class for Model.

	This class instance will behave pretty much
	like a dictionary except that the content is constant.

	This class also provide convenient access to the
	default Model when multiple models exist.

	Usage:
	models.get() 
	models['default']
	
	The two functions will retrieve the default model for 
	most of the cases.

	Also other models can be directly accessed 
	as in normal dictionary instances.
	'''

	def __init__( self, pdb=None, no_sequence_mapping=False ) :
		'''
		Initialize Models using given PDB class instance.

		If PDB is not given, the models will not contain
		any data.

		If PDB contains correctly formatted non-coordinate information 
		(like Headers and remarks as in PDB file downloaded from PDB database),
		the primary sequence information will be used to map missing residues
		and modified residues.
		
		If PDB does not contain any non-coordinate information 
		or no_sequence_mapping is True,
		it will skip the mapping of primary sequence mapping using SEQRES record
		and MODRES record.
		
		Skipping of the sequence mapping might help the parsing of the PDB file faster,
		but it is not recommended.
		'''
		if pdb :
			self.fp = pdb.get_fp()

		self.pdb = pdb
		self.models = {}
		self.default_model_id = -1

		#################################
		#parse cooridnate section information.
		if self.fp :
			self.parse( self.fp )
		#################################

		if verbose :
			print("List of Models")
			print(self.models)

		#################################
		#Adding information from
		#non-cooridnate section
		#PDB class instance is passed 
		#for this purpose
		#################################
		if self.pdb and self.pdb.pdbinfo :

			#Setting default model
			self.default_model_id = self.pdb.pdbinfo.get_best_representative_model()
			#checking for the error proof default model setting
			if self.models and self.default_model_id not in self.models :
				print("WARNING: Parsing of representative NMR model failed!", self.default_model_id, file=sys.stderr)
				model_ids = list(self.keys())
				model_ids.sort()
				self.default_model_id = model_ids[0]
		#################################

		if verbose :
			print("BEST model:", self.default_model_id)
	
	def get( self, model_id=-1 ) :
		'''
		returns the model of the model_id.
		
		If no model_id is given, this function will return
		default model for the NMR PDB.
		'''
		
		if model_id == 'default' or model_id==-1 :
			return self.models.get( self.default_model_id )
		elif model_id :
			return self.models.get( model_id )
		else :
			return self.models.get( self.default_model_id )

	def __getitem__( self, model_id ) :
		if model_id == None :
			return self.get( self.default_model_id )
		elif model_id == 'default' or model_id== -1  :
			return self.get( self.default_model_id )
		else :
			return self.models[model_id]

	def __contains__( self, model_id ) :
		return model_id in self.models

	def iterkeys( self ) :
		return iter(self.models.keys())

	def keys( self ) :
		return list(self.models.keys())

	def iteritems( self ) :
		return iter(self.models.items())

	def values( self ) :
		return list(self.models.values())

	def parse( self, fp ) :
		'''
		This function parses the input stream fp.
		
		1. adjust for the start point!
			read fp until the (MODEL, ATOM, HETATM, ANISOU, TER or ENDMDL)
			if the TER or ENDMDL is the first, raise Error.

		2. read model (or models)
			and save them into the dictionary, models.
		'''

		line = fp.readline()
		while line[:6] not in ( "MODEL ", "ATOM  ", "HETATM", "ANISOU", "TER   ", "ENDMDL" ) :
			line = fp.readline()
			if not line :
				raise EOFError( "Cooridnate section was not found!" )
		else :
			fp.seek( -len(line), 1 ) #put the line back!

		while not is_eof( fp ) :
			model = Model( self.pdb )
			if not model.has_atom() :
				break
			model_id = model.get_model_id()
			self.models[model_id] = model

			if verbose :
				print("MODEL %s is parsed!"% model_id)
 
class Model :
	'''
	Container for polymer chains.
	'''
	def __init__( self, pdb=None, no_sequence_mapping=False ) :
		self.pdb = pdb
		self.model_id = None
		self.chain_ids = [] #keeps appearing order
		self.chains = {}
		self.atoms = [] #probably not necessary.
		#self.non_atoms = []

		self.no_sequence_mapping = no_sequence_mapping
		if pdb :
			self.parse(  pdb.get_fp() )
			if not no_sequence_mapping :
				for chain in self.chains.values() :
					chain.map_seqres_record()
					chain.set_biological_sequence()

	def has_atom( self ) :
		'''
		returns True if the Model instance has Atoms.
		Otherwise, returns False.

		This is intended for checking if the NULL model is generated
		at the end of MODEL Parsing.
		'''
		return len(self.atoms) != 0

	def __getitem__( self, chainid ) :
		'''
		x['A'] -> y

		Chain ID 'A' will be returned.
		'''
		return self.chains[chainid]


	def get_model_id(self) :
		return self.model_id

	def get_chain( self, chain_id ) :
		if chain_id in self.chains :
			return self.chains[chain_id]
		else :
			print("WARNING: No chain", chain_id, "found in the pdb file!", file=sys.stderr)
			#raise ValueError( chain_id, "No Chain found." )
			return None

	def add_chain( self, chain_id ) :
		'''
		adds chain to the Model.
		'''
		chain = None
		if chain_id in self.chains :
			chain = self.chains[ chain_id ]
		else :
			chain = Chain( chain_id, model=self )
			self.chains[chain_id] = chain
			self.chain_ids.append( chain_id )

		return chain


	def add_atom( self, atom ) :
		'''
		Adds the ATOM line to the MODEL and
		register Chains and Residues according to the 
		parsed information from ATOM line.
		'''
		self.atoms.append( atom )
		residue_id = atom.get_residue_id()
		if residue_id :
			chain_id = atom.get_chain_id()
			chain = self.add_chain( chain_id )
			chain.add_residue( residue_id, atom )

			return True
			

	def parse_MODEL( self, model_line ) :
		'''
		parses Model ID from MODEL line.
		'''
		self.model_id = int(model_line[10:14].strip())
			
	def parse( self, fp ) :
		'''
		parses the content of a Model contained in fp.
		
		If multiple models exist (NMR case), 
		only single model is parsed.

		If the MODEL format is not right,
		ModelFormatError will be raised.
		'''
		model_record_found = 0
		endmdl_record_found = 0
		line = "1"
	 	while line :
			line = fp.readline()
			if not line :
				break

			if line.startswith( "CONECT" ) or line.startswith("MASTER") or line.startswith("END   ") :
				fp.seek( -len(line), 1 ) #put the line back!
				break

			elif line.startswith( "MODEL "  ) :
				model_record_found = 1
				self.parse_MODEL( line )
				pass

			elif line.startswith( "ENDMDL" ) :
				endmdl_record_found = 1
				break

			else :
				atom = Atom()
				if not atom.parse( atom_line=line ) :
					continue
				self.add_atom( atom )
				continue

		if model_record_found == 1 and endmdl_record_found == 0 :
			raise ModelFormatError( "MODEL was found but ENDMDL was not found!" )
		elif model_record_found == 0 and endmdl_record_found == 1 :
			raise ModelFormatError( "MODEL was not found but ENDMDL was found!" )

		#Single MODEL PDB
		if model_record_found ==0 and endmdl_record_found == 0 :
			self.model_id = -1

		

			
