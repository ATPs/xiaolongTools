import sys

from evdblib.Utils.Parsers.PDB import verbose
from .MoleculeInfo import MoleculeInfo
from .SequenceInfo import SequenceInfo

class PDBInfo :
	'''
	PDBInfo class manages information parsed from header section
	in PDB file. 

	Currently this class is dealing "list of molecules in PDB",
	"experimental methods", "resolution", "seqres", "ligand type and status".
	'''
	def __init__( self, pdb = None ) :
		'''
		Initialize the PDBInfo class instance.
		
		If file or fp is given,
		the parse() function will automatically be called.
		'''
		self.pdb = pdb
		if pdb :
			self.fp = pdb.get_fp()
			
		self.experiment = ''
		self.resolution = -1
		self.pdb_id = ''
		
		self.HEADER = []
		self.COMPND = [] #COMPND record strings one line per item. :)
		self.EXPDTA = [] #EXPDTA 
		self.RESOLN = [] #REMARK   2 RESOLUTION.
		self.MISSNG = [] #REMARK   
		self.SEQADV = [] #necessary to map mutations or changes introduced in the sequence
		self.SEQRES = []
		self.MODRES = []

		#best model from NMR experimental detail
		self.BSTMDL = [] #REMARK 210 BEST REPRESENTATIVE
		
		self.DBREF = [] #necessary to map biologically meaningful part of sequence 
		self.DBREF1= [] #necessary to map biologically meaningful part of sequence 

		self.moleculeinfo = MoleculeInfo()
		self.sequenceinfos = {}
		self.active_sites = [] #not using right now..
		self.dbref_regions = {} #chain_id mapped to the start residue_id and end residue_id

		if pdb :
			self.parse()

	def get_sequence_info( self, chainid ) :
		'''
		returns SequeceInfo instance for the given polymer Chain ID.
		'''
		return self.sequenceinfos.get(chainid)

	def get_chain_ids( self ) :
		'''
		returns list of chain ids.
		'''
		return self.moleculeinfo.get_chain_ids()
	

	def get_molecule_name( self, chainid ) :
		'''
		returns molecule name for the given chain id.

		If the correct header is not formatted,
		it will return None
		'''
		return self.moleculeinfo.get_molecule_name( chainid )

	def get_molecule_id( self, chainid ) :
		'''
		returns molecule id for the given chain id.

		If header does not exist, returns None.
		'''
		return self.moleculeinfo.get_molecule_id( chainid )


	def get_best_representative_model( self ) :
		'''
		returns best representative NMR model specified
		at the REMARK 210.
		
		if PDB does not contain multiple model, or 
		REMARK 210 is not parsed correctly, 
		it returns "default" which is the reserved keyword
		for the default model.
		'''
		if self.BSTMDL :
				if verbose :
					print("returning best model", self.BSTMDL)
					print("parsed number:",self.BSTMDL[0].split()[-1])

				s = self.BSTMDL[0].split()[-1]
				if s and ':' in s :
					model_num = (s.split(':')[-1]).strip()
				else :
					model_num = s

				#some NMR record cases the model number is 
				#NULL... strange!
				if model_num == 'NULL' :
					return 1
				else :
					try :
						return int(model_num)
					except ValueError :
						print("WARNING: Best Model Number is not integer.")
						return 1

		else :
			return -1
		
				
		
	def get_resolution( self ) :
		'''
		returns resolution of the structure when the information is available.
		Otherwise returns -1.
		'''
		return self.resolution
	
	def get_experiment( self ) :
		return self.experiment

	def parse( self ) :
		'''
		Parses non-coordinate information from PDB.

		This function is the main parser 
		for the non-coordinate information
		and calls many other functions.
		'''

		fp = self.pdb.get_fp()
		####################
		#reading each records
		self.HEADER = self.parse_record(fp, "HEADER")
		self.COMPND = self.parse_record(fp, "COMPND") #get initial chain information setting from COMPND record
		self.EXPDTA = self.parse_record(fp, "EXPDTA")
		self.RESOLN = self.parse_record(fp, "REMARK   2 RESOLUTION." )
		self.BSTMDL = self.parse_record(fp, "REMARK 210 BEST REPRESENTATIVE" )
		self.MISSNG = self.parse_record(fp, "REMARK 465" )
		self.SEQADV = self.parse_record(fp, "SEQADV")
		self.SEQRES = self.parse_record(fp, "SEQRES")
		self.MODRES = self.parse_record(fp, "MODRES")
		self.DBREF  = self.parse_record(fp, "DBREF " )
		if not self.DBREF :
			self.DBREF1 = self.parse_record(fp, "DBREF1" )

		
		##########################
		#start parsing records
		self.parse_HEADER( self.HEADER )
		self.moleculeinfo.parse_COMPND( self.COMPND )
		self.parse_EXPDTA( self.EXPDTA )
		self.parse_resolution( self.RESOLN)
		chain_ids = self.moleculeinfo.get_chain_ids() # all chain ids in the PDB file
		if verbose :
			print("chain_ids:", chain_ids) 

		for chain_id in chain_ids :
			if verbose :
				print(self.SEQRES)

			sequenceinfo = SequenceInfo()
			sequenceinfo.set_name( self.moleculeinfo.get_molecule_name( chain_id ) )
			sequenceinfo.set_chain_id( chain_id )
			sequenceinfo.parse_SEQRES( self.SEQRES )
			sequenceinfo.parse_MISSNG( self.MISSNG )
			sequenceinfo.parse_MODRES( self.MODRES )
			sequenceinfo.parse_SEQADV( self.SEQADV )
			if self.DBREF :
				sequenceinfo.parse_DBREF( self.DBREF )
			else :
				sequenceinfo.parse_DBREF( self.DBREF1 )

			self.sequenceinfos[ chain_id ] = sequenceinfo
		
	def parse_HEADER( self, header_list ) :
		line = header_list[0]
		self.pdb_id = line[62:66].lower()

	def parse_EXPDTA( self, expdta_list )  :
		common_types = {'X-RAY DIFFRACTION':"X-RAY", 'SOLUTION NMR':"NMR"}
		
		s = expdta_list[0][10:79].strip()
		if s in common_types :
			self.experiment = common_types[s]
			
	def parse_resolution( self, resolution_list ) :
		#parse self.RESOLN to get the resolution
		if not len(resolution_list) :
			print("WARNING! the protein resolution is not found!!", file=sys.stderr)
			return -1.0 
		elif len(resolution_list) > 1 :
			print("WARNING! the protein resolution line is more than expected!", resolution_list, file=sys.stderr)
		try :
			self.resolution = float ( resolution_list[0][23:30].strip() ) 
		except :
			try :
				self.resolution = float ( resolution_list[0].split()[-2] )
			except :
				pass

	def parse_record( self, fp, marker ) :
		'''
		currently rewind each time this function calls.
		"it should be improved to be more efficient by using orderings in PDB"
		'''
		fp.seek(0) #rewinding
		line = fp.readline()
		record = []
		while not line.startswith( marker ) :
			line = fp.readline()
			if not line :
				#print >> sys.stderr, "WARNING!", marker ,"records was not found!", self.pdb_fn
				return []
		while line.startswith( marker ) :
			record.append(line) #save compound line
			line = fp.readline() #get next compound line
			if not line :
				#print >>sys.stderr, "WARNING!", marker, "records was not parsed completely!", self.pdb_fn
				break
		
		fp.seek( -len(line), 1 ) 
		return record

