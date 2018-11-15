import sys, math

###############
#from evdblib.Utils.Parsers.PDB import *
##this is same as the previous line!
#importing all package utility functions defined in "__init__.py"
###############
from .. import PDB #importing PDB package!
from .PDB import * 
###############

class Atom :
	'''
	Atom class contains common information about a line in PDB coordinate section.
	'''
	def __init__( self ) :
		self.atom_line = None
		self.type = None #ATOM, HETATM  TER ?
		self.number = 0
		self.name = ""
		self.insertion_code = None
		self.alternative_position = None
		self.coordinate = None
		self.bfactor = None
		self.occupancy = None
		self.element = ""
		self.charge = ""

		self.resname = None 
		self.resnum = None
		self.chainid = None
		
		

	def write_record( self, fp ):
		'''
		writes atom line into the file descriptor fp.
		'''
		fp.write( self.atom_line )

	def get_ter_record( self, number=None, resname=None, chainid=None, resnum=None, icode=None ) :
		
		type = 'TER   '

		#if no change should be made, why not just what I already parsed from PDB file?
		if self.is_terminal() and not any( [number,resname,chainid,resnum,icode] ) :
			return self.atom_line
			
		if number == None :
			number = self.number

		if resname == None :
			resname = self.get_residue_name()

		if chainid == None :
			chainid = self.get_chainid()

		if resnum  == None :
			resnum = self.get_residue_number()

		if icode == None :
			icode = self.get_insertion_code()


		record = '%(type)6s%(number)5d      %(resname)3s %(chainid)s%(resnum)4d%(icode)s'%locals()
		ter_record_length = 27
		full_record_length = 80
		if len(record) == ter_record_length :
			record += ' '*(full_record_length-ter_record_length) + '\n'
			return record
		else :
			raise RecordBuildingError( record )


	def get_atom_record( self, type=None, number=None, name=None, altloc=None, resname=None, chainid=None, resnum=None, icode=None, coordinate=None, bfactor=None, element=None, occupancy=None, charge=None ) :
		'''
		returns an coordinate record line (ATOM, or HETATM) 

		Note that this function can be used for rewriting PDB file
		without HETATM or chaning those HETATM records to ATOM records.
		'''
		if type == None :
			type = self.type

			if not self.has_coordinate() :
				return None

		if not self.has_coordinate()  :
			return None

		#if no change should be made, why not just what I already parsed from PDB file?
		if not any( [type,number,name,altloc,resname,chainid,resnum,icode,coordinate,bfactor,element,occupancy,charge] ) :
			return self.atom_line

		if number == None :
			number = self.number

		if name == None :
			name = self.element

		if altloc == None :
			altloc = self.get_alternative_position()

		if resname == None :
			resname = self.get_residue_name()

		if chainid == None :
			chainid = self.get_chainid()

		if resnum == None :
			resnum = self.get_residue_number()

		if icode == None :
			icode = self.get_insertion_code()

		if coordinate == None :
			coordinate = self.get_coordinates()

		if bfactor == None :
			bfactor = self.get_bfactor()

		if element == None :
			element = self.get_element()

		if occupancy == None :
			occupancy = self.get_occupancy()

		if charge == None :
			charge = self.get_charge()

			
		coord_str = '%8.3f%8.3f%8.3f'%coordinate

		record = '%(type)6s%(number)5d %(name)4s%(altloc)s%(resname)3s %(chainid)s%(resnum)4d%(icode)s   %(coord_str)s%(occupancy)6.2f%(bfactor)6.2f          %(element)-2s%(charge)2s'%locals()

		record_length = 80
		if record_length == len(record) :
			return record+'\n'

		else :
			raise RecordBuildingError( record )


	def get_element( self ) :
		if self.has_coordinate() and self.atom_line :
			return self.atom_line[76:78] 
		

	def get_charge( self ) :
		if self.has_coordinate() and self.atom_line :
			return self.atom_line[78:80]

	def get_residue_id( self ) :
		try: 
			n = int(self.atom_line[22:26])
			resnum = self.atom_line[22:26]
			chain_id = self.atom_line[21]
			insertion_code = self.atom_line[26]
			#return resnum +'\t'+ chain_id +'\t'+ insertion_code
			return build_residue_id( chain_id, resnum, insertion_code )
		except :
			return None

	def get_residue_number( self ) :
		if self.resnum == None :
			try :
				self.resnum = int(self.atom_line[22:26])
			except :
				return

		return self.resnum

	def get_insertion_code( self ) :
		if self.insertion_code == None :
			try :
				self.insertion_code = self.atom_line[26]
			except :
				return

		return self.insertion_code

	def is_terminal( self ) :
		if self.type == 'TER   ' :
			return True
		else :
			return False

	def get_chainid( self ) :
		if self.chainid == None :
			self.chainid = self.atom_line[21]

		return self.chainid
		

	def get_residue_name( self ) :
		if self.resname == None :
			self.resname = self.altome_line[22:26]
		return self.resname

	def distance( self, atom ) :
		coord1 = self.get_coordinates()
		coord2 = atom.get_coordinates()
		if coord1 and coord2 :
			return math.sqrt( (coord1[0]-coord2[0] )**2 + (coord1[1]-coord2[1])**2 + (coord1[2]-coord2[2])**2 )
		else :
			print("WARNING! Coordinates are not set!", file=sys.stdout)

	def get_coordinates( self ) :
		if not self.has_coordinate() :
			print("WARNING! Coordinates cannot be extracted!", self.type)
			return 

		if self.coordinate :
			return self.coordinate
		else :
			self.set_coordinates()
			return self.coordinate

	def set_coordinates( self, atom_line='' ) :
		if atom_line == '' :
			atom_line = self.atom_line
		if atom_line :
			try :
				self.coordinate = ( float( atom_line[30:38] ), float( atom_line[38:46] ), float( atom_line[46:54] ) ) 
			except TypeError :
				print("Error, Extracting Coordinates failed!", file=sys.stdout)
				return false 
		else :
			return False
				
	def get_type( self ) :
		return self.type

	def get_number( self ) :
		return self.number

	def parse( self, fp=None, atom_line=None ) :
		if fp :
			atom_line = fp.readline()
		elif atom_line :
			pass
		else :
			return

		self.atom_line=atom_line
		self.type = atom_line[:6]

		self.element = atom_line[12:16]
		
		try :
			self.number = int(atom_line[6:11])
		except :
			return
		return 1
		
	
	def get_chain_id( self ) :
		return self.atom_line[21]

	def get_residue_name( self ) :
		if self.resname == None :
			self.resname = self.atom_line[17:20]

		return self.resname


	def is_nitrogen( self ) :
		if self.element == ' N  ' :
			return True
		else :
			return False

	def is_ca_atom( self ) :
		if self.element == ' CA ' and ( self.type == "ATOM  " or self.type == "HETATM" ) :
			return True
		else :
			return False 

	def is_carbon( self ) :
		if self.element == ' C  ' :
			return True
		else :
			return False

	def is_oxygen( self ) :
		if self.element == ' O  ' :
			return True
		else :
			return False

	def is_alternative_location( self ) :
		if self.has_coordinate() :
			if self.get_alternative_location == ' ' :
				return False 
		return True

	def get_alternative_position( self ) :
		if self.has_coordinate() :
			pass
		else :
			return None

		if self.alternative_position == None :
			self.alternative_position = self.atom_line[16]

		return self.alternative_position

	def has_coordinate( self ) :
		if self.type == "ATOM  " or self.type == "HETATM" :
			return True
		else :
			return False

	def get_bfactor( self ) :
		if self.atom_line :
			if not self.has_coordinate() :
				return 

			if self.bfactor == None :
				try :
					self.bfactor = float( self.atom_line[60:66].strip() )
				except:
					print("WARNING: parsing B-factor failed!", file=sys.stderr)
					print(self.atom_line, file=sys.stderr)
					return None

			return self.bfactor
		else :
			return None
			
	def get_occupancy( self ) :
		if self.atom_line :
			if not self.has_coordinate() :
				return

			if self.occupancy == None :
				try :
					self.occupancy = float( self.atom_line[54:60].strip() )
				except:
					print("WARNING: occupancy failed!", file=sys.stderr)
					print(self.atom_line, file=sys.stderr)
					return None

			return self.occupancy
		else :
			return None
		

class RecordBuildingError( Exception ) :
	'''
	Exceptions after checking simple test of new record built.
	'''
	pass
