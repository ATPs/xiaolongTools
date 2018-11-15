'''
This module contains classes that deal with information saved in a PDB file.
'''
import sys
from . import verbose

#verbose = 1

class MoleculeInfo :
	'''
	contains information about the polymers defined in the PDB file.
	This class can provide information about molecules 
	by the MOL_ID or CHAIN ID in PDB file.
	'''
	def __init__( self ) :
		self.names = {}
		self.mol_id2chains = {} #one molecule can be defined duplicated in multiple chains
		self.chain2mol_id = {}  #one chain can only be mapped to one moledule.. #maybe not???

	def print_info( self, fp=sys.stderr ) :
		'''
		debugging function
		'''
		print(self.names, file=fp)
		print(self.mol_id2chains, file=fp)
		print(self.chain2mol_id, file=fp)


	def get_molecule_id_and_chain_groups( self )  :
		'''
		returns molecule_id list of lists of chain IDs and chain IDs.
		The chains from same molecules are grouped into same list.
		'''
		return list(zip( *list(self.mol_id2chains.items()) )) #unzipping.. :)

	def get_chain_ids( self ) :
		'''
		returns list of all chain IDs.
		'''
		return list(self.chain2mol_id.keys())

	def get_molecule_name( self, chainid ) :
		'''
		returns name of the molecule of the chain ID.
		'''
		if chainid in self.chain2mol_id :
			return self.names[ self.chain2mol_id[chainid] ]
		else :
			print("WARNING: No chain id found @MolculeInfo.get_chain_name", chainid, file=sys.stderr)
			return None

	def get_molecule_id( self, chainid ) :
		if chainid in self.chain2mol_id :
			return self.chain2mol_id[chainid]
		else :
			print("WARNING: No chain id found @MolculeInfo.get_chain_name", chainid, file=sys.stderr)
			return None


	def get_name( self, mol_id ) :
		'''
		returns name of the molecule speicified by the molecule ID.
		'''
		if mol_id in self.mol_id2chains  :
			return self.names[ mol_id ]
		else :
			print("WARNING: No mol_id found @MolculeInfo.get_name", mol_id, file=sys.stderr)
			return ""

	def add_info( self, mol_id, name, chainids ) :
		'''
		Adds information. 
		
		Internal use only.
		'''
		assert mol_id not in self.names, 'Duplicated mol_id %s'%mol_id
	
		self.names[ mol_id ] = name
		self.mol_id2chains[mol_id] = chainids
		if verbose :
			print('@MoleculeInfo.add_info', chainids)
		for chainid in chainids :
			self.chain2mol_id[chainid] = mol_id

	def nucleotide_checker_from_name_string( self, name ) :
		'''
		checks if the given name describing a protein or not.
		Returns true if the name is a protein, false for everything else.
		
		Probably not the best way to find the identity of the molecule.
		Possibly the SEQRES records (residues) can solve this problem.
		'''
		return (name.startswith( 'DNA' ) or name.startswith('RNA') ) and name.find("5'")!=-1 and name.find("3'")!=-1 and name.find("*")!=-1 

	def parse_COMPND( self, compnd_list ) :
		'''
		Main parsing function for the COMPND record in PDB.

		Internal use only.
		'''
		######################
		#internal note:
		#This function is hastely written.
		#It might need to rewritten to make it more robust and 
		#maintainable code.
		######################

		mol_id_marker = 'MOL_ID:'
		name_marker = 'MOLECULE:'
		chain_marker = 'CHAIN:'
		mol_id = None
		name = ""
		chainids = []
		compound_lines = compnd_list[:]
		while compound_lines :
			line = compound_lines.pop(0) #get the first line
			line = line[10:].strip() #removing junks
			if line.startswith( mol_id_marker ) :
				if mol_id != None :
					assert name, 'problem in finding name\n%s'% compnd_list
					assert chainids, 'problem in finding chain ids\n%s'% compnd_list
					self.add_info( mol_id, name, chainids )
				
				chainids = []
				name = ""
				mol_id = line[len(mol_id_marker):].strip() ##MOL_ID: 1;
				if mol_id[-1] == ';' :
					mol_id = mol_id[:-1]
				else :
					print("WARNING: COMPND, semi-colons are expected at the end of the line!", line, file=sys.stderr)
				continue
			elif line.startswith( name_marker ) :
				name = line[len(name_marker):].strip()
				#if the name part is emptry or does not end with ";",
				#usually mean name is so long..
				if not name or name[-1] != ';' :
					line = compound_lines.pop(0)[10:].rstrip()
					while line[-1] != ';' :
						name += line
						line = compound_lines.pop(0)[10:].rstrip()
					else :
						name += line
				if name[-1] == ';': 
					name = name[:-1]
				else :
					print("WARNING: COMPND name should end with ';'", name, file=sys.stderr)
				
			elif line.startswith( chain_marker ) :
				for s in line[len(chain_marker):].split() :
					if len(s) == 2 and ( s[-1] == ',' or s[-1] == ';' ) :
						chainids.append( s[:-1] )
					elif len(s) == 1 :
						chainids.append( s )
					else :
						print("WARNING: COMPND chains should end with , or ;", s, file=sys.stderr)
						print("WARNING1 Trying to recover...", end=' ', file=sys.stderr) 
						recovery = [x for x in s if x.isalnum()]
						
						if len(recovery) == 1 :
							chainids.append( recovery )
							print("WARNING2 Recovery succeded...", recovery, file=sys.stderr)
						
							
						continue
				else :
					if s[-1] == ';' :
						pass #correctly parsed!!
					elif len(s) == 1 :
						pass #probably the chain id line is the last line of the COMPND..
					else :
						line = compound_lines.pop(0)[10:].rstrip()
						while line[-1] != ';' :
							for s in line.split() :
								if len(s) == 2 and ( s[-1]==',' or s[-1]==';' ) :
									chainids.append( s[:-1] )
								elif len(s) == 1 :
									chainids.append( s )
								else :
									print("WARNING: COMPND chains should end with , or ;",s, file=sys.stderr)
									print("WARNING1 Trying to recover...", end=' ', file=sys.stderr) 
									recovery = [x for x in s if x.isalnum()]
									
									if len(recovery) == 1 :
										chainids.append( recovery )
										print("WARNING2 Recovery succeded...", recovery, file=sys.stderr)
								
									continue
							else :
								if len(s) == 1 :
									break
							line = compound_lines.pop(0)[10:].rstrip()
						else :
							for s in line.split() :
								if len(s) == 2 and ( s[-1]==',' or s[-1]==';' ) :
									chainids.append( s[:-1] )
								elif len(s) == 1 :
									chainids.append( s )
								else :
									print("WARNING: COMPND chains should end with , or ;",s, file=sys.stderr)
									print("WARNING1 Trying to recover...", end=' ', file=sys.stderr) 
									recovery = [x for x in s if x.isalnum()]
									
									if len(recovery) == 1 :
										chainids.append( recovery )
										print("WARNING2 Recovery succeded...", recovery, file=sys.stderr)
					
		else :
			#assert name, 'problem in finding name \n%s'% compnd_list
			assert chainids, 'problem in finding chain ids \n%s'% compnd_list
			self.add_info( mol_id, name, chainids )

		if verbose :
			self.print_info()


