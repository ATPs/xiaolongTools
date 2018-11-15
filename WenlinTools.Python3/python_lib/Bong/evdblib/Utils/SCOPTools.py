
import io
import weakref

class SCOPNode :
	'''
	Represents each element in the SCOP classification.
	'''
	def __init__( self, suid=None, sccs=None, description=None, comment=None, ordering_number=None ) :
		self.suid = suid
		self.sccs = sccs
		self.description = description
		self.comment = comment
		self.ordering_number = ordering_number 
		
		#ordering number is a new feature that will help
		#to ordering nodes within each scop level
		#in structural similarities.

		self.children = None
		self.parent = None

		#additional info in SCOP hierarchy.
		#these parts are most historical but
		#more intuitive for human uage

		#'cl', 'cf', etc. two letter identifier for the node type in scop classification file..
		self.type = None 

		#domain only information
		self.scopid = None #scopid

	def __iter__( self ) :
		if self.children :
			return self.children.__iter__()

	def add( self, child ) :
		if not self.children :
			self.children = []

		self.children.append( child )
		child.set_parent( self )

	def set_parent( self, parent ) :
		self.parent = weakref.ref( parent )

	def get_parent( self ) :
		if self.parent :
			return self.parent()

	def is_domain( self ) :
		if self.type == 'px' :
			return True
		else :
			return False

	def is_root( self ) :
		if self.parent :
			return False
		else :
			return True

	def get_domains( self ) :
		'''
		retuns list of domains in the hierarchy of this node.
		'''
		if self.is_domain() :
			return [ self ]
		else :
			domains = []
			for child in self :
				domains.extend( child.get_domains() )

			return domains

	def get_parents( self, include_root=True ) :
		'''
		returns list of parent domains up to the root
		'''
		if self.is_root() :
			if include_root :
				return [ self ]
			else :
				return []

		parents = []
		parent = self.get_parent()
		while parent :
			parents.insert(0, parent)
			parent = parent.get_parent()

		return parents

	def is_matching_sccs( self, sccs ) :
		'''
		returns True if given sccs partially matches the node's sccs value.
		Otherwise returns False.

		For exameple sccs, "a" will be True fo all "All alpha proteins".
		'''
		if sccs.endswith( '.' ) or sccs.count('.') == 3 :
			return self.sccs.startswith(sccs)

		#need to add dot to avoid confusion!
		#and sccs does not have dots!!
		elif sccs.count('.') < 3 :
			return self.sccs.startswith( sccs+'.' )

		else :
			raise ValueError( sccs, "SCOP sccs value cannot have more than three dots!" )

			

class SCOPTree :
	'''
	Parses and stores pretty much all info from SCOP Parsable files.
	'''
	def __init__( self, hie, des, com, cla=None ) :
		#key->suid
		#value-> SCOPNode
		self.nodes = {} 

		if hie :
			self.build_tree( hie )
		if des :
			self.import_description( des )
		if com :
			self.import_comment( com )

		#currnetly
		#cla file is not used!

		#cla and hie contain essentially
		#same information in different format.
		#and cla is actually more strict format.

	def get_root_node( self ) :
		if self.nodes :
			root = self.nodes['0']
			return root

	def print_domain_hierarchy( self, selected_domain_ids=None, format='scopid cl cf sf fa dm sp px', field_seperator='\t', record_seperator='\n' ) :
		'''
		returns string of domain hierarchy.

		This function is useful in printing out generic string 
		for SCOP description in hierarchy.
		'''
		root = self.get_root_node()

		fp = io.StringIO()

		domains = root.get_domains()
		format_list = format.split()

		for domain in domains :
	
			if selected_domain_ids :
				if not domain.scopid in selected_domain_ids :
					continue

			#Getting all nodes between the domain and the root
			hierarchy = domain.get_parents()
			hierarchy.append( domain )
			descriptions = []

			for key in format_list :
				if key == 'scopid' :
					descriptions.append( domain.scopid )
					continue

				if key == 'suid' :
					descriptions.append( domain.suid )
					continue

				if key == 'range' :
					descriptions.append( domain.description.split()[1] )
					continue

				if key == 'pdbid' :
					descriptions.append( domain.description.split()[0] )
					continue

				if key == 'sccs' :
					descriptions.append( domain.sccs )
					continue

				if key.endswith( '.com' ) :
					stripped_key = key[:-4] 
					for node in hierarchy :
						if node.type == stripped_key :
							if node.comment :
								descriptions.append( node.comment )
							else :
								descriptions.append( '' )
							break
					else :
						raise ValueError( "Keyword not found!", key )

					continue

				for node in hierarchy :
					if node.type == key :
						descriptions.append( node.description )
						break
				else :
					raise ValueError( "Keyword not found!", key )

			fp.write( field_seperator.join( descriptions ) )
			fp.write( record_seperator )

		return fp.getvalue()


	def import_comment( self, com ) :
		'''
		parses comment file and imports 
		comments into each node.
		'''

		fp = open( com )
		for l in fp :
			if l.startswith( '#' ) :
				continue

			suid, comment = l.strip().split( ' ', 1 )
			if suid in self.nodes :
				node = self.nodes[suid]
				if not node.comment :
					node.comment = comment
				else :
					node.comment += comment

	def import_description( self, des ) :
		'''
		parses des file
		and add the information into each node.
		'''

		fp = open( des )
		for l in fp :
			if l.startswith( '#' ) :
				continue

			suid, nodetype, sccs, scopid, description = l.strip().split('\t')
			if suid in self.nodes :
				node = self.nodes[suid]
				node.description = description
				if scopid != '-' :
					node.scopid = scopid
				node.type = nodetype
				node.sccs = sccs
			else :
				print("The node to import description not found!", suid) 

	def build_tree( self, hie ) :
		'''
		parses hie, SCOP dir.hie.scop.txt file,
		or hierarchy file.
		And this function builds the entire tree structure
		based on the links defined in hie file.
		'''
		fp = open( hie ) 

		for l in fp :

			#comment
			if l.startswith( '#' ) :
				continue

			nodeid, parentid, childrenid_list = l.split()


			if nodeid in self.nodes :
				node = self.nodes[nodeid]
			else :
				node = SCOPNode( suid=nodeid )
				self.nodes[nodeid] = node


			if parentid == '-' :
				parent = None
			else :
				parent = self.nodes[parentid]


			if node.get_parent() == parent :
				pass
			else :
				print("Error! Parent information is missing or inconsistent!")

			
			if childrenid_list == '-' :
				children = None
			else :
				for childid in childrenid_list.split(',') :
					if childid in self.nodes :
						child = self.nodes[childid]
					else :
						child = SCOPNode( suid=childid )
						self.nodes[ childid ] = child

					node.add( child )

if __name__ == '__main__' :
	import sys
	hie = sys.argv[1]
	des = sys.argv[2]
	com = sys.argv[3]
	rep = sys.argv[4]
	scoptree = SCOPTree( hie, des, com )

	representatives = open(rep).read().split()
	
	#sys.stdout.write( scoptree.print_domain_hierarchy(selected_domain_ids = set(representatives) ) )

	#scoptree.print_domain_hierarchy( self, selected_domain_ids=None, format='scopid cl cf sf fa dm sp px', field_seperator='\t', record_seperator='\n' ) :
	sys.stdout.write( scoptree.print_domain_hierarchy(selected_domain_ids = set(representatives), format='suid sccs scopid pdbid range cl cl.com cf cf.com sf sf.com fa fa.com dm dm.com sp sp.com px px.com' ) )
