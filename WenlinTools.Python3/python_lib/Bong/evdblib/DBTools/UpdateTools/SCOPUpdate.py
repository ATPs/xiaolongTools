'''
This module contains formatting tools for 
SCOP parsible files into the add_domains input files.
'''
import sys, os, re

#dir_files = ["dir.cla.scop.txt_1.69",  'dir.cla.scop.txt_1.71',  'dir.cla.scop.txt_1.75']
#rep_id_files = ['scop40.1.69', 'scop40.1.71', 'scop40.1.75']
#versions = [1.69, 1.71, 1.75]
	

class SCOPRange :
	def __init__( self, range_string ) :
		#regular expression type 1
		#full range starting with chain id
		self.scoprange_re1 = r'^(.):(-?[0-9]+)([a-zA-Z]?)-(-?[0-9]+)([a-zA-Z]?)$'
		#full range without chain id
		self.scoprange_re2 = r'^(-?[0-9]+)([a-zA-Z]?)-(-?[0-9]+)([a-zA-Z]?)$'
		#full chain with chainid
		self.scoprange_re3 = r'^(.):$'
		self.scoprange_re4 = r'^-$'

		self.range_string = range_string
		self.starts = [] #list of ResidueID class instances
		self.ends = [] #list of ResidueID class instances
		self.pdbid = ''

		self.parse( range_string )
			
	class ResidueID :
		def __init__( self, chain, residue_number, insertion_code ) :
			self.chain = chain
			self.resnum = residue_number
			self.icode = insertion_code

		def has_chainid( self ) :
			if self.chain :
				return True
			else :
				return False

		def change_null_chainid_to_A( self ) :
			if not self.has_chainid() :
				self.chain = 'A'
				return True
			return False

		def is_whole_chain( self ) :
			if self.resnum == '' and self.icode == '' :
				return True
			else :
				return False

		def __eq__ ( self, r ) :
			if self.chain == r.chain and self.resnum == r.resnum and self.icode == r.icode :
				return True
			else :
				return False

		def __str__( self ) :
			return "ResidueID <%s,%s,%s>" % (self.chain, self.resnum, self.icode)

		def is_same_chain( self, r ) :
			if self.chain == r.chain :
				return True
			else :
				return False

		def get_chain_string ( self ) :
			''' returns SCOP range chain id string, 
			for the compatibility for the old ids, the missing chain id ('') chains will return null string '', 
			otherwise the chain string (e.g. "A:") will be returned.
			'''
			if self.has_chainid() :
				return '%s:'%self.chain
			else :
				return ''

		def get_residue_number_string( self ) :
			''' returns residue number + icode, e.g. "102A" for insertion code A with residue number 102.
			'''
			return '%s%s'%(self.resnum, self.icode)


		def get_range_string( self, r ) :
			''' returns string of ranges between the self and ResidueID class instance, r.
			'''
			if self.is_same_chain( r ) :
				pass
			else :
				print("WARNING!, self and r do not belong to the same chain!", r, file=sys.stderr)

			s = ""

			if self == r :
				return self.get_chain_string()
			else :
				return self.get_chain_string() + self.get_residue_number_string() + '-' + r.get_residue_number_string()
		
		
	def parse( self, range_string ) :
		#first split into fragments using comma, ,
		fragments = range_string.strip().split(',')

		for s in fragments :
			#matching template #1
			match = re.match( self.scoprange_re1, s )
			if match :
				chain, resnum1, icode1, resnum2, icode2 = match.groups()

				self.starts.append( self.ResidueID( chain, resnum1, icode1 ) )
				self.ends.append( self.ResidueID( chain, resnum2, icode2 ) )
				continue

			#matching template #2
			match = re.match( self.scoprange_re2, s )
			if match :
				chain = ''
				resnum1, icode1, resnum2, icode2 = match.groups()

				self.starts.append( self.ResidueID( chain, resnum1, icode1 ) )
				self.ends.append( self.ResidueID( chain, resnum2, icode2 ) )
				continue

			match = re.match( self.scoprange_re3, s ) 
			if match :
				chain = match.group(1)
				resnum1 = icode1 = resnum2 = icode2 = ''
				
				self.starts.append( self.ResidueID( chain, resnum1, icode1 ) )
				self.ends.append( self.ResidueID( chain, resnum2, icode2 ) )
				continue
			
			match = re.match( self.scoprange_re4, s )
			if match :
				chain = resnum1 = icode1 = resnum2 = icode2 = ''
				if len(fragments) == 1 :
					self.starts.append( self.ResidueID( chain, resnum1, icode1 ) )
					self.ends.append( self.ResidueID( chain, resnum2, icode2 ) )

				else :
					print("WARNING! fragment parsing error!", fragments, file=sys.stderr)

				continue
			
			print("WARNING! A fragment did not parsed correctly!", s, file=sys.stderr)

	def get_original_range_string( self ) :
		return self.range_string


	def is_chain_name_correct( self ) :
		''' returns true if there is chain name
		or if no chain name is defined it returns false.
		'''

		#assumes that starts and ends should have same chain ids.
		for residue_id in starts :
			if residue_id.has_chainid() :
				continue
			else :
				return False
		return True
			


	def get_range_string( self ) :
		''' retruns reconstructed range string.
		'''

		fragments = []
		for startid, endid in zip( self.starts, self.ends ) :
			fragments.append( startid.get_range_string( endid ) )
			

		range_string = ','.join( fragments)
	
		#correction for single chain PDB with the whole structure
		if range_string == '' :
			range_string = '-'

		return range_string


	def get_number_of_fragments( self ) :
		if not starts :
			return 1

		if len(starts) == len(ends) :
			return len(starts)
	

	def change_null_chainid_to_A( self ) :
		for startid, endid in zip(self.starts, self.ends) :
			startid.change_null_chainid_to_A()
			endid.change_null_chainid_to_A()

class SCOPnode :
        def __init__( self, id='', sccs = '', description='', parent=None, range=None, version=''  ) :
                self.id = id

		#debugging
		#print id, "is made!"

		if not self.id :
			print("Error! Null ID is given!", file=sys.stderr)
			sys.exit()

                self.type = self.determine_type( id )
                self.description = description
		self.sccs = sccs
		
		self.version=version
		self.range = range

                self.container = False
                if self.type == 'cl' or self.type=='fo' or self.type == 'sf' or self.type=='fa' :
                        self.container = True
                        self.children = []

                self.parent = parent
		self.ecoddomain = None

	def register_ecodnode( self, ecoddomain ) :
		'''connect ECODnode to the SCOPnode'''
		self.ecoddomain = ecoddomain

	def is_scop_domain_id( self, id ) :
		'''quick check if the given id is scop id or not'''
		if id[:5].isalnum() and len(id) == 7 and id.count('.') <= 1 :
			return True
		else :
			return False

	is_scop_id = is_scop_domain_id

        def determine_type( self, id ) :
                if not id :
                        return None

		if self.is_scop_id( id ) :
                        self.level = 4
                        return 'do'

		self.level = id.count('.')
		if self.level == 0 :
			return 'cl'
		elif self.level == 1 :
			return 'fo'
		elif self.level == 2 :
			return 'sf'
		elif self.level == 3 :
			return 'fa'


        def find( self, id, level=4 ) :
		
		#debugging
		#print "looking for id", id ,"in", self.id

                if self.id == id :
                        return self

                elif self.container and level > self.level and self.children :
			#debugging
			#print "looking for id", id ,"in", self.id, "checking its children"
                        for node in self.children :
                                result = node.find(id, level)
                                if result :
                                        return result
			return None
		else :
			return None

        def get_description( self ) :
                return self.description

        def get_domains( self, version='' ) :
                domains = []
                if self.container :
                        for node in self.children :
                                domains.extend( node.get_domains( version=version) )
			return domains
                else :  
			if version  :
				if self.version == version :
					return ( self, )
				else :
					return []
			else :
                        	return (self, )

        def get_nodes( self, level ) :
                ''' return all nodes in the subtree with the matching level!
                level 0: cl
                level 1: fo
                level 2: sf
                level 3: fa
		level 4: do
                '''
                if self.level == level : # family level is two!
                        return [self ]
                elif self.level < level :
                        nodes = []
                        for child in self.children :
                                nodes.extend( child.get_nodes( level ) )
			return nodes
		else :
			return []

	def get_parental_node( self, level ) :
		'''return the upper level parental node according to the level'''
		node = self
		for i in range( self.level - level ) :
			node = node.parent
		return node
		

        def add( self, node ) :
                if self.container :
			if node.parent :
				node.parent.remove( node )

                        node.parent = self
                        self.children.append( node )

	add_child=add


        def remove( self, node ) :
                if self.container and node in self.children :
                        node.parent = None
                        self.children.remove( node )

        def get_root( self ) :
                if self.parent :
                        return self.parent
                else :  
                        return self

        def deduce_parent_id( self ) :
		if self.container :
			#internal nodes
                	return '.'.join( self.id.split('.')[:-1] )
		else :
			#scopdomain nodes
                	return self.sccs #family id


	def get_sf( self ) :
		if not self.container :
			return '.'.join( self.sccs.split('.')[:3] )
		return ''


	def get_pdbid( self ) :
		if self.type == 'do' :
			return self.id[1:5]


	def get_id( self ) :
		if self.type=='do' and not self.has_chainid() :
			self.change_null_chainid_to_A() 

		return self.id

	def get_range( self ) :
		if self.type=='do' and not self.has_chainid() :
			self.change_null_chainid_to_A()
		
		return self.get_pdbid() + " " + self.range.get_range_string()

	def has_chainid( self ) :
		if self.id[-2] == '_' :
			return False
		else :
			return True

	def change_null_chainid_to_A( self ) :
		if self.type=='do' and not self.has_chainid() :
			self.range.change_null_chainid_to_A() #change range definition
			self.id = self.id[:-2] + 'a' + self.id[-1]

		return self.id

	def same_id( self, dom ) :
		if self.get_id() == dom.get_id() :
			return True
		else :
			return False
			

	def same_range( self, dom ) :
		if self.type=='do' and dom.type=='do' and self.get_range() == dom.get_range() :
			return True
		else :
			return False

	def get_version( self ) :
		return self.version

	def is_upper( self, node ) :
		''' returns true if the given node is upper level node of self node
		'''
		parent = self.parent
		while parent :
			if parent == node :
				return True
			parent = parent.parent

		return False

	def is_same_or_upper( self, node ) :
		parent = self
		while parent :
			if parent == node :
				return True
			parent = parent.parent

		return False

	'''
	def find_nearest_common_ancestor( self, node ) :
                ancestor = self
                while 1 :
                        if ancestor.find( node.get_id() ) :
                                return ancestor
                        else :
                                ancestor = self.parent
                                if ancestor == None :
                                        return None
	'''

	def find_nearest_common_ancestor( self, node ) :
		ancestor1 = self
		ancestor2 = node

		while ancestor1 and ancestor2 :
			if ancestor1 == ancestor2 :
				return ancestor1

			if ancestor1.level > ancestor2.level :
				ancestor1 = ancestor1.parent
			elif ancestor1.level < ancestor2.level :
				ancestor2 = ancestor2.parent
			else :
				ancestor1 = ancestor1.parent
				ancestor2 = ancestor2.parent

		return None

		

	def get_similarity( self, node ) :

		self_index = -900
		node_index = -200
		if node == None :
			return self_index - node_index

		#print >>sys.stderr, "getting similarity between", self.get_id(), node.get_id(), "..."
		
		ancestor = self.find_nearest_common_ancestor( node )

		#print >>sys.stderr, "ancesotr node has child nodes:", ' '.join( [ n.get_id() for n in ancestor.children ] )
		if ancestor == None :
			return -100
		
		if ancestor == node or ancestor == self :
			return ancestor.level
		else :
			for i, n in enumerate( ancestor.children ) :
				if self.is_same_or_upper(n) :
					self_index = i
				elif node.is_same_or_upper(n) :
					node_index = i
			else :
				if not ancestor.children :
					print("WARNING! ancestor's children is empty!", ancestor.get_id(), file=sys.stderr)
				if self_index < 0 or node_index < 0 :
					print("Error! distance calculation has problem!", end=' ', file=sys.stderr)
					print(self_index, node_index, file=sys.stderr)
					sys.exit()

			return ancestor.level - 0.0001*abs(self_index-node_index)


class SCOPtree :
	def __init__(self, dir_fn, rep_fn='', des_fn='') :
		self.dir_fn = dir_fn
		self.rep_fn = rep_fn

		self.roots = set()
		self.domains = []
		self.id2dom = {}
		self.range2dom = {}
		self.descriptions = {}

		#####################
		self.representatives = set()
		
		self.read_description_file( des_fn )
		self.read_representative_file( rep_fn )
		self.read_scop_dir_file( dir_fn )

	def find_closest_neighbor( self, node_list, neighbor_list_list ) :
		''' calculate neighbor distances 
		'''
		neighbor_similarity = []
		new_neighbor_list = []
		for neighbor_list in neighbor_list_list :
			if not neighbor_list :
				continue

			temp_list = []
			for neighbor in neighbor_list :
				if neighbor :
					for node in node_list :
						temp_list.append( node.get_similarity( neighbor ) )
			if temp_list :
				neighbor_similarity.append(  max(temp_list) )
			else :
				neighbor_similarity.append( -100 )
			new_neighbor_list.append( neighbor_list )

		a = list(zip( neighbor_similarity, new_neighbor_list ))
		a.sort()
		return a[-1][1], a[-1][0]
		

	def is_scop_domain_id( self, id ) :
		'''quick check if the given id is scop id or not'''
		if not id :
			return False

		if id[:5].isalnum() and len(id) == 7 and id.count('.') <= 1 :
			return True
		else :
			return False

	def read_representative_file( self, fn ) :
		if not os.path.exists( fn ) :
			return 
		fp = open( fn )
		representatives = set( [ convert_null_chainid_to_A(d) for d in fp.read().split() ] )
		fp.close()

		#convert representatives starts with 'g' astral domains
		#to 'd' as in SCOP
		for id in set(representatives) :
			if id[0] == 'g' :
				representatives.remove( id ) 
				representatives.add( 'd' + id[1:] )
			
		self.representatives = representatives

	def read_description_file( self, fn ) :
		if not os.path.exists( fn ) :
			return 
		fp = open( fn )
		descriptions = {}
		for line in fp.readlines() :
			if line and line[0]=='#' :
				continue

			l = line[:-1].split('\t')
			if not l :
				continue

			id_type = l[1]
			if id_type in ['cl','cf','sf','fa','dm','sp'] :
				id = l[2]
				des = l[4]
				descriptions[id] = des

		fp.close()
		self.descriptions = descriptions


	def get_domain_from_id( self, id ) :
		if id in self.id2dom :
			return self.id2dom[id]

	def find_node( self, id, level=4 ) :
		if self.is_scop_domain_id( id ) :
			return self.get_domain_from_id( id )
			
                for root in self.roots :
                        node = root.find( id, level )
                        if node :
                                return node

	def get_domain_from_range( self, range ) :
		if range in self.range2dom :
			return self.range2dom[range]

	def exists( self, dom_info ) :
		''' check if there is domain matches to the given dom_info.
		dom_info can be an id or a range.
		'''
		if self.get_domain_from_id( dom_info ) :
			return True

		if self.get_domain_from_range( dom_info ) :
			return True

		if self.find_node( dom_info ) :
			return True

		return False

	def get( self, dom_info ) :
		''' returns domain that matches to the given dom_info.
		Otherwise, it returns None.
		'''
		if not dom_info :
			return list(self.id2dom.values())
			
		dom = self.get_domain_from_id( dom_info )
		if dom :
			return dom
	
		dom = self.get_domain_from_range( dom_info )
		if dom :
			return dom

		dom = self.find_node( dom_info )
		if dom :
			return dom
		
		return None

	def get_description( self, id ) :
		if id in self.descriptions :
			return self.descriptions[id]
		else :
			return 'NO DESCRIPTION'

	#def get_domain( self, dom_info ) :
	def read_scop_dir_file( self, fn, version='', five_classes_only=True ) :
		'''reads in scop dir file
		'''

		if not os.path.exists( fn ) :
			return
		
		fp = open( fn )
		representative_id_set = self.representatives

		for l in fp.readlines() :
			if l.startswith( '#' ) or l.startswith( '=' ) :
				continue
			l = l.split()
			#print l

			sccs = l[3]
	 		class_type = sccs[0]

			if five_classes_only and class_type in ['a','b','c','d','e','f','g', 'h','i','j','k'] :
				pass
			else :
				continue
			
			range = SCOPRange( l[2] )
			dnode = SCOPnode( id=l[0], sccs=sccs, range=range, version=version )
			scop_id = dnode.get_id()
			scop_range = dnode.get_range()

			if representative_id_set and not scop_id in representative_id_set :
				continue

			#########################
			#Register in the SCOPtree
			#just add since the changed one is removed from the system

			if scop_id in self.id2dom :
				old_dnode = self.id2dom[scop_id]
				if old_dnode.sccs != dnode.sccs :
					print("WARNING! The node already exists!", file=sys.stderr)
					print(old_dnode.sccs, old_dnode.id, old_dnode.get_range(), file=sys.stderr)
					print(dnode.sccs, dnode.id, dnode.get_range(), file=sys.stderr)

			self.range2dom[scop_range] = dnode
			self.id2dom[scop_id] = dnode

			#putting into scop tree structure
			fanode_id = dnode.deduce_parent_id()
			#print fanode_id
			fanode = self.find_node( fanode_id, level=3 )
			if fanode :
				fanode.add_child( dnode )
				continue
			else :
				fanode = SCOPnode( id=fanode_id, description=self.get_description(fanode_id) )
				fanode.add_child( dnode )

			sfnode_id = fanode.deduce_parent_id()
			#print sfnode_id
			sfnode = self.find_node( sfnode_id, level=2 )
			if sfnode :
				sfnode.add_child( fanode )
				continue
			else :
				sfnode = SCOPnode( id=sfnode_id, description=self.get_description(sfnode_id) )
				sfnode.add_child( fanode )

			fonode_id = sfnode.deduce_parent_id()
			#print fonode_id
			fonode = self.find_node( fonode_id, level=1 )
			if fonode :
				fonode.add_child( sfnode )
				continue
			else :
				fonode = SCOPnode( id=fonode_id, description=self.get_description(fonode_id) )
				fonode.add_child( sfnode )
	
			clnode_id = fonode.deduce_parent_id()
			#print clnode_id
			clnode = self.find_node( clnode_id, level=0 )
			if clnode :
				clnode.add_child( fonode )
				continue
			else :
				clnode = SCOPnode( id=clnode_id, description=self.get_description(clnode_id) )
				clnode.add_child( fonode )
				self.roots.add( clnode )

        def get_nodes( self, level ) :
                nodes = []
                for root in self.roots :
                        nodes.extend( root.get_nodes( level ) ) #superfamily level = 2

                return nodes

	def get_domains( self, version ='' ) :
		nodes = []
		for root in self.roots :
			nodes.extend( root.get_domains(version) )
		return nodes

	def get_node_ids( self, level ) :
		return [ node.get_id() for node in self.get_nodes(level) ]
	
	def get_node_members( self, id, level=4 ) :
		node = self.find_node(id, level )
		return node.get_domains()

	def get_node_member_ids( self, id, level=4 ) :
		return [ dom.get_id() for dom in self.get_node_members(id,level) ]


	def get_superfamilies( self ) :
		return self.get_nodes( 2 )

	def get_set_of_superfamily_ids( self ) :
		return set([ node.get_id() for node in self.get_superfamilies() ])
		

	def get_superfamily_members( self, sf_id ) :
		'''returns list of domains in the given superfamily by sccs sf id (e.g. a.100.1).
		If the sccs sf id does not exists, returns None
		'''
		sfnode = self.find_node( sf_id, level=2 )
				
		if sfnode :
			return sfnode.get_domains()
		else :
			return []

	def get_superfamily_member_ids( self, sf_id) :
		doms = self.get_superfamily_members( sf_id )
		return [dom.get_id() for dom in doms]
	

def convert_null_chainid_to_A( scopid ) :
	if scopid[-2] == '_' :
		return scopid[:-2] + 'a' + scopid[-1]
	else :
		return scopid


def old_main() :
#if __name__ == '__main__' :
	dir_files = ["dir.cla.scop.txt_1.69",  'dir.cla.scop.txt_1.75']
	rep_id_files = ['scop40.1.69', 'scop40.1.75']
	versions = [1.69, 1.75]

	r2dom, id2dom = {}, {}
	errored_ids = []
	for dir_fn, rep_fn, v in zip(dir_files, rep_id_files, versions) :
		fp = open( rep_fn )
		rep_ids = set([ convert_null_chainid_to_A(scopid) for scopid in  fp.read().split()] )

		fp = open( dir_fn )
		
		read_scop_dir_file( fp, r2dom, id2dom, errored_ids, v, rep_ids )

	for errorid in errored_ids :
		scop_ids = list(id2dom.keys())
		for scopid in sorted(scop_ids) :
			if scopid[1:5] == errorid[1:5] :
				print(scopid, id2dom[scopid].get_range(), id2dom[scopid].get_version())

	#for id, dom in id2dom.iteritems() :
		#print id, dom.get_range()

def main() :
	dir_files = ["dir.cla.scop.txt_1.69",  'dir.cla.scop.txt_1.75']
	rep_id_files = ['scop40.1.69', 'scop40.1.75']
	versions = [1.69, 1.75]

	old_scop = SCOPtree( dir_files[0], '' )
	new_scop = SCOPtree( dir_files[1], '' )

	#old_sf_ids = old_scop.get_set_of_superfamily_ids()
	old_sf_ids = set(old_scop.get_node_ids(1)) #folds
	#new_sf_ids = new_scop.get_set_of_superfamily_ids()
	new_sf_ids = set(new_scop.get_node_ids(1))

	new_and_old = new_sf_ids.intersection( old_sf_ids )

	print(old_sf_ids)
	print(new_sf_ids)


	print("####################################")
	print("#Common SF in both old and new scop")
	print("####################################")
	for sf in new_and_old :
		#n = set(new_scop.get_superfamily_member_ids( sf ))
		n = set(new_scop.get_node_member_ids( sf, 1 )) #actually a fold!
		#o = set(old_scop.get_superfamily_member_ids( sf ))
		o = set(old_scop.get_node_member_ids( sf, 1 )) #actually a fold!

		n_and_o = n.intersection(o)
		if not n_and_o :
			print()
			print(sf)

			for scopid in o :
				if new_scop.exists( scopid ) :
					dom = new_scop.get(scopid)
					print(dom.get_id(), dom.get_sf())

			print()

	print("###############################")
	print("#SF only in new scop")
	print("###############################")

	new_minus_old_sf = new_sf_ids - old_sf_ids
	#old_minus_new_sf = old_sf_ids - new_sf_ids
	for sf in new_minus_old_sf :
		#n = set(new_scop.get_superfamily_member_ids( sf ))
		n = set(new_scop.get_node_member_ids( sf, 1 )) #actually a fold!

		print() 
		print(sf)

		for scopid in n :
			if old_scop.exists( scopid ) :
				dom = old_scop.get(scopid)
				print(dom.get_id(), dom.get_sf())

		print()
	

	print("###############################")
	print("#SF only in old scop")
	print("###############################")

	old_minus_new_sf = old_sf_ids - new_sf_ids
	for sf in old_minus_new_sf :
		#o = set(old_scop.get_superfamily_member_ids( sf ))
		o = set(old_scop.get_node_member_ids( sf, 1 )) #actually a fold!

		print() 
		print(sf)

		for scopid in o :
			if new_scop.exists( scopid ) :
				dom = new_scop.get(scopid)
				print(dom.get_id(), dom.get_sf())

		print()
	

if __name__ == '__main__' :
	main()
