'''
This module contains Range classes
that parse multiple formats of range definition strings.
and facilitate the extraction of residues from
PDB or FASTA.
'''
import re
from . import ParseError
from .PDB import build_residue_id

verbose = 0
debug = 0
	
class Range :
	'''
	Range class is designed to provide conveninet interfaces
	to parse domain definitions of different fomats.

	Domain boundaries or the "range" can be divided into 
	one or more contiguous fragments in sequence.
	
	Each contiguous fragment is modelled by ContiguousRange objects.
	
	Note that Range is a abstract class provides ideas about interfaces
	provided to the specific type of dervied classes of Range
	like SequenceRange or StructureRange.
	'''

	def __init__( self )  :
		self.contiguous_ranges = []
		self.parsing_function_list = []

	def __getitem__( self, i ) :
		return self.contiguous_ranges[i]

	def __iter__( self ) :
		return self.contiguous_ranges.__iter__()

	def __str__( self ) :
		return ','.join( [ str(fragment) for fragment in self ] )

	def __len__( self ) :
		return len( self.contiguous_ranges )
	 
	def get_breaking_positions( self ) :
		'''
		This function will return a list of breaking points,
		i.e. from range "1-100,101-220,221-300,310-400"
		the returned list will be
		[101,221].

		Note that neither 300 nor 310 will not be marked as breaking
		points becase there are insertions between 300 to 310.
		'''
		breaking_points = []
		for r, l in zip( self[:-1], self[1:] ) :
			if verbose :
				print('checking breaking points')
				print(r, l)
				print(r.end, l.start)

			if r.end +1 == l.start :
				breaking_points.append( l.start )
				if verbose :
					print(l.start, "added!")

		return breaking_points
	
	def parse( self, string_of_certain_range_format, stride=1 ) :
		'''
		trying to parse the given range input format
		using the list of parsing functions defined in parsing_function_list.

		Note that the parsing function list has priority in parsing,
		if the earlier one is successful in parsing, the later functions will 
		not be even tried.
		'''
		for parse_fn in self.parsing_function_list :
			try :
				parse_fn( string_of_certain_range_format, stride=stride )
			except RangeParseError :
				#if error happens check the next function
				continue
			else :
				#if no error happens
				#accept the result.
				break
		else :
			#visiting this block means that the
			#parsing function list is empty or 
			#all functions in the funtion list 
			#failed to parse the given input range format.

			raise RangeParseError( 'All tried parsing methods failed; %s'
				%str(self.parsing_function_list) )

class ContiguousRange :
	'''
	Pretty much like slice object defined in standard library.

	Abstract class.
	'''

	def __str__( self ) :
		raise NotYetImplemented
	
	def parse( self, start_stop_string, stride=1 ) :
		raise NotYetImplemented

	def get_start( self ) :
		raise NotYetImplemented

	def get_end( self ) :
		raise NotYetImplemented
		

class SequenceContiguousRange( ContiguousRange ) :
	'''
	ContiguousRange for Sequences.

	Start and end class variable can be used for the start and end positions
	for the contiguous region.
	'''
	def __init__(self, start=None, end=None ) :
		self.start = start
		self.end = end

	def __str__( self ) :
		return '%s-%s'%(self.start, self.end)

	def parse( self, start_stop_string, stride=1) :
		if not start_stop_string :
			return

		i, j = start_stop_string.split( '-' )
		self.start, self.end = int(i), int(j)

	def get_start( self ) :
		return self.start
	def get_end( self) :
		return self.end
 
	def get_overlap( contiguous_range ) :
		'''
		get overlap bewteen SequenceContiguousRanges.
		'''
		i = self.get_start()
		j = self.get_end()

		n = contiguous_range.get_start()
		m = contiguous_range.get_end()

		x = max(i,n) #max of starts
		y = min(j,m) #min of ends

		if x > y :
			return None #no overlap found!

		new_contig = SequenceContiguousRange()
		new_contig.start = x
		new_contig.end = y

		return new_contig


class SequenceRange( Range ) :
	'''
	Range class that deal with Sequence.
	'''
	def __init__( self ) :
		Range.__init__(self)

		############################################
		#the function should be added into the list
		#for parsing!!
		#further parsing functions should be added very simply.
		self.parsing_function_list.append( self.parse_scoplike_sequence_format )

	def add_contiguous_range( self, contig ) :
		self.contiguous_ranges.append( contig )

		
		

	def is_sorted( self ) :
		'''
		returns True if the fragement is in sorted order in the sequence.
		Otherwise, returns False.

		Note that this function does not check if the ranges are overlapping or not.
		'''
		for contig1, contig2 in zip( self[:-1], self[1:] ) :
			if contig1.get_start() < contig2.get_start() and contig1.get_end() < contig2.get_end() :
				pass
			else :
				return False
		else :
			return True

	def get_inversed_range( self, whole_length=None ) :
		'''
		get inversed SequenceRange object.
		New SequenceRange object will be made.

		Due to the design issue, whole_length should be given
		to get the inversion after the last fragment.
		If whole_length is not given, the last inverted fragment might be omitted!

		Note that the inverted ranges will be ordered from low to high indices.
		'''
		sorted_contigs = [ (contig.get_start(), contig) for contig in self ]
		sorted_contigs.sort()
		sorted_range = SequenceRange()
		for i, contig in sorted_contigs :
			sorted_range.add_contiguous_range( contig )

		if debug :
			for i, contig in sorted_contigs :
				print(contig)

		new_range = SequenceRange()
		if not len(self) :
			return None

		#inverted fragment before first contig
		new_start = 1
		new_end = sorted_range[0].get_start() - 1
		if new_start == None or new_end == None :
			raise SequenceRangeError( "Inversion failed before first contig." )

		if new_end >= new_start :
			new_contig = SequenceContiguousRange( start=new_start, end=new_end )
			new_range.add_contiguous_range( new_contig )

		for contig1, contig2 in zip( sorted_range[:-1], sorted_range[1:] ) :
			new_start = contig1.get_end()
			new_end = contig2.get_start()

			if new_start == None or new_end == None :
				raise SequenceRangeError( "Inversion failed!" )
			
			new_start = new_start + 1
			new_end = new_end - 1
			
			if new_start - new_end == 1 :
				#just breaking point
				#cannot make any contig...
				continue
			elif new_end < new_start :
				raise SequenceRangeError( "Fragment cannot be generated!" )
			
			new_contig = SequenceContiguousRange( start=new_start, end=new_end )
			new_range.add_contiguous_range( new_contig )
		else :
			if whole_length :
				last_contig = self[-1]
				new_start = last_contig.get_end()+1
				new_end = whole_length

				if new_start <= new_end :
					new_contig = SequenceContiguousRange( start=new_start, end=new_end )
					new_range.add_contiguous_range( new_contig )
			

		if len(new_range) :
			return new_range
		else :
			return None
	

	def parse_scoplike_sequence_format(self, range_string, stride=1) :
		'''
		parses SCOP like sequence range format parser.

		Example:
		"1-100,112-200"

		The above example means the sequence boundary is defined from
		position 1 to 100 and then start from 112 to 200.
		Each contiguous fragment is divided by ','.
		'''
		contiguous_divider = ','
		for contiguous_string in range_string.split(contiguous_divider) :
			contig = SequenceContiguousRange()
			try :
				start, end = contiguous_string.split('-')
				contig.start = int(start)
				contig.end = int(end)
				contig.stride = stride

				#print contig, start, end

				self.add_contiguous_range( contig )
			except :
				raise RangeParseError( "Error in parsing simple scoplike format" )

class PDBContiguousRange( ContiguousRange ) :
	'''
	ContiguousRange for Sequences.

	Start and end class variable can be used for the start and end positions
	for the contiguous region.
	'''

	def __init__(self) :
		self.chain_id = None

		self.start_resnum = None
		self.start_icode = ''

		self.end_resnum = None
		self.end_icode = ''

	def __str__( self ) :

		s = ''
		if self.chain_id != ' ' and self.chain_id :
			s = '%s:' % self.chain_id
		elif self.chain_id == ' ' :
			s = ''
		elif not self.chain_id :
			s = ''
	

		if s :
			if (self.start_resnum==None and self.end_resnum==None ):
				return s
			else :
				return s + '%s%s-%s%s'%(self.start_resnum,self.start_icode,self.end_resnum, self.end_icode)
		else :
			return '%s%s-%s%s'%(self.start, self.end)

	def get_start( self ) :
		'''
		returns start residue id.
		
		if start_residue_id was not defined,
		returns None
		'''

		if self.chain_id !=None and self.start_resnum !=None :
			if self.start_icode :
				return build_residue_id( self.chain_id, self.start_resnum, self.start_icode )
			else :
				start_icode = ' '
				return build_residue_id( self.chain_id, self.start_resnum, start_icode )

	def get_end( self ) :
		'''
		returns end residue id.
		
		if end was not defined,
		returns None
		'''

		if self.chain_id !=None and self.end_resnum !=None :
			if self.end_icode :
				return build_residue_id( self.chain_id, self.end_resnum, self.end_icode )
			else :
				end_icode = ' '
				return build_residue_id( self.chain_id, self.end_resnum, end_icode )




class PDBRange( Range ) :
	'''
	Range class that deal with structure Range.
	'''
	def __init__( self ) :
		Range.__init__(self)

		#PDBRange specific settings
		#full range starting with chain id
		self.scoprange_re1 = r'^(.):(-?[0-9]+)([a-zA-Z]?)-(-?[0-9]+)([a-zA-Z]?)$'
		#full range without chain id
		self.scoprange_re2 = r'^(-?[0-9]+)([a-zA-Z]?)-(-?[0-9]+)([a-zA-Z]?)$'
		#full chain with chainid
		self.scoprange_re3 = r'^(.):$'
		self.scoprange_re4 = r'^-$'

		############################################
		#the function should be added into the list
		#for parsing!!
		#further parsing functions should be added very simply.
		self.parsing_function_list.append( self.parse_scoprange_format )

	def get_unique_chain_ids( self ) :
		uniq_ids = []
		for contig in self :
			if not contig.chain_id in uniq_ids :
				uniq_ids.append( contig.chain_id )
		return uniq_ids

	def add_contiguous_range( self, contig ) :
		self.contiguous_ranges.append( contig )
	
	def parse_scoprange_format( self, range_string, stride=1 ) :
		'''
		string format parser for the scop range
		'''
		#first split into fragments using comma, ,
		fragments = range_string.strip().split(',')

		for s in fragments :
			contig = PDBContiguousRange()
			#matching template #1
			match = re.match( self.scoprange_re1, s )
			if match :
				(contig.chain_id, contig.start_resnum, contig.start_icode,
					contig.end_resnum, contig.end_icode ) = match.groups()
				self.add_contiguous_range( contig )
				continue

			#matching template #2
			match = re.match( self.scoprange_re2, s )
			if match :
				#no chain id need to be set!
				(contig.start_resnum, contig.start_icode,
					contig.end_resnum, contig.end_icode) = match.groups()
				self.add_contiguous_range( contig )
				continue

			match = re.match( self.scoprange_re3, s ) 
			if match :
				contig.chain_id = match.group(1)
				#no start and end residue information is needed!
				self.add_contiguous_range( contig )
				continue
			
			match = re.match( self.scoprange_re4, s )
			if match :
				if len(fragments) == 1 :
					#nothing to be done!
					pass 
				else :
					raise RangeParseError( "SCOP range re4 has problem!" )

				self.add_contiguous_range( contig )
				continue
			
			raise RangeParseError( "SCOP range parsing has problem!", range_string )


class RangeParseError( ParseError ) :
	'''
	This error means that all the available parsers in the 
	Range class failed to correctly parse the given range format..
	'''
	pass

class SequenceRangeError( ParseError ) :
	'''
	This error means that the SequenceRange is not right.
	'''
	pass
