verbose = 1

from evdblib.Utils import rewind

class PairwiseAlignmentParser :
	def __init__( self, fp ) :
		'''
		Parase an alignment from result file.

		This class is abstract class for the 
		specific alignment parsers.
		'''
		self.fp = fp
		self.parse_start_point = fp.tell()
		#self.parsed_record = self.parse()
		#self.parse_end_point = fp.tell()

	def parse( self ) :
		'''
		Abstract function for parsing
		'''
		raise NotYetImplementedError() 

	def parse_all( self ) :
		'''
		Abstract fucntion for parsing all available records
		in one file.
		'''

class HHsearchParser( PairwiseAlignmentParser ) :
	'''
	HHsearch alignment parsing tool.
	'''
	def __init__( self, record_start='No ', record_end='Done!', queryid=None, **kwargs ) :
		PairwiseAlignmentParser.__init__( self, **kwargs )

		self.record_start = record_start
		self.record_end = record_end
		self.queryid = queryid

	def _parse_query_id( self, rewind_position=0 ) :
		'''
		read the beginning of the HHsearch result file pointer
		to get the query ID in case where the query id is not given and
		cannot be safely inferred from the alignment line.
		'''
		self.fp.seek(rewind_position)
		header = self.readline()
		query_id = header.split()[1]
		self.fp.seek( self.parse_start_point )

		return queryid


	def _parse_hit_id( self, hit_header ) :
		'''
		parse the hit header line to get the 
		'''
		return hit_header.split()[1][1:]
		

	def _read_a_record( self, fp ) :
		'''
		read an alignment and retun
		list of lines like "readlines()"
		function.

		Note that this function will raise an exception
		NoMoreRecordError when it reaches end of file
		or record end marker is reached.
		'''

		#reads from a record start marker to another record start marker
		#or record end maker. 
		record_start = self.record_start
		record_end = self.record_end

		line = fp.readline()
		while line and ( line.startswith( record_start ) or line.startswith( record_end ) ) :
			line = fp.readline()
			if verbose :
				print("checking line:", line[:-1])
		else :
			if not line :
				raise InvalidFormatError( "End of file is reached erroneously before the start mark.") 
			if line.startswith( record_end ) :
				raise NoMoreRecord( "End of HHsearch result is reached." )

		content = [line]
		line = fp.readline()
		while line :
			if line.startswith( record_start ) or line.startswith( record_end ) :
				rewind( fp, len(line) )
				break
		
			content.append( line )
		else :
			#null line means end of file is reached
			#before end of record marker is reached.
			if not line :
				raise InvalidFormatError( "End of file is reached erroneously before end of record.") 
				
		return content
			
	def _parse_score_line( self, line ) :
		key_val_pairs = line.split()
		scores = {}

		for pair in key_val_pairs :
			key, value = pair.split('=')
			if key == 'Probab' :
				scores[key] = float( value ) / 100.0
			elif key == 'identities' :
				scores[key] = float( value[-1] ) / 100.0
			else :
				scores[key] = float( value )

	def _parse_alignment( self, contents ) :
		query_lines = [ l for l in contents if l.startswith('Q') ]
		hit_lines = [ l for l in contents if l.startswith('T') ]

		#HHsearch lines are repeating by 4
		line_units = 4
		alignment_line_index = 2

		#################################
		#checking alignment lines
		#################################

		if len(query_lines) != len(hit_lines) :
			if verbose :
				print("Error in _parse_alignment:")
				print("Qeury:")
				print("\n".join( query_lines ))
				print("Target:")
				print("\n".join( hit_lines ))
			
			raise InvalidFormatError( "Parsing alignment lines has problems." )

		if len( query_lines ) and not len(query_lines)%4 :
			pass
		else :
			if verbose :
				print("Error in query lines _parse_alignment.")
				print("Query:")
				print("\n".join( query_lines ))

			raise InvalidFormatError( "Parsing query alignment lines has problems." )
	
		if len( hit_lines ) and not len(hit_lines)%4 :
			pass
		else :
			if verbose :
				print("Error in hit lines _parse_alignment.")
				print("Hit:")
				print("\n".join( hit_lines ))

			raise InvalidFormatError( "Parsing hit alignment lines has problems." )

		#############################
		#parsing query line
		#############################
		
		query_aln = ''
		for i, line in enumerate( query_lines ) :
			ri = i%4
			if ri == alignment_line_index and i == alignment_line_index :
				l = line.split()
				query_start_index = int(l[-4]) - 1

				if not line[-3].isalpha() : 
					raise InvalidFormatError( "Invalid alignment line: %s"%line)
				query_aln = l[-3]
			elif ri == alignment_line_inex :

				if not line[-3].isalpha() : 
					raise InvalidFormatError( "Invalid alignment line: %s"%line)
				query_aln += line[-3]
		else :
			if not query_aln :
				raise InvalidFormatError( "Query alignment was not found!" )
				
		#############################
		#parsing hit line
		#############################
		
		hit_aln = ''
		for i, line in enumerate( hit_lines ) :
			ri = i%4
			if ri == alignment_line_index and i == alignment_line_index :
				l = line.split()
				hit_start_index = int(l[-4]) - 1

				if not line[-3].isalpha() : 
					raise InvalidFormatError( "Invalid alignment line: %s"%line)
				hit_aln = l[-3]
			elif ri == alignment_line_inex :

				if not line[-3].isalpha() : 
					raise InvalidFormatError( "Invalid alignment line: %s"%line)
				hit_aln += line[-3]
		else :
			if not query_aln :
				raise InvalidFormatError( "Hit alignment was not found!" )

		return query_start_index, hit_start_index, query_aln, hit_aln
	
		
	
	def parse( self ) :
		'''
		Parse a record out of the file pointer.
		'''
		#read whole record first!
		content = self._read_a_record( self.fp )

		#parse a line  "No <hit_number>\n"
		hit_number = int( content.pop(0).split()[-1] )
		#parse a line ">hit_name\n"
		hit_header = content.pop(0)
		if hit_header[0] != '>' :
			raise InvalidFormatError( "HHsearch hit header line is wrong: %s" %hit_header )
		hitid = self._parse_hit_header( hit_header )

		score_line = content.pop(0)
		try :
			scores = self._parse_score_line( score_line )
		except IndexError as ValueError :
			raise InvalidFormatError( score_line )
			
		blank_line = content.poo(0)
		if blank_line != '\n' :
			raise InvalidFormatError( "Blank Line after header is not correct!" )

		query_start_index, hit_start_index, query_aln, hit_aln = _parse_alignment( content )
		
		raw_score = scores['Score']
		normalized_score = score['Probab']

		aln_method_record = PairwiseAlignmentMethodRecord( start_postion1=query_start_index, start_position2=hit_start_index, alignment1=query_aln, alignment2=hit_aln, raw_score=raw_score, norm_score=normalized_score, queryid=self.queryid, hitid=hitid, additional_scores=scores )

		return aln_method_record


	def parse_all( self ) :
		'''
		Parse all alignments in the given file.
		'''
		self.parse_start_point = fp.tell()
		if not self.queryid :
			self.qeuryid = self._parse_query_id()

		alignments = PairwiseAlignmentRecords()
		
		while 1:
			try :
				aln_method_record = self.parse()
			except NoMoreRecordError :
				break
			alignments.add( aln_method_record )

		return alignments
		
		
		
		
		

##################################
#Record Parsing Exceptions
##################################
class NoMoreRecordError( Exception ) :
	'''
	Indicator class for end of file.
	'''
	pass

class InvalidFormatError( Exception ) :
	pass

