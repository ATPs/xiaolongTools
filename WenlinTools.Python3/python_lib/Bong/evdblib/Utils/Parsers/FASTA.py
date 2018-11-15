'''
This module reads a FASTA file and
returns a FASTA class instance.

Calling parse function will conveniently generate FASTA instance containing all
FASTA data in it.

ToDo...
1. Maybe we need some convinient FASTA record writing wrapper outside the FASTA class.
especially for generating multiple FASTA records.
2. Maybe we need seperate class for more OOP way of managing multiple FASTA records.
currently it is simply implemented as a list of FASTA class instances.

'''
import sys, copy
from io import StringIO
from evdblib.Utils import is_eof, rewind

from .Range import SequenceRange

debugging = 0
verbose = 0

def parse( file='', descriptor=None ) :
	'''
	returns a FASTA record read from the file or descriptor.

	Note that this function will only read one record at a time,
	but the descriptor will maintain its position after reading a FASTA record.
	'''
	return FASTA( file=file, descriptor=descriptor )


def parse_multiple( file='', descriptor=None, checksequence=True ) :
	'''
	returns a list of FASTA records read from the file or descriptor.

	Note that this function will read all records in the file 
	or descriptor (following the current position ignoring what were already read).
	'''
	fp = descriptor
	if file :
		fp = open( file )

	fastas = []
	while not is_eof(fp) :
		fastas.append( FASTA( descriptor=fp, checksequence=checksequence ) )
	else :
		if not fastas :
			raise FASTANullRecordError( 'File does not contain any FASTA records.' )

	return fastas


def parse_a3m( file='', descriptor=None, checksequence=False ) :
	'''
	returns a list of FASTA records read from the file or descriptor.

	Note that this function will read all records in the file 
	or descriptor (following the current position ignoring what were already read).
	'''
	fp = descriptor
	if file :
		fp = open( file )

	fastas = []
	while not is_eof(fp) :
		fastas.append( A3MFASTA( descriptor=fp, checksequence=checksequence ) )
	else :
		if not fastas :
			raise FASTANullRecordError( 'File does not contain any FASTA records.' )
	return fastas

	
class FASTA :
	''' 
	reads given file or descriptor if they are passed.
	Otherwise, just make an empty FASTA instance.
	
	If there are mulitple fasta records in a single file,
	the first record is read. For descriptor, the first record
	relative to the start point is read.

	For multiple fasta reading use parse_multiple_fasta method
	in this module.

	'''

	def __init__( self, file='', descriptor=None, sequencerange=None, checksequence=True ) :
		'''
		checksequence: a flag for purging non alphabetical fasta sequence
		'''
		self.header = ''
		self.sequence = ''
		self.gi = ''
		self.file = file
		self.fp = descriptor

		self.checksequence = checksequence

		if file :
			fp = open(file) 
			self.read( fp )
			fp.close()
		elif descriptor :
			self.read( descriptor )

		self.sequencerange = sequencerange

	def __len__( self ) :
		return len(self.sequence)


	def get_range_string_for_whole_sequence( self ) :
		return '1-%d'%len(self)


	def get_header( self ) :
		'''
		Returns the defline of the FASTA record.
		
		If no FASTA record is read, 
		this function will return NULL string or ''.
		'''
		return self.header

	def get_sequence( self ) :
		'''
		Returns the sequence of the FASTA record.
		
		If no FASTA record is read, 
		this function will return NULL string or ''.
		'''
		return self.sequence

	def read_from_string( self, fasta_string ) :
		'''
		reads fasta from string 
		instead of file.

		Note that this function reads only one record.
		'''
		fp = StringIO( fasta_string )
		self.read( fp )
		fp.close()

	def read( self, fp ) :
		'''
		Reads a FASTA record from the file pointer.
		'''

		################################
		#reading FASTA defline
		annotation_line = fp.readline().strip()
		while annotation_line and annotation_line[0] != '>' :
			annotation_line = fp.readline().strip()
		else : 
			self.header = annotation_line
		################################

		if debugging :
			print("annotation_line:", annotation_line.strip(), file=sys.stderr)
			print("self.header:", self.header, file=sys.stderr)

		################################
		#reading FASTA sequence
		sequence_line = ''
		line = fp.readline()
		while line and line[0] != '>' :
			sequence_line += line.strip()
			line = fp.readline()
		else :
			if line and line[0] == '>' :
				rewind( fp, len(line) ) 
		################################
		

		#checking the sequence in FASTA
		##############
		#the following code will remove all non alphabet characters
		if self.checksequence :
			self.sequence = ''.join( a for a in sequence_line if a.isalpha() )
		#do not remove non-alphabetical characters
		else :
			self.sequence = sequence_line
		##############

		if debugging :
			print("sequence_line:", sequence_line, file=sys.stderr)
			print('sequence checking test:', ''.join( a for a in sequence_line if a.isalpha() ), file=sys.stderr)


		#########################################
		#chekcing if the fasta header contains 
		#NCBI GI number in it.
		if self.header.startswith( '>gi|' ) :
			splitted_header = self.header.split( '|' )
			if len(splitted_header) > 2 :
				self.gi = splitted_header[1]
		#########################################
	
		###############################################
		#Start to check the integrity of FASTA record!
		if not self.sequence and not self.header :
			msg = '%s does not contain any FASTA record!'
			if self.file :
				raise FASTANullRecordError( msg%self.file )
			else :
				raise FASTANullRecordError( msg%str(fp) )
				
		elif not self.sequence or not self.header :
			msg = '%%s does not have correct FASTA format. header:%s sequence:%s'%(self.sequence, self.header)
			#print '%%s does not have correct FASTA format. \nheader:%s \nsequence:%s'%(self.sequence, self.header)
			if self.file :
				raise FASTAFormatError( msg%self.file )
			else :
				raise FASTAFormatError( msg%str(fp) )
		#End of integrity check!
		###############################################

	def is_empty( self ) :
		'''
		returns True when both of the annotation property and sequence
		are NULL.
		'''
		if not self.sequence and not self.header :
			return True
		else :
			return False


	def is_protein_fasta( self ) :
		'''
		returns True when the fasta record contains protein sequence
		'''

		#secondary structure test
		sec_str_elements = set( 'HCE' )
		if len(self.sequence) >= 3 and set(self.sequence).issubset( sec_str_elements ) :
			return False

		#DNA test
		dna_elements = set( 'GATCgatc' )
		if len( self.sequence ) >= 4 and set(self.sequence).issubset( dna_elements ) :
			return False

		#Non alaphabet test
		s = [ a for a in self.sequence if a.isalpha() ]
		if s :
			return True
		else :
			return False

		
			

	def is_correct_format( self ) :
		'''
		returns True either of seqeunce or annotation is NULL
		'''
		if not self.sequence or not self.header :
			return True
		else :
			return False

	def __len__(self ) :
		return len(self.sequence)

	def __str__( self ) :
		'''
		returns String of FASTA record reconstructed by this class.
		'''

		ostream = StringIO()
		self.save( fp=ostream )
		s = ostream.getvalue()
		ostream.close()

		return s

	def save( self, file='', mode='a', fp=None, char_per_line=-1 ) :
		'''
		Saves the FASTA into the file or the file descriptor, fp.

		Note that the file can be added as write ('w', by removing all the content)
		or addition ('a', by adding the FASTA record at the end of the file).
		For saving into file, the default is addition mode.

		For saving into stream or descriptor,
		the mode is always addition.

		char_per_line -1 means all sequence appear in a line.
		Otherwise, the number of characters defined in char_per_line 
		will appear in a line!
		'''

		if not file and not fp :
			raise FASTASaveRecordError( "FASTA record cannot be saved!" )
		elif self.is_empty() :
			raise FASTANullRecordError( "FASTA record is empty" )
		elif self.is_correct_format() :
			raise FASTAFromatError( "FASTA record is not correctly read." )

		if file :
			ostream = open( file, mode )
		else :
			ostream = fp

		istream = StringIO( self.sequence )

		ostream.write( self.header )
		ostream.write( '\n' )

		seq_per_line = istream.read( char_per_line )
		while seq_per_line :
			ostream.write( seq_per_line )
			ostream.write( '\n' )

			seq_per_line = istream.read( char_per_line )

		#if the output stream is opend 
		#by this function!!
		if file :
			ostream.close()

	def extract_fasta_string( self, sequencerange=None, split_fragments=False, inverse=False, change_header=False  ) :
		'''
		returns FASTA string or List of FASTA strings,
		extracted according to SequenceRange object.

		For split_fragments = True, all sequence fragments will be concatenated 
		and returns single FASTA string.
		For split_fragments = False, each fragments will be seperated into each
		FASTA string.

		Inverse = False do not change anything in the sequencerange.
		Inverse = True, makes the fragments inverted,
		i.e. 30-60,120-150 is inverted into 1-29,61-119

		Change_header = False, do not change new fasta header.
		Change_header = True, adds range information to the fasta header.
		Default: False, to prevent the errorneous parsing of query id.

		If no sequencerange is given, sequencerange == None,
		whole FASTA will be given.
		'''
		if sequencerange ==None :
			sequencerange = self.sequencerange

		if not sequencerange :
			return str(self)

		if inverse :
			sequencerange = sequencerange.get_inversed_range( len(self) )

		sequence = self.get_sequence()
		header = self.get_header()

		if change_header :
			new_header = header + ' ' + str(sequencerange)
		else :
			new_header = header

		new_sequence = '' 
		new_splitted_sequences = []
		new_splitted_headers = []

		for contig in sequencerange :
			start = contig.get_start()
			end = contig.get_end()
			if start != None :
				start = start -1
			#if end != None :
				#end = end - 1

			if split_fragments :
				new_splitted_sequences.append( sequence[start:end] )
				if change_header :
					new_splitted_headers.append( self.get_header() + ' ' + str(contig) )
				else :
					new_splitted_headers.append( self.get_header() )
			else :
				new_sequence += sequence[start:end]

		if new_sequence :
			return new_header + '\n' + new_sequence + '\n'
		else :
			new_fastas = []
			for h, s in zip( new_splitted_headers, new_splitted_sequences ) :
				new_fastas.append( h + '\n' + s + '\n' )
			return new_fastas


	def extract_fastas( self, sequencerange=None, split_fragments=False, inverse=False, change_header=False  ) :
		'''
		returns List of FASTA objects,
		extracted according to SequenceRange object.

		For split_fragments = True, all sequence fragments will be concatenated 
		and returns single FASTA
		For split_fragments = False, each fragments will be seperated into each
		FASTA and a list of FASTA objects will be returned.

		Inverse = False do not change anything in the sequencerange.
		Inverse = True, makes the fragments inverted,
		i.e. 30-60,120-150 is inverted into 1-29,61-119

		#change header=True adds range information into the header
		However, the default is False, since many cases the added 
		range hinder to correctly parsing the query id.

		If no sequencerange is given, or the object do not have
		sequencerange variable defined,
		new fasta record copy of the self will be returned.

		Note that the range is full and inverse cannot be made,
		None will be returned.
		'''

		#check input
		if sequencerange == None :
			sequencerange = self.sequencerange

		if sequencerange == None :
			return copy.deepcopy( self )

		if inverse :
			if verbose :
				print("getting inverted range...")

			sequencerange = sequencerange.get_inversed_range(len(self) )
			if not sequencerange :
				return None

			if verbose:
				print("invertedrange:", sequencerange)


		sequence = self.get_sequence()
		header = self.get_header()

		new_splitted_sequences = []
		new_splitted_headers = []


		if split_fragments :
			new_fastas = []
		else :
			#prepation for split_fragment
			if change_header :
				new_header = self.get_header() + ' ' + str(sequencerange)
			else :
				new_header = self.get_header()

			new_sequence = ''
			new_sequencerange_strings = []


		for contig in sequencerange :
			start = contig.get_start()
			end = contig.get_end()
			if start != None :
				start = start -1
			#if end != None :
				#end = end - 1

			if split_fragments :
				#build new fasta object and sequencerange
				#then save the fasta into list
				new_seq = sequence[start:end] + '\n'
				if change_header :
					new_header = header + ' ' + str(contig) +'\n'
				else :
					new_header = header + '\n'

				new_fasta = FASTA()
				new_fasta.read_from_string( new_header+new_seq )

				new_sequencerange = SequenceRange() 
				new_sequencerange.parse( '1-%d'%len(new_seq) )
				new_fasta.sequencerange = new_sequencerange

				new_fastas.append( new_fasta )
			else :
				start_num = len(new_sequence) + 1
				new_sequence += sequence[start:end]
				end_num = len(new_sequence)
				new_sequencerange_strings.append( '%s-%s'%(start_num, end_num) )
				
				

		if split_fragments :
			return new_fastas
		else :
			new_fasta = FASTA()
			new_fasta.read_from_string( new_header +'\n'+new_sequence + '\n' )
			new_fasta.sequencerange = SequenceRange()
			new_fasta.sequencerange.parse( ','.join( new_sequencerange_strings ) )

			return new_fasta

		

	def save_header_sequence( self, header='', sequence='', file='', mode='a', fp=None, char_per_line=-1 ) :
		'''
		Saves the given fasta header and sequence into the file or the file descriptor, fp.

		Note that the file can be added as write ('w', by removing all the content)
		or addition ('a', by adding the FASTA record at the end of the file).
		For saving into file, the default is addition mode.

		For saving into stream or descriptor,
		the mode is always addition.

		char_per_line -1 means all sequence appear in a line.
		Otherwise, the number of characters defined in char_per_line 
		will appear in a line!
		'''

		if not file and not fp :
			raise FASTASaveRecordError( "FASTA record cannot be saved!" )
		
		#cleanup the sequence
		if not sequence.isalpha():
			sequence = ''.join( a for a in sequence if a.isalpha() )

		if not sequence :
			raise FASTASaveRecordError( "Sequence record is empty!" )

		if not header :
			raise FASTASaveRecordError( "Header record is empty!" )

		if header[0] != '>' :
			header = '>' + header

		if header[-1] == '\n' :
			header.strip()

		if file :
			ostream = open( file, mode )
		else :
			ostream = fp

		istream = StringIO( self.sequence )

		ostream.write( self.header )
		ostream.write( '\n' )

		seq_per_line = istream.read( char_per_line )
		while seq_per_line :
			ostream.write( seq_per_line )
			ostream.write( '\n' )

			seq_per_line = istream.read( char_per_line )

		#if the output stream is opend 
		#by this function!!
		if file :
			ostream.close()



class FASTANullRecordError( Exception ) :
        def __init__( self, value ) :
                self.value = value

        def __str__( self ) :
                return repr( self.value )

class FASTAFormatError( Exception ) :
        def __init__( self, value ) :
                self.value = value

        def __str__( self ) :
                return repr( self.value )

class FASTASaveRecordError( Exception ) :
        def __init__( self, value ) :
                self.value = value

        def __str__( self ) :
                return repr( self.value )


class A3MFASTAParseError( Exception ) :
	pass


import re
range_length_regex = re.compile( "\((.*)-(.*):(.*)\)" )
class A3MFASTA( FASTA ) :
	'''
	FASTA record to contain A3M data
	'''

	def parse_hit_range( self ) :
		'''
		parse out the hit range in defline of A3M fasta record
		made by buildali program.
		For example,
		>gi|116619592|ref|YP_821748.1|(16-336:606) ABC transporter-like protein [Candidatus Solibacter usitatus Ellin6076]  E=0.0006 s/c=0.15 id=24% cov=99%
		'''
		range_length_match = range_length_regex.search( self.header.split()[0] )
		if range_length_match :
			a, b, c = range_length_match.groups() #a,b,c : start, end, hitlength
			start, end = int(a)-1, int(b)
			hitregionlen = len([k for k in self.sequence if k.isalpha()])
			if hitregionlen == end-start :
				return start, end
			else :
				return 1, hitregionlen
		else :
			return 1, len([ i for i in self.sequence if i.isalpha()])

	def parse_hit_length( self ) :
		'''
		parse out the hit length
		'''
		range_length_match = range_length_regex.search( self.header.split()[0] )
		if range_length_match :
			a, b, c = range_length_match.groups() #a,b,c : start, end, hitlength
			return int(c)
		else :
			return len([i for i in self.sequence if i.isalpha()])


	def parse_query_range( self, queryseq ) :
		'''
		calculate query start index and query end index.
		'''
		query_length = len(queryseq )
		query_start_index = 1
		for i, a in enumerate(self.sequence) :
			if a.isalpha() :
				query_start_index = i
				break

		query_end_index = query_length
		for i, a in enumerate( reversed(self.sequence) ) :
			if a.isalpha() :
				query_end_index = query_length-i
				break

		return query_start_index, query_end_index

	
	def build_query_and_hit_sequences( self, queryseq, query_start_index, query_end_index ) :
		'''
		returns query and hit sequence in BLAST result.
		'''
		query_temp = queryseq[query_start_index:query_end_index]
		hit_temp = self.sequence[query_start_index:len(self.sequence)-len(queryseq)+query_end_index]

		query_fp = StringIO()
		hit_fp = StringIO()
		
		qi = 0
		for hi in range( len(hit_temp) ) :
			if hit_temp[hi].isupper() :
				query_fp.write( query_temp[qi] )
				hit_fp.write( hit_temp[hi] )
				qi += 1
			elif hit_temp[hi].isalpha() :
				query_fp.write( '-' )
				hit_fp.write( hit_temp[hi].upper() )
			else :
				query_fp.write( query_temp[qi] )
				hit_fp.write( hit_temp[hi].upper() )
				qi += 1

		
		if qi == len(query_temp) :
			pass
		else :
			A3MFASTAParseError( 'query and hit length does not match' )

		return query_fp.getvalue(), hit_fp.getvalue()


