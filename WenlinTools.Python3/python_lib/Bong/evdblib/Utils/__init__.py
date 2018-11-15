import os, io, glob, sys

#Valid pathname characters is defined for safe conversion of a string
#into a pathname.
valid_pathname_characters = '.-+ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

verbose=0

def string2pathname( input_value ) :
	'''
	converts input_value into valid pathname 
	by replacing all non valid pathname characters
	into underscores.
	
	valid pathname characters: ".+-[A-Z][a-z][0-9]"

	Note that the pathname starting with '-' is prohited
	to prevent parameter confusion by leading '-'.
	'''

	str_fp = io.StringIO()

	#removing the leading '-' if it is there.
	if input_value and input_value[0] == '-' :
			input_value = '_' + input_value[1:] 

	#convert all characters do not belong to valid characters
	#to underscores.
	for c in input_value :
		if c in valid_pathname_characters :
			str_fp.write( c )
		else :
			str_fp.write( '_' )

	return str_fp.getvalue()
	

def parse_sequence_filename( filename, known_extentions=('.fa', '.fasta', '.seq', '.sequence', '.msa', '.aln', '.hhm', '.pnp', '.cnp', '.cnp.len','.pdb', '.ca', '.dat', '.ent','.psi', '.bla', '.blast', '.psiblast', '.a3m', '.a2m', '.bb') ) :
	'''
	returns tuple of strings;
	( <directory>, <base_filename_without_extension>, <known extension> )

	This is useful function for analyzing the structure of filename.
	When known extension is not found, all base
	'''
	directory, basename = os.path.split( filename )
	for ext in known_extentions :
		if basename.endswith( ext ) :
			return directory, basename[:-len(ext)], ext

	return directory, basename, ''


def build_sequence_filename( directory, basename, ext ) :
	return os.path.join( directory, basename + ext )



def parse_profile_filename( filename, known_extentions=( '.fa', '.fasta', '.seq', '.msa', '.aln', '.hhm', '.pnp', '.cnp', '.cnp.len','.psi', '.bla', '.blast', '.psiblast', '.a3m', '.a2m' ) ) :
	'''
	returns <directory>, <base_filename>, <iteration>, <known_extention>
	'''
	directory, basename = os.path.split( filename )

	for ext in known_extentions :
		if basename.endswith( ext ) :
			basename = basename[:-len(ext)]
			try :
				basenames = basename.split('.')
				iteration = str( int(basenames[-1]) )
				basename = '.'.join( basenames[:-1] )
			except :
				iteration = ''
	
			return directory, basename, iteration, ext
	else :
		return directory, basename, '', ''


def build_profile_filename( directory, basename, iteration, ext ) :
	'''
	returns a profile filename.

	Note that the iteration can be "*" for wild card matching. :)
	'''
	return os.path.join( directory, '%s.%s%s'%(basename, str(iteration), ext) )


def build_profile_filename( directory, basename, iteration='', ext='') :
	'''
	retrurns iteration added filename
	'''
	filename = os.path.join( directory, basename )
	if iteration :
		return filename + '.' + str(iteration) + ext
	else :
		return filename + ext

def find_max_iteration( directory, basename, ext ) :
	'''
	returns maximum number of iterations found for the given profile
	filename signature.
	'''
	profiles = glob.glob( build_profile_filename( directory, basename, '*', ext ) )
	iterations = [int(parse_profile_filename(fn)[2]) for fn in profiles ]

	if not iterations :
		raise ValueError( "No matching files found!", build_profile_filename( directory, basename, '*', ext ) )

	return max( iterations )
		

	

def peek( file_descriptor, bytes_to_peek=1 ) :
	'''
	returns a character read from the filedescriptor but
	the filedescriptor position did not change.
	
	If the EOF is met, the peek function will return Null string.
	'''

	read_bytes = file_descriptor.read(bytes_to_peek)
	file_descriptor.seek(-len(read_bytes), 1)#os.SEEK_CUR)

	return read_bytes

def rewind( file_descriptor, bytes_to_rewind=None ) :
	'''
	rewind the file descriptor by the amount of bytes_to_rewind.

	simple wrapper for the seek function of python
	file object.

	This function is useful for the case 
	that requires rewinding the last read line.

	Usage:
	last_read_line = fp.readline()
	rewind( fp, len(last_read_line) )
	'''
	file_descriptor.seek( -bytes_to_rewind, 1 )#os.SEEK_CUR )

def is_eof( file_descriptor ) :
	'''
	returns True if the file_decriptor is reached EOF.
	otherwise returns False.

	Note that this function simply use peek function to see
	if there are remaining bytes to read more or not.
	'''
	peeked = peek( file_descriptor )

	if peeked :
		return False
	else :
		return True

def find_command_in_path( cmd, path_list=None ) :
	'''
	returns the command combined with the first system default path
	containing the given cmd.
	If the search fails, it returns None.
	'''

	base, cmd = os.path.split(cmd) #split dir and filename

	if path_list == None :
		path_list = os.getenv( 'PATH' ).strip().split(':')

	for path in path_list :
		cmd_path = os.path.join( path, cmd )
		if os.path.exists( cmd_path  ) :
			return cmd_path
	else :
		return None
		

##modified from
## {{{ http://code.activestate.com/recipes/577187/ (r9)



"""
def is_sorted( lst ):
    '''
    Returns True for increasing order.
    Otherwise, returns False.
    '''
    for a, b in izip(lst[:-1], lst[1:]):
        if a >= b :
            return False
    return True
"""
def is_sorted( lst) :
	for i in range( len(lst)-1 ) :
		if lst[i] >= lst[i+1] :
			return False
	return True


def find_all_indices( lst, el ) :
	'''
	Returns a list of indices
	of all elements same as el.
	'''
	ret = []
	try :
		idx = lst.index(el)
		while True :
			ret.append( idx )
			idx = lst.index(el, idx+1)
	except ValueError :
		pass

	return ret

import re

class TagChecker :
	def __init__( self, check_histag=True, histag_start_position=30, user_tag_regex=[], user_tag_start_positions=[] ) :
		'''
		This is simple sequence Tag checking class.

		By default, the following types of histindine tag is checked.
		special characters are following the convention in perl regular expressions.

		HHHHH* within 30 AA from N-term or C-term.

		HHHHH.{0,10}ENLYFQ
		ENLYFQ.{0,10}HHHHH
		within 30 AA from N-term or C-term

		HHHHH.{0,10}LVPRGS
		LVPRGS.{0,10}HHHHH
		within 30 AA from N-term or C-term 

		Additional tag can be added using regular expressions.
		Also disable the default HisTag checking by check_histage=False.
		'''

		self.tags = []
		self.tag_positions = []


		if check_histag == True :
			#add histidine tags
			self.tags.append( re.compile( "HHHHH.{0,10}ENLYFQ" ) )
			self.tag_positions.append( histag_start_position )

			self.tags.append( re.compile( "ENLYFQ.{0,10}HHHHH" ) )
			self.tag_positions.append( -histag_start_position )
	
			self.tags.append( re.compile( "HHHHH.{0,10}LVPRGS" ) )
			self.tag_positions.append( histag_start_position )

			self.tags.append( re.compile( "LVPRGS.{0,10}HHHHH" ) )
			self.tag_positions.append( -histag_start_position )

			self.tags.append( re.compile( "HHHHH*" ) )
			self.tag_positions.append( histag_start_position )

			self.tags.append( re.compile( "HHHHH*" ) )
			self.tag_positions.append( -histag_start_position )


		for tag, starting_position in zip(user_tag_regex, user_tag_start_positions) :
			#need to check if they are regular expressions.
			if tag :
				self.tags.append( tag )
				self.tag_positions.append( starting_position )


	def check( self, protein_sequence ) :
		'''
		Checks if the given protein sequence contains histidine tag
		or other user defined tag sequences.

		It returns the match object (which can be equivalent to true).
		Otherwise, it returns None.
		'''
		matches = []
		for tag, position in zip(self.tags, self.tag_positions) :
			match = tag.search( protein_sequence )
			if match : #match found!!
				if position > 0 and match.start() < position :
					matches.append( match )
				elif position < 0 and match.start()- len(protein_sequence) > position :
					matches.append( match )
			
		else :
			return matches
					

def compare_sequences( seq1, seq2 ) :
	'''
	returns True if the two sequences are same.
	Otherwise, returns False.

	Note that this function check rather promiscuously
	considering X is matching to every other Alphabets.

	Also this function requires both seq1 and seq2 
	are in Upper cases.
	'''
	if seq1 == seq2 :
		#fastest checking.
		return True
	elif len(seq1) == len(seq2 ) :
		#when lengths are Same there are chances.
		for a, b in zip( seq1, seq2 ) :
			if a == b :
				continue
			elif a == 'X' or b == 'X' :
				continue
			else :
				return False
		return True

	else :
		return False


def get_cpu_number() :
	'''
	returns number of CPUs.

	Note that this method is written only linux machine!
	Since it relies on the /proc/cpuinfo.
	'''
	fp = open( '/proc/cpuinfo' )
	fp.read()
	return fp.count( 'processor\t:' )


if __name__ == '__main__' :
	print(parse_profile_filename( "junk/12asB.1.aln" ))
	tagchecker = TagChecker()
	match= tagchecker.check( "MHHHHHHHABDECCDE" )
	print(match)
	print(match.group(), match.span(), match.expand("YYYY"))
