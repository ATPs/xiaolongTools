'''
This module contains simple parsers not following loaded OOP concept 
but rather focused on fast and simple parse.

Most of the parsers in this module are simple
functions rather than classes.
'''
import os, sys, io

class SimpleParserError( Exception ) :
	'''
	Simple parsing functions 
	raise this exception when they ran into problems.
	'''
	pass


def parse_palsse_output( file, sequence_length=None ) :
	'''
	Palssee output file reader.
	Reads the palssee program output.
	
	If sequence_length is given,
	parsed length can be checked with expected length.
	'''
        atom = 0
        pssd = {}

	if os.path.exists( file ) :
        	fp = open(file)
	else :
		raise SimpleParserError( file, "Not found!" )

        for l in fp :
                if l[:4] == 'ATOM' :
                        atom += 1

                elif l[:5] == 'HELIX' :
                        start = int(l[21:25])
                        end = int(l[33:37])
                        for i in range(start, end+1) :
                                pssd[i] = 'H'

                elif l[:5] == 'SHEET' :
                        start = int(l[22:26])
                        end = int( l[33:37])
                        for i in range( start, end+1) :
                                pssd[i] = 'S'
        fp.close()

	if sequence_length != None :
		if atom == sequence_length :
			pass #Good!
		else :
			raise SimpleParserError( "Input sequence_length does not match with parsed length.", sequence_length, atom )

	pssd_array = io.StringIO()
        for i in range(1, atom+1) :
                if i in pssd :
                        pssd_array.write( pssd[i] )
                else :
                        pssd_array.write( 'L' )


        return pssd_array.getvalue()



def parse_coordinates_from_ca_only_pdb( fn ) :
        '''
        Parses coordinates from a CA Only PDB 
        that corresponds to the sequence 1:1.
        '''
        fp = open( fn )
        coord = []
        for l in fp :
                if ( l.startswith( 'ATOM  ' ) or l.startswith( 'HETATM' ) ) and l[12:16]==' CA ' :
                        coord.append ( [float(l[30:38]), float( l[38:46] ), float( l[46:54] )] )
        return coord


class EOFScoreFile( Exception ) :
	pass


def parse_score_file_record( fp ) :
	'''
	Parses profile score file pointer (or descriptor)
	and returns 
	1. id1
	2. id2
	2. methodname
	3. scores dictionary {'length:':['47','320','327'], 'id:':['6','47','47', '3.21277'],...}

	raises EOFScoreFile at the end of the file.

	Note that, for simplicity in parsing, 
	this function assumes one method (or score set) per pair!
	'''

	#profile scores example
	'''
	bash-3.00$ cat result.selected_alignment | ../../html/profile_scores_stdin 
	## 1dom aaaa
	# dali
	length: 47 320 327
	id: 6 47 47 3.21277
	blosum: -31 241 244 -48.9149
	cpssm: 0 0 83.4929 0
	pssm: 0 0 166.986 0
	dotproduct: 0 0 1141.55 0
	crossproduct: 0 0 -232.818 0
	compass: 0 0 799.679 0
	picasso: 0 0 76.136 0
	pearson: 0 0 47 0
	--
	'''
	l = fp.readline()
	if not l :
		raise EOFScoreFile()

	if l.startswith( '##' ) :
		junk, id1, id2 = l.split()
	else :
		print(l, file=sys.stderr)
		raise SimpleParserError( "Score Format is different from expectations!" )

	pairid = id1 + ' ' + id2 #^## 1dom aaaa: 1dom aaaa is the hit id.
	methodname = ''
	raw_scores = {}
	l = fp.readline()
	while l :
		if l.startswith( '--' ) :
			break
		elif l.startswith( '#' ) :
			methodname = l.split()[1]
		else :
			score_line = l.split()
			if not score_line :
				l = fp.readline()
				continue 
			
			raw_scores[ score_line[0] ] = score_line[1:]

		l = fp.readline()
		

	return id1, id2, methodname, raw_scores
		
def read_svm_scores( score_fp, total_score_number ) :
        '''
        reads svm input scores and annotations from file-like descritor.
        And this function returns a dictionary 

        {'<score_index1>': [<score11>,<score12>,..., <score1N>],
         '<score_index2>': [<score21>, ... <score2N> ]
         ...
         '<score_indexM>': [<scoreM1>, ... <scoreMN> ] }

        M, N is length of feature vector and size of inputs respectively.
        and an annotation list for all input vector.
        '''
        scores = {}
        annotations = []

        for line in score_fp :

		if line.startswith( '#' ) :
			continue

                # process annotation part!
                if '#' in line :
                        score_string, annotation_string = line.split('#', 1)
                else :  
                        score_string = line
                        annotation_string = ''

                score_string = score_string.strip()
                annotation_string = annotation_string.strip()

                #annoation string is expected to be "1abc 2edf".
                annotations.append( annotation_string )

                l = score_string.split()

                #checking integrity
                if not l :
                        continue

                if l[0] == '+1' :
                        pass
                else :  
                        print("warning: possibly wrong line, skip this line", file=sys.stderr)
                        print(line[:-1], file=sys.stderr)
                        continue

                for name_and_score in l[1:]:
                        name, score = name_and_score.split(":")
			if name not in scores :
				scores[name] = []

                        scores[name].append( float(score) )

        return scores, annotations

class SavingSVMScoresError( Exception ) :
	pass

def write_svm_scores( scores, annotations=[], save_fp=None, sort=True ) :
	'''
	Writes SVM scores into a file-like descriptor.

	If no save_fp is given,
	cStringIO.StringIO object will be returned.
	'''

	if save_fp == None :
		save_fp = io.StringIO()

	#checking if all of feactures has same length in scores.
	#scores is a dictionary of score list.
	lengths = [ len(score_list) for score_list in scores.values() ]
	if len( set(lengths) ) == 1:
		pass #good!
	else :
		raise SavingSVMScoresError( "Feature lengths are not same!", lengths )
	feature_length = lengths[0]

	for i, annotation in enumerate(annotations) :

		print('+1', end=' ', file=save_fp)
		integer_pairs = []
		non_integer_pairs = []
		for name, score_list in scores.items() :
			try : integer_pairs.append( (int(name), score_list[i]) )
			except ValueError: non_integer_pairs.append( (name, score_list[i]) )
			
		for name, score in sorted(integer_pairs) :
			print("%s:%5f"%(name,score), end=' ', file=save_fp)
		for name, score in sorted(non_integer_pairs) :
			print("%s:%5f"%(name,score), end=' ', file=save_fp)
		
		print('#', annotation, file=save_fp)

	annotation = ''
	for i in range(len(annotations),feature_length) :
		print('+1', end=' ', file=save_fp)

		integer_pairs = []
		non_integer_pairs = []
		for name, score_list in scores.items() :
			try : integer_pairs.append( (int(name), score_list[i]) )
			except ValueError: non_integer_pairs.append( (name, score_list[i]) )
			
		for name, score in sorted(integer_pairs) :
			print("%s:%5f"%(name,score), end=' ', file=save_fp)
		for name, score in sorted(non_integer_pairs) :
			print("%s:%5f"%(name,score), end=' ', file=save_fp)

		print('#', annotation, file=save_fp)

	return save_fp
