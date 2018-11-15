'''
This module contains classes related 
generic pairwise alignment parsers and objects.

Pairwise alignments are contained in a hierarchy
grouped by a query then a hit and finally a method.
'''
import os, sys, io, bz2
from evdblib.Utils import rewind
from evdblib.Utils import compare_sequences
from evdblib.Utils.Algorithms.SequenceAlignment import align_sequences

verbose = 0

##################################
#Record object
##################################
class PairwiseAlignmentMethodRecord :
        def __init__(self, id1=None, id2=None, 
		method_name=None,
		alignment1=None, alignment2=None, 
                start_position1=0, start_position2=0,
                additional_scores=None, raw_score=None, norm_score=None ) :
                '''
                Container class for a pairwise alignment.
                To accomodate local alignments, full length sequences or 
                starting positions might also be contained.

                Also the similarity also can be saved as raw_score, normalized_score,
                or additional scores as a dictionary.
                '''

		self.method_name = method_name
                self.alignment1 = alignment1
                self.alignment2 = alignment2

                self.start_position1 = start_position1
                self.start_position2 = start_position2
        

		#This is accounting purpose only
		#I put these variables for checking..
		if self.alignment1 :
			self.seq1 = ''.join([a for a in self.alignment1 if a.isalpha()]).upper() 
		else :
			self.seq1 = None

		if self.alignment2 :
			self.seq2 = ''.join([a for a in self.alignment2 if a.isalpha()]).upper() 
		else :
			self.seq2 = None

                self.raw_score = raw_score
                self.normalized_score = norm_score        

		self.id1 = id1
		self.id2 = id2

                #optional info
                self.additional_scores = additional_scores

		#contants for formatting
		self.method_mark = '#'
		self.score_mark = 'scores_from_program:'
		self.alternative_score_mark = 'scores from program:'
		self.null_record_mark = '\n'
		self.additional_score_mark = '#'
		self.hit_record_start = '##'
		self.hit_record_end = '--\n'

	def __str__( self ) :
		fp = io.StringIO()
		self.write( fp ) 
		return fp.getvalue()

	def set_ids( self, id1, id2 ) :
		self.id1 = id1
		self.id2 = id2

	def get_sequence_length_in_alignment1( self ) :
		'''
		returns query sequence length in the aligned region.
		'''
		if self.alignment1 :
			return len([ a for a in self.alignment1 if a.isalpha() ])
		else :
			return 0

	def get_sequence_length_in_alignment2( self ) :
		'''
		returns hit sequence length in the aligned region.
		'''
		return len([ a for a in self.alignment2 if a.isalpha() ])

	def get_equivalent_map( self, use_case=1 ) :
		'''
		returns equivalent positions in both alignment1 and alignment2,
		mapped to the original position in seq1 and seq2.

		It is quite useful to use this equivalence.
		'''
		i = self.start_position1-1
		j = self.start_position2-1
		eq1 = []
		eq2 = []
		for a, b in zip(self.alignment1, self.alignment2) :
			if a.isalpha() :
				i += 1

			if b.isalpha() :
				j += 1

			if use_case :
				if a.isupper() and b.isupper() :
					eq1.append(i)
					eq2.append(j)
			else :
				if a.isalpha() and b.isalpha() :
					eq1.append(i)
					eq2.append(j)
					
		return [eq1, eq2]


	def post_process( self, sec1, sec2 ) :
		'''
		modifies current alignments by making 
		the positions of incompatible secondary structures
		to be lower cases.

		It returns number of modified residues.
		'''
		modification = 0

		if self.is_null() :
			return modification

		fp1 = io.StringIO() 
		fp1.write( self.alignment1 )
		fp2 = io.StringIO()
		fp2.write( self.alignment2 )

		i = self.start_position1-1
		j = self.start_position2-1

		null_record = True
		for n, (a, b) in enumerate( zip(self.alignment1, self.alignment2) ) :
			if a.isalpha() :
				i += 1

			if b.isalpha() :
				j += 1

			if a.isupper() and b.isupper() :
				if (sec1[i] == sec2[j] ) or sec1[i] == 'L' or sec2[j] == 'L' :
					if null_record :
						null_record = False
				else :
					fp1.seek(n)
					fp1.write( a.lower() )

					fp2.seek(n)
					fp2.write( b.lower() )

					modification += 1

		if modification :
			self.alignment1 = fp1.getvalue()
			self.alignment2 = fp2.getvalue()

		if null_record :
			self.alignment1 = None
			self.alignment2 = None

		return modification


					
        def get_normalized_score( self ) :
		if self.is_null() :
			raise ValueError( "Null record do not have a normalized score." )
                return self.normalized_score

        def get_raw_score( self ) :
		if self.is_null() :
			raise ValueError( "Null record do not have a raw score." )
                return self.raw_score

	def is_null( self ) :
		if self.alignment1 :
			return False
		else :
			return True

	def write( self, fp ) :
		'''
		make output of the class into the stream fp.
		'''
		if not self.method_name :
			return 

		if self.alignment1 and self.alignment2 :
			pass
		else :
			#print >>fp, "%s %s\n" % (self.method_mark,self.method_name)
			fp.write( "%s %s\n%s" % (self.method_mark, self.method_name, self.null_record_mark ) )
			return
			
		print(self.method_mark, self.method_name, file=fp)

		if self.additional_scores :
			print(self.score_mark, self.raw_score, self.normalized_score, self.additional_score_mark, self.additional_scores, file=fp)
		else :
			print(self.score_mark, self.raw_score, self.normalized_score, file=fp)
			
		print('%s\t%s'% (self.start_position1, self.alignment1), file=fp)
		print('%s\t%s'% (self.start_position2, self.alignment2), file=fp)


	def parse( self, fp ) :
		line = fp.readline()
		if line == self.hit_record_end :
			rewind( fp, len(line) )
			raise NoMoreMethodRecord()

		while not line.startswith( self.method_mark )  :
			line = fp.readline()
			if not line :
				raise NoMoreMethodRecord( "EOF reached within a record." )
			if line == self.hit_record_end or line.startswith( self.hit_record_start ) :
				rewind( fp, len(line) )
				raise NoMoreMethodRecord()

		#parse the method header line
		columns = line.split()
		if len( columns ) == 2 :
			self.method_name = columns[1]
		elif len( columns ) == 3 :
			self.method_name = columns[1]
			self.real_iteration = columns[2]

		else :
			raise PairwiseAlignmentFormatError( "Method header format is wrong.%s"%line )

		#parse the score line
		line = fp.readline()
		if not line  :
			raise PairwiseAlignmentFormatError( "EOF reached where the score line is expected!" )

		if line.startswith( self.score_mark ) :
			score_line = line[len(self.score_mark):]

		elif line.startswith( self.alternative_score_mark ) :
			score_line = line[len(self.alternative_score_mark):]
		
		elif line == self.null_record_mark :
			#stops parsing if we got the null record!!
			#this is a placeholder for marking a tried run of 
			#alignment generation.
			#but no significant alignment were found!
			return
		else :
			raise PairwiseAlignmentFormatError( "Score line expected:%s"%line)

		if self.additional_score_mark in score_line :
			#when additional scores are there
			try :
				two_score_line, additional_score_line = score_line.split( '#' )
			except :
				raise PairwiseAlignmentFormatError( "Score line with additional score is expected:%s"%score_line )
		else :
			#when only two column of raw and normaalized scores are expected.
			two_score_line = score_line
			additional_score_line = None

		#parse the two score line
		try :
			self.raw_score, self.normalized_score = two_score_line.split()
			self.raw_score = float(self.raw_score)
			self.normalized_score = float(self.normalized_score)
		except :
			raise PairwiseAlignmentFormatError( "Two score line has problem:%s"%two_score_line )

		
		#setting additional scores
		if additional_score_line :
			self.additional_scores = eval( additional_score_line )
			
		#parse the two lines of alignments
		alignment_line = fp.readline()
		try :
			start_pos, aln = alignment_line.split()
			self.start_position1 = int(start_pos)
			self.alignment1 = aln
		except :
			raise PairwiseAlignmentFormatError( "Error in parsing the alignment line1:%s"%alignment_line )

		alignment_line = fp.readline()
		try :
			start_pos, aln = alignment_line.split()
			self.start_position2 = int(start_pos)
			self.alignment2 = aln
		except :
			raise PairwiseAlignmentFormatError( "Error in parsing the alignment line2:%s"%alignment_line )

		#fix some stupid FAST record has None as alignment!
		if self.start_position1 == 0 and self.start_position2 == 0 and self.alignment1== 'None' and self.alignment2 == 'None' :
			#need to be reset!
			self.alignment1 = None
			self.alignment2 = None
			return
		

		self.seq1 = ''.join([a for a in self.alignment1 if a.isalpha()]).upper() 
		self.seq2 = ''.join([a for a in self.alignment2 if a.isalpha()]).upper() 


	def is_query_sequence_correct( self, seq ) :
		'''
		returns True if the query sequence is 
		same as the seq.
		'''
		if self.alignment1==None :
			return True

		seq = seq[self.start_position1: self.start_position1+len(self.seq1)]
		if compare_sequences( seq, self.seq1 ) :
			return True
		else :
			if verbose :
				print("Query is not correct..")
				print(self.id1, self.id2, self.method_name, ":")
				print(seq) 
				print(self.seq1)
			return False


	def is_hit_sequence_correct( self, seq ) :
		'''
		returns True if the hit sequence is 
		same as the seq.
		'''
		if self.alignment2==None :
			return True

		if self.alignment1 == self.alignment2 == 'None' :
			return True

		seq = seq[self.start_position2: self.start_position2+len(self.seq2)]
		if compare_sequences( seq, self.seq2 ) :
			return True
		else :
			if verbose :
				print("Hit is not correct..")
				print(self.id1, self.id2, self.method_name, ":")
				print(seq) 
				print(self.seq2)

			return False


		
	def is_query_sequence_consistent( self, record ) :
		'''
		returns True if the two method records 
		query lines are consistent (meaning the overlapping part is same )
		Otherwise returns False.
		'''
		if self.alignment1 and record.alignment1 :
			pass
		else :
			return True

		#origianlly I thought the '-''s are the only gap character.
		#but figure compass is using "=", too.
		#seq1 = self.alignment1.replace( "-", '' )
		#seq1 = ''.join([a for a in self.alignment1 if a.isalpha()])
		#seq2 = record.alignment1.replace( '-', '' )
		#seq2 = ''.join([a for a in record.alignment1 if a.isalpha()])
		#the above regularization part is now moved to initializer.

		seq1 = self.seq1
		seq2 = record.seq1

		if self.start_position1 < record.start_position1 :
			offset = record.start_position1 - self.start_position1 
			seq1 = seq1[offset:]
		else :
			offset = self.start_position1 - record.start_position1
			seq2 = seq2[offset:]

		minlen = min( len(seq1), len(seq2) )

		if compare_sequences( seq1[:minlen], seq2[:minlen] ) :
			return True
		else :
			if verbose :
				print(self.method_name, record.method_name)
				print(seq1)
				print(seq2)

			return False

	def compare_query_sequences_detail( self, record, print_flag=True ) :
		'''
		prints out the alignment between two records.
		This function is for designed for printing out
		debugging information.
		'''
		if self.alignment1 and record.alignment1 :
			pass
		else :
			if print_flag :
				print("Error: Missing alignments!", self.method_name, record.method_name)

		#origianlly I thought the '-''s are the only gap character.
		#but figure compass is using "=", too.
		#seq1 = self.alignment1.replace( "-", '' )
		#seq1 = ''.join([a for a in self.alignment1 if a.isalpha()])
		#seq2 = record.alignment1.replace( '-', '' )
		#seq2 = ''.join([a for a in record.alignment1 if a.isalpha()])
		#the above regularization part is now moved to initializer.

		seq1 = self.seq1
		seq2 = record.seq1

		if self.start_position1 < record.start_position1 :
			offset = record.start_position1 - self.start_position1 
			seq1 = seq1[offset:]
		else :
			offset = self.start_position1 - record.start_position1
			seq2 = seq2[offset:]

		minlen = min( len(seq1), len(seq2) )

		aln1, aln2, score = align_sequences( seq1[:minlen], seq2[:minlen] )

		if print_flag :

			print("*"*3, self.id1, self.id2, '*'*3) 
			print(self.method_name + '\t'+ aln1)
			print(record.method_name + '\t' + aln2)


	def compare_hit_sequences_detail( self, record, print_flag=True ) :
		'''
		prints out the alignment between two records.
		This function is for designed for printing out
		debugging information.
		'''
		if self.alignment2 and record.alignment2 :
			pass
		else :
			if print_flag :
				print("Error: Missing alignments!", self.method_name, record.method_name)

		#origianlly I thought the '-''s are the only gap character.
		#but figure compass is using "=", too.
		#seq1 = self.alignment1.replace( "-", '' )
		#seq1 = ''.join([a for a in self.alignment1 if a.isalpha()])
		#seq2 = record.alignment1.replace( '-', '' )
		#seq2 = ''.join([a for a in record.alignment1 if a.isalpha()])
		#the above regularization part is now moved to initializer.

		seq1 = self.seq2
		seq2 = record.seq2

		if self.start_position2 < record.start_position2 :
			offset = record.start_position2 - self.start_position2
			seq1 = seq1[offset:]
		else :
			offset = self.start_position2 - record.start_position2
			seq2 = seq2[offset:]

		minlen = min( len(seq1), len(seq2) )

		aln1, aln2, score = align_sequences( seq1[:minlen], seq2[:minlen] )

		if print_flag :

			print("*"*3, self.id1, self.id2, '*'*3)
			print(self.method_name + '\t'+ aln1)
			print(record.method_name + '\t' + aln2)


	
	def is_hit_sequence_consistent( self, record ) :
		'''
		returns True if the two method records 
		hit lines are consistent (meaning the overlapping part is same )
		Otherwise returns False.
		'''
		if self.alignment2 and record.alignment2 :
			pass
		else :
			return True

		#origianlly I thought the '-''s are the only gap character.
		#but figure compass is using "=", too.
		#seq1 = self.alignment1.replace( "-", '' )
		#seq1 = ''.join([a for a in self.alignment1 if a.isalpha()])
		#seq2 = record.alignment1.replace( '-', '' )
		#seq2 = ''.join([a for a in record.alignment1 if a.isalpha()])
		#the above regularization part is now moved to initializer.

		seq1 = self.seq2
		seq2 = record.seq2

		if self.start_position2 < record.start_position2 :
			offset = record.start_position2 - self.start_position2
			seq1 = seq1[offset:]
		else :
			offset = self.start_position2 - record.start_position2
			seq2 = seq2[offset:]

		minlen = min( len(seq1), len(seq2) )

		if compare_sequences( seq1[:minlen], seq2[:minlen] ) :
			return True
		else :
			if verbose :
				print(self.method_name, record.method_name)
				print(seq1)
				print(seq2)

			return False

	def count_all( self ) :
		if self.alignment1 and self.alignment2 :
			return 1
		else :
			return 0

		
class PairwiseAlignmentHitRecord :
	def __init__( self, id1=None, id2=None, sequence2=None, fasta2=None ) :
		'''
		Container class for pairwise alignments
		modelled by PairwiseAlignmentMethodRecords class.

		id1 is the identifier for a protein (or domain) 1 AKA query.
		id2 is the identifier for a protein (or domain) 2 AKA hit.

		sequence2 is the sequence of the hit.
		Sequence2 is optional.

		fasta2 is the fasta object of the hit.
		fasta2 is optional.

		strict 
		'''
		self.id1 = id1
		self.id2 = id2
		
		self.sequence2 = sequence2
		self.fasta2 = fasta2

		self.method2record = {}

		#constants for writing and reading the format
		self.record_start = '##'
		self.record_end = '--\n'

	def __len__( self ) :
		return len( self.method2record )

	def __delitem__(self, method_name) :
		del self.method2record[ method_name ]

	def __getitem__( self, method ) :
		return self.method2record.get( method )

	def __iter__( self ) :
		return self.method2record.__iter__()

	def __contains__( self, method ) :
		return method in self.method2record

	def __str__( self ) :
		fp = io.StringIO()
		self.write( fp ) 
		return fp.getvalue()


	def get_method_records( self, method_name, exact=True ) :
		'''
		Returns a tuple of method records that matches
		to the given method_name.

		If exact flag is True, only exactly matching
		record will be returned.
		Otherwise, all method records have method name
		starting with the method_name will be returned.

		Note that if no matching method recrod is found,
		None will be returned.
		'''
		if exact :
			if methodname == method_name :
				return ( self[mehod_name], )
			else :
				return None

		method_records = []
		for methodname in self :
			if methodname.startswith( method_name ) :
				method_records.append( self[methodname] )

		return tuple( method_records )


	def change_id( self, conversion_dictionary, strict=True ) :
		'''
		Changes IDs according to the convresion dictionary.
		If strict is True, the ID missing in the conversion dictionary will
		raise an exception.
		'''
		if self.id1 in conversion_dictionary :
			newid1 = conversion_dictionary[self.id1]
			self.id1 = newid1
		elif strict :
			raise IDNotFoundError( self.id1 )
		else :
			newid1 = self.id1

		if self.id2 in conversion_dictionary :
			newid2 = conversion_dictionary[self.id2]
			self.id2 = newid2
		elif strict :
			raise IDNotFoundError( self.id2 )
		else :
			newid2 = self.id2
		
		for method in self :
			methodrec = self[method]
			methodrec.set_ids( newid1, newid2 )
			


	def parse_method_name_and_iteration( self, method_name ) :
		splitted_list = method_name.split('_', 1)
		if len(splitted_list) == 1 :
			return splitted_list[0], 0
		else :
			return splitted_list[0], int(splitted_list[1])


	def get_highest_iteration_method_record( self, method_name ) :
		'''
		Returns the highest iteration "compass" or "hhsearch"
		record. Basically the digits after an underscore
		will be treated as iteration number, e.g. "compass_1" or 
		"hhsearch_8".

		When the parsing iteration is not an integer,
		ValueError will be raised.
		'''
		method_records = self.get_method_records( method_name, exact=False )
		if len(method_records) == 1:
			return method_records[0]

		else :
			record_list = []
			for methodrec in method_records :
				basename, iteration = self.parse_method_name_and_iteration( methodrec.method_name )
				record_list.append( (iteration, methodrec) )

			iteration, methodrec = max(record_list)
			return methodrec

					
	def update( self, hit_record, overwrite=False ) :
		'''
		Adds the hit records in the query record variable.
		'''
		for methodid in hit_record :
			if methodid in self :
				if overwrite :
					self.method2record[methodid] = hit_record[methodid]
				else :
					if verbose :
						print("New:")
						print(hit_record)
						print("Old:")
						print(self)

					raise AlignmentAddingError( "MethodID %s already exists."%methodid )
			else :
				self.add( hit_record[methodid], overwrite=overwrite )
	
	def add( self, methodrec, overwrite=False ) :
		'''
		adds the method records.
		'''
		method_name = methodrec.method_name
		if not method_name :
			if verbose :
				print(methodrec)
			raise NullMethodNameError( "An alignment for a method should have a defined method name." )
		self.method2record[ method_name ] = methodrec


	def write( self, fp ) :
		'''
		print out the contents of the alignment hit records
		'''
		if not len(self ) :
			return

		#print >>fp, '##', self.id1, self.id2
		print(self.record_start, self.id1, self.id2, file=fp)

		for rec in self.method2record.values() :
			rec.write( fp )

		#print >>fp, '--'
		fp.write( self.record_end )
		

	def parse( self, fp, overwrite=False ) :
		'''
		parse the content in the fp.
		The format is assumed to be the "## --" format.
		'''

		line = fp.readline()

		#finding the initial point
		#this is a little bit loose formatting.
		#I allow the fp actually contains some junk before 
		# the record start maker or "##"
		while not line.startswith( self.record_start ) :
			line = fp.readline()
		
			if not line  :
				raise NoMoreHitRecord()

		try :
			start_marker, id1, id2 = line.split()
		except :
			raise AlignmentFormatError( "Wrong header line:" + line )


		self.id1 = id1
		self.id2 = id2
		
		#read whole record first!
		method_fp = io.StringIO()
		line = fp.readline()
		while line != self.record_end :
			method_fp.write( line )
			line = fp.readline()
			if not line :
				break

		#after hit record is read!
		#put the pointer in the front!
		method_fp.seek(0)
			
		while 1 :
			methodrecord = PairwiseAlignmentMethodRecord()
			try :
				methodrecord.parse( method_fp )
			except NoMoreMethodRecord :
				break

			methodrecord.id1 = self.id1
			methodrecord.id2 = self.id2

			self.add( methodrecord, overwrite=overwrite )

		#line = fp.readline() #should have be read before!
		if line != self.record_end :
			raise AlignmentFormatError( "Wrong end of the hit record.%s"%line )

	def is_query_sequence_correct( self, seq ) :
		'''
		returns True if all methods' query sequences are correct.
		Otherwise returns False.
		'''
		for rec in self.method2record.values() :
			if rec.is_query_sequence_correct( seq )  :
				continue
			else :
				return False
		return True


	def compare_query_sequences_detail( self ) :
		rec_list = list(self.method2record.values())

		for i in range( len(rec_list) ) :
			rec1 = rec_list[i]
			for j in range( i+1, len(rec_list) ) :
				rec2 = rec_list[j]

				if not rec1.is_query_sequence_consistent( rec2 ) :
					rec1.compare_query_sequences_detail( rec2 )
		return True


	def compare_hit_sequences_detail( self ) :
		rec_list = list(self.method2record.values())

		for i in range( len(rec_list) ) :
			rec1 = rec_list[i]
			for j in range( i+1, len(rec_list) ) :
				rec2 = rec_list[j]

				if not rec1.is_hit_sequence_consistent( rec2 ) :
					rec1.compare_hit_sequences_detail( rec2 )
		return True


	def is_hit_sequence_correct( self, seq ) :
		'''
		returns True if all methods' hit sequences are correct.
		Otherwise returns False.
		'''
		for rec in self.method2record.values() :
			if rec.is_hit_sequence_correct( seq )  :
				continue
			else :
				return False
		return True


	def is_query_sequence_consistent( self ) :
		rec_list = list(self.method2record.values())

		for i in range( len(rec_list) ) :
			rec1 = rec_list[i]
			for j in range( i+1, len(rec_list) ) :
				rec2 = rec_list[j]

				if rec1.is_query_sequence_consistent( rec2 ) :
					continue
				else :
					if verbose :
						print("get false from a rec", rec1.id1, rec1.id2, rec1.method_name, rec2.method_name)
					return False

		return True
			

	def is_hit_sequence_consistent( self ) :
		rec_list = list(self.method2record.values())
		
		for i in range( len(rec_list) ) :
			rec1 = rec_list[i]
			for j in range( i+1, len(rec_list) ) :
				rec2 = rec_list[j]

				if rec1.is_hit_sequence_consistent( rec2 ) :
					continue
				else :
					if verbose :
						print("get false from a rec", rec1.id1, rec1.id2, rec1.method_name, rec2.method_name)
					return False
		return True
			

	def filter( self, method_list=None ) :
		'''
		removes the records that does not contain methods 
		not belong to the id and method lists.
		'''
		if method_list == None :
			return

		for methodname in list(self.method2record.keys()) :
			if methodname in method_list :
				#do nothing.
				pass
			else :
				del self[methodname]

	def count_all( self ) :
		number = 0
		for methodname in self :
			number += self[methodname].count_all()
		return number


class PairwiseAlignmentQueryRecord :
	'''
	Container class for Pairwise alignments
	modelled by PariwiseAlignmentHitRecord class.

	This class is the highest level container class
	of paiwsie alignments.
	'''
	def __init__( self, queryid=None ) :
		self.queryid = queryid

		self.hit2record = {}

	def __len__( self ) :
		return len( self.hit2record )
	
	def __delitem__(self, hitid) :
		del self.hit2record[hitid]

	def __getitem__( self, hitid ) :
		return self.hit2record.get( hitid )

	def __iter__( self ) :
		return self.hit2record.__iter__()

	def __contains__( self, hitid ) :
		return hitid in self.hit2record

	def __str__( self ) :
		fp = io.StringIO()
		self.write( fp ) 
		return fp.getvalue()

	def update( self, query_record, overwrite=False ) :
		'''
		Adds the hit records in the query record variable.
		'''
		for hitid in query_record :
			if hitid in self :
				self[hitid].update( query_record[hitid], overwrite=overwrite )
			else :
				self._add_hit_record( query_record[hitid], overwrite=overwrite )
	

	def change_id( self, conversion_dictionary, strict=True ) :
		'''
		Changes IDs according to the convresion dictionary.
		If strict is True, the ID missing in the conversion dictionary will
		raise an exception.
		'''
		for hitid in self :
			hitrec = self[hitid]
			if hitid in conversion_dictionary :
				new_hitid = conversion_dictionary[hitid]
				hitrec = self[hitid]

				del self[hitid]
				self.add(hitrec)

			elif strict :
				raise IDNotFoundError( hitid )

			hitrec.change_id( conversion_dictionary, strict=strict )



	def add( self, record, overwrite=False ) :
		'''
		Adds the record into the class.
		If the class is an instance of PairwiseAlignmentHitRecord,
		it will try to add it as is.
		If the class is an instance of PairwiseAlignmentMethodRecord,
		this function will try to add into an existing hit record.
		'''
		if isinstance( record, PairwiseAlignmentQueryRecord ) :
			self.update( record, overwrite )
		if isinstance( record, PairwiseAlignmentHitRecord )  :
			self._add_hit_record( record, overwrite )
		else :
			if record.id2 in self :
				hit_record = self[record.id2]
				hit_record.add( record )
			else :
				hit_record = PairwiseAlignmentHitRecord( id1=record.id1, id2=record.id2 )
				hit_record.add( record, overwrite=overwrite )
				self._add_hit_record( hit_record, overwrite )



	def _add_hit_record( self, record, overwrite=False ) : 
		'''
		Adds a hit record.
		If the overwrite value is true, 
		the pre-existing record will be overwritten.
		'''
		if record.id2 in self :
			self[record.id2].update(record, overwrite=overwrite)
		else :
			self.hit2record[record.id2] = record

			
	def write( self, fp ) :
		for rec in self.hit2record.values() :
			rec.write( fp )
			

	def is_query_sequence_correct( self, seq ) :
		'''
		returns True if all hit records query sequences are correct.
		Otherwise returns False.
		'''
		for rec in self.hit2record.values() :
			if rec.is_query_sequence_correct( seq )  :
				continue
			else :
				return False
		return True

	def is_hit_sequence_correct( self, seq_dic ) :
		'''
		returns True if all hit records hit sequences are correct.
		Otherwise returns False.

		seq_dic is a dictionary of sequences 
		with hit id keys.
		'''
		for hit_id, rec in self.hit2record.items() :
			if rec.is_hit_sequence_correct( seq_dic[hit_id] )  :
				continue
			else :
				return False
		return True


	def is_query_sequence_consistent( self ) :
		'''
		returns True if all hit records query sequences are consistent.
		Otherwise returns False.
		'''
		for rec in self.hit2record.values() :
			if rec.is_query_sequence_consistent()  :
				continue
			else :
				if verbose :
					print("got false from a rec", rec.id1, rec.id2)
				return False
		return True


	def is_hit_sequence_consistent( self ) :
		'''
		returns True if all hit records hit sequences are consistent.
		Otherwise returns False.
		'''
		for rec in self.hit2record.values() :
			if rec.is_hit_sequence_consistent()  :
				continue
			else :
				if verbose :
					print("got false from a rec", rec.id1, rec.id2)
				return False
		return True

	
	def filter( self, id_list=None, method_list=None ) :
		'''
		removes the records that does not contain ids and methods 
		not belong to the id and method lists.
		'''
		for hitid in list(self.hit2record.keys()): #need to copy!!
			hit_record = self[hitid]
			if hitid in id_list :
				hit_record.filter( method_list )

				if len(hit_record ) == 0:
					del self[hitid]

			elif id_list == None :
				hist_record.filter( id_list, method_list )

				if len(hit_record ) == 0:
					del self[hitid]

			else :
				del self[hitid]

	def count( self, queryid, id_list, methodname, strict=True ) :
		number = 0
		for hitid in list(self.hit2record.keys()) :
			if hitid in id_list :
				hit_record = self[hitid]
				if strict :
					if hit_record and hit_record[methodname] and hit_record[methodname].alignment1:
						number+=1
				else :
					if hit_record and hit_record[methodname] :
						number+=1
		return number
	
	def count_all( self ) :
		'''
		Count all non-null method records.
		'''
		number = 0
		for hitid in self :
			number += self[hitid].count_all()

		return number
			


class PairwiseAlignmentRecords :
	'''
	Container class for Pairwise alignments
	modelled by PariwiseAlignmentQueryRecord class.

	This class is the highest level container class
	of paiwsie alignments.
	'''
	def __init__( self, alignment_filename=None ) :

		self.query2record = {}

		if alignment_filename :
			self.parse( alignment_filename )

	def __len__( self ) :
		return len( self.query2record )

	def __getitem__( self, queryid ) :
		return self.query2record.get( queryid )

	def __iter__( self ) :
		return self.query2record.__iter__()

	def __delitem__(self, queryid) :
		del self.query2record[queryid]

	def __contains__( self, queryid ) :
		return queryid in self.query2record

	def __str__( self ) :
		fp = io.StringIO()
		self.write( fp ) 
		return fp.getvalue()

	def get_query_ids( self ) :
		return list(self.query2record.keys())

	def change_id( self, conversion_dictionary, strict=True ) :
		'''
		Changes IDs according to the convresion dictionary.
		If strict is True, the ID missing in the conversion dictionary will
		raise an exception.
		'''
		for queryid in self :
			queryrec = self[queryid]
			if queryid in conversion_dictionary :
				new_queryid = conversion_dictionary[queryid]

				del self[queryid]
				self.add(queryrec)

			elif strict :
				raise IDNotFoundError( queryid )

			queryrec.change_id( conversion_dictionary, strict=strict )


	def update( self, parecords, overwrite=False ) :
		'''
		Adds the query records in the parecords variable.
		'''
		for queryid in parecords :
			if queryid in self :
				self[queryid].update( parecords[queryid], overwrite=overwrite )
			else :
				self._add_query_record( parecords[queryid] )
			

	def add( self, record, overwrite=False ) :
		'''
		Adds the record into the class.
		If the class is an instance of PairwiseAlignmentHitRecord,
		it will try to add it as is.
		If the class is an instance of PairwiseAlignmentMethodRecord,
		this function will try to add into an existing hit record.
		'''

		if verbose :
			print("Adding New Record in the PairwiseAlignmentRecords class")

		if isinstance( record, PairwiseAlignmentRecords ) :
			if verbose :
				print("The new record type is PairwsieAlignmentRecords")
				print("starting to update...")

			self.update( record, overwrite )

		elif isinstance( record, PairwiseAlignmentQueryRecord )  :
			if verbose :
				print("The new record type is PairwsieAlignmentQueryRecords")
				print("starting to add the Query Record...")

			self._add_query_record( record, overwrite )

		elif isinstance( record, PairwiseAlignmentHitRecord ) :
			if verbose :
				print("The new record type is PairwsieAlignmentHitRecords")
				print("starting to add the Hit Record...")

			if record.id1 in self :
				if verbose :
					print("The Query Record of the Hit Record is already present.")
					print("adding into the existing query record.")
				query_record = self[record.id1]
				query_record.add( record, overwrite )
			else :
				if verbose :
					print("The Query Record of the Hit Record is not present.")
					print("building new query record")
				query_record = PairwiseAlignmentQueryRecord( queryid=record.id1 )
				query_record.add( record, overwrite )
				self._add_query_record( query_record, overwrite )

		elif isinstance( record, PairwiseAlignmentMethodRecord ) :
			if verbose :
				print("Thie new record type is PairwiseAlignmentMethodRecord")
				print("starting to add the Method record...")

			if record.id1 in self :
				query_record = self[record.id1]
				query_record.add( record, overwrite )
			else :
				query_record = PairwiseAlignmentQueryRecord( queryid=record.id1 )
				query_record.add( record, overwrite )
				self._add_query_record( query_record, overwrite )

		else :
			raise ValueError( record, "Not allowed record type" )


	def _add_query_record( self, record, overwrite=False ) :
		'''
		Adds a query record.
		If the overwrite value is true, 
		the pre-existing record will be overwritten.
		'''
		if record.queryid in self :
			queryrecord = self[ record.queryid ]
			queryrecord.update( record, overwrite=overwrite )
			
		else :
			self.query2record[record.queryid] = record
				

	def save( self, fn, overwrite=False ) :
		'''
		Saves all content in the current alignment records
		into a file, fn.
		'''

		if not overwrite and os.path.exists( fn ) :
			raise IOError( "file already exits!", fn )

		if fn.endswith( ".bz2" ) :
			fp = bz2.BZ2File( fn, 'w' )
		else :
			fp = open( fn, 'w' ) 

		self.write( fp )
		fp.close()
			

	def write( self, fp ) :
		'''
		Writes all contents in the current alignment
		into the file pointer fp.
		'''
		for rec in self.query2record.values() :
			rec.write( fp )
	
	def is_query_sequence_consistent( self ) :
		'''
		returns True if all qeury records query sequences are consistent.
		Otherwise returns False.
		'''
		for rec in self.query2record.values() :
			if rec.is_query_sequence_consistent() :
				continue
			else :
				if verbose :
					print("got false from a rec", rec.queryid)
				return False
		return True

	def is_hit_sequence_consistent( self ) :
		'''
		returns True if all hit records query sequences are consistent.
		Otherwise returns False.
		'''
		for rec in self.query2record.values() :
			if rec.is_hit_sequence_consistent() :
				continue
			else :
				if verbose :
					print("got false from a rec", rec.queryid)
				return False
		return True


	def is_sequence_correct( self, seq_dict ) :
		'''
		returns True if all sequences in the alignments are 
		correct.
		'''
		for queryid, rec in self.query2record.items() :
			if ( rec.is_query_sequence_correct(seq_dict[queryid]) and rec.is_hit_sequence_correct(seq_dict) ) :
				continue
			else :
				return False
		return True


		
	def parse( self, filename=None, fp=None, overwrite=False ) :
		'''
		This function parse pairwise alignment format from the file.
		'''
		if filename :
			if filename.endswith( '.bz2' ) :
				zfp = bz2.BZ2File( filename )
				fp = io.StringIO( zfp.read() )
				zfp.close()
			else :
				fp = open( filename )

		self.fp = fp
		if not self.fp :
			return

		while 1 :
			hitrecord = PairwiseAlignmentHitRecord()
			try :
				hitrecord.parse( fp, overwrite=True )
				#if verbose :
					#print "reading a hit record in alignment class..."
					#print hitrecord
			except NoMoreHitRecord :
				break

			self.add( hitrecord, overwrite )


	def filter( self, id_list=None, method_list=None ) :
		'''
		removes the records that does not contain ids and methods 
		not belong to the id and method lists.

		Note that the filtering only applies to the Hit and Method records.
		'''
		for queryid in list(self.query2record.keys()) : #need to copy
			
			query_record = self[queryid]

			if queryid in id_list :
				query_record.filter( id_list, method_list )
			elif id_list == None :
				query_record.filter( id_list, method_list )
			else :
				del self[queryid]


	def count( self, queryid, id_list, methodname=None, strict=True ) :
		query_record = self[queryid]
		
		if not query_record :
			return 0

		return query_record.count( queryid, id_list, methodname, strict )


	def count_all( self ) :
		number = 0
		for queryid in self :
			number += self[queryid].count_all()
			
		return number


	def save_seperate_queries( self, save_dir, overwrite=False ) :
		'''
		Saves each queries seprerately into the save_dir 
		using the <queryid>.aln format.
		'''

		for query in self :
			save_fn = save_dir + "/" + query + '.aln' 
			if os.path.exists( save_fn ) :
				if overwrite :
					pass
				else :
					print("Cannot write save_fn", save_fn, file=sys.stderr)
					continue
			
			rec = self[query]
			fp = open(save_fn, 'w' )
			rec.write( fp )
			fp.close()

	def count_inconsistencies( self, exclude_set=None, title='', print_flag=True ) :
		'''
		a utility function to help diagnosing differences in query and hit sequences.
		'''
	
		if print_flag :
			print("*"*30)
			print("Start diagnosis:", title)

		count = 0
		for queryid in self :
			if print_flag :
				print("queryid:", queryid)

			queryrec = self[queryid]
			for hitid in queryrec :
				if exclude_set and hitid in exclude_set :
					continue

				hitrec = queryrec[hitid]

				query_consistency = hitrec.is_query_sequence_consistent()
				hit_consistency = hitrec.is_hit_sequence_consistent()
				if not query_consistency or not hit_consistency :
					if print_flag:
						print("############")
						print("#consistency:", query_consistency, end=' ') 
						print(hit_consistency)
						print("############")
						
						if exclude_set and queryid in exclude_set or query_consistency :
							pass
						else :
							print("Query Problem Details:")
							hitrec.compare_query_sequences_detail()
							
						if exclude_set and hitid in exclude_set or hit_consistency :
							pass
						else :
							print("Hit Problem Details:")
							hitrec.compare_hit_sequences_detail()

				
					count += 1
				
		if print_flag :
			print("Total inconsistencies:", count)
			print("End diagnosis:", title) 
			print("*"*30) 

		return count

class IDNotFoundError( Exception ) :
	pass
		
class AlignmentAddingError( Exception ) :
	pass

class PairwiseAlignmentFormatError( Exception ) :
	pass

class NoMoreHitRecord( Exception ) :
	pass

class NoMoreMethodRecord( Exception ) :
	pass
