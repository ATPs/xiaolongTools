import copy
from evdblib.Utils.Parsers.BLAST import *
from io import StringIO

verbose=0

class Alignment:
	
	def __init__ ( self, offset=0 ) :
		self.query_start_index = 0
		self.query_end_index = 0
		self.hit_start_index = 0
		self.hit_end_index = 0
		
		self.query = None
		self.hit = None
		self.evalue = 0.0
		self.bitscore = 0.0
		self.score = 0
		self.identity = 0.0
		self.positive = 0.0
	
		self.hsp_start = hsp_start #' Score ='
		self.query_line_start = 'Query:'
		self.hit_line_start = 'Sbjct:'
		self.modified = 0
		self.overlapping_mark = '='
		self.orignal_hit = ''
		self.index_offset = offset #hit index offset


	def is_short( self, cutoff=5 ) :
		'''
		Returns True when the alignment is shorter than the cutoff.
		Otherwise, returns False.
		'''
		if self.modified :
			hitlength = len( [ a for a in self.hit if a.isupper() ] )
			if hitlength >= cutoff :
				return False
			else :
				return True
		else :
			if self.query_end_index - self.query_start_index + 1 >= cutoff:
				return False
			else :
				return True
				

	def __str__ (self):
		fp = StringIO()
		fp.write( "Query (%4s-%4s): %s\n"%( self.query_start_index, self.query_end_index, self.query) )
		fp.write( "Hit   (%4s-%4s): %s\n"%(self.hit_start_index, self.hit_end_index, self.hit)  )
		return fp.getvalue()

	def to_string( self, query=None, hit=None, query_start_index=None, query_end_index=None, hit_start_index=None,hit_end_index=None ):
		
		fp = StringIO()
		fp.write( "Query (%4s-%4s): %s\n"%( query_start_index, query_end_index, self.query) )
		fp.write( "Hit   (%4s-%4s): %s\n"%(hit_start_index, hit_end_index, self.hit)  )
		return fp.getvalue()

	def write_hit_fasta( self, fp, gi, annot='' ) :
		'''
		Writes FASTA record of Hit sequence.
		'''
		hit = self.hit.replace( '-', '' )
		hit = hit.replace( 'U', 'X' )
		hit = hit.replace( 'u', 'x' )

		str = ">gir|%s_%d..%d|%s\n%s\n" % (gi, self.hit_start_index, self.hit_end_index, annot, hit )
		fp.write( str )


	def write_a3m_fasta( self, fp, header, queryseq, remove_gap=False ) :
		'''
		Writes A3M FASTA record of Hit sequence
		'''
		if self.is_short() :
			return 

		hit_outstream = StringIO()

		#writing the leading -
		hit_outstream.write( (self.query_start_index)*'-' )
		if remove_gap :
			for a,b in zip(self.query, self.hit) :
				if a.isalpha() :
					hit_outstream.write( b )
		else :
			for a, b in zip( self.query, self.hit ) :
				if a.isalpha() :
					hit_outstream.write( b )
				else :
					hit_outstream.write( b.lower() )
		hit_outstream.write( (len(queryseq)-self.query_end_index)*'-' )

		if len([a for a in hit_outstream.getvalue() if a.isupper() or (not a.isalpha())])!= len(queryseq) :
			raise Exception( "Built string length and Query sequence length does not match!", hit_outstream.getvalue(), queryseq )

		print(header, file=fp)
		print(hit_outstream.getvalue(), file=fp)
		
		
	def parse_a3m_fasta( self, hitfasta, queryfasta ) :
		'''
		parse a3m fastas and fill in alignment content.
		'''
		self.query_start_index, self.query_end_index = hitfasta.parse_query_range( queryfasta.sequence )
		self.hit_start_index, self.hit_end_index = hitfasta.parse_hit_range()
		self.query, self.hit = hitfasta.build_query_and_hit_sequences( queryfasta.sequence, self.query_start_index, self.query_end_index )

		return 0

	
	def parse( self, fp ) :
		'''
		Returns 0 for normal end 
		Returns -2 for PSI-BLAST end
		Returns 1 for successful parsing of BlastHSP
		Returns -1 for non-normal end or errorneous parsing.
		'''
		line = fp.readline()
		if verbose :
			print('@@ BlastHSP.parse() first line:', line)
		while not line.startswith(self.hsp_start) :
			if line.startswith( psiblast_end_of_result_marker ) :
				fp.seek(-len(line), 1 )
				if verbose :
					print("End of BLAST result. @HSP paring", line)
				raise NoMoreIteration( "Alignment parsing reached at the end of BLAST result!" )
				return -2
					
			elif line.startswith(psiblast_start) :
				#fp.seek(-len(line), 1 ) #rewind one line back for the parameters for blast run
				if verbose :
					print("End of BLAST iteration reached! @HSP paring:", line) 
				raise NoMoreHit( "Alignment parsing reached at the end of Iteration")
				return -2

			line = fp.readline()
			if verbose :
				print('@@ BlastHSP.parse() skipping lines:', line)
			
			if not line :
				raise NoMoreAlignment("All Alingments in the HIT was read!")
				return 0
		if verbose :
			print('@@ BlastHSP.parse() after initial line check', line)
			
		score_line1 = line.split()
		score_line2 = fp.readline().split()
		self.bitscore = float( score_line1[2] )
		self.score = float( score_line1[4][1:-2] )

		if score_line1[7][0] == 'e' :
			score_line1[7] = '1' + score_line1[7]
		self.evalue = float(score_line1[7][:-1])
		n1, n2 = score_line2[2].split('/')
		self.identity = float( int(n1)*1.0/int(n2) )
		n1, n2 = score_line2[2].split('/')
		self.positive = float( int(n1)*1.0/int(n2) )
		#push aray one line
		line = fp.readline() #white line
		
		#query line
		query_line = fp.readline()
		line = query_line.split()
		if verbose :
			print('@@ BlastHSP.parse() query_line:', query_line)
		if line and line[0] == self.query_line_start :
			pass
		else :	
			print("Error in parsing Query line of HSP", query_line, file=sys.stderr)
			raise FormatError( "query line:" + query_line )
			return -1

		if len(line) == 4 :
			#save query and hit info temporarily
			query_start_index = int(line[1]) - 1
			query_end_index = int(line[3])
			query = line[2]
			#positive line
			line = fp.readline()
			
			#hit line
			hit_line = fp.readline()
			line = hit_line.split()
			if verbose :
				print('@@ BlastHSP.parse() hit_line:', hit_line)
			if line[0] == self.hit_line_start :
				pass
			else  :
				print("Error in parsing Query line of HSP", query_line, file=sys.stderr)
				raise FormatError( "query line:" + query_line )
				return -1

			if len(line) != 4 :
				raise FormatError( 'hit line:' + hit_line )
				return -1

			hit_start_index = int(line[1]) - 1
			hit_end_index = int(line[3])
			hit = line[2]
			fp.readline() #consume one line
		else :
			print("Error in parsing Query line of HSP", query_line, file=sys.stderr)
			raise FormatError( "query line:" + query_line )
			return -1
			
		while True :
			#query line
			query_line = fp.readline()
			if verbose :
				print('@@ BlastHSP.parse() in while query_line:', query_line)
			if not query_line :
				break
			if not query_line.startswith(self.query_line_start) :
				break
				
			line = query_line.split()
			if len(line) != 4 :
				print("Error in parsing Query line of HSP", file=sys.stderr)
				raise FormatError( "query line:" + query_line )
				return -1

			query_end_index = int(line[3])
			query += line[2]
			
			#positive line
			line = fp.readline()
			
			#hit line
			hit_line = fp.readline()
			if verbose :
				print('@@ BlastHSP.parse() in while hit_line:', hit_line)

			if not hit_line :
				raise FormatError( "hit line:" + hit_line )
				return -1
			if not hit_line.startswith(self.hit_line_start) :
				raise FormatError( 'hit line:' + hit_line )
				return -1

			line = hit_line.split()
			if len(line) != 4 :
				raise FormatError( 'hit line:' + hit_line )
				return -1
			hit_end_index = int(line[3])
			hit += line[2]
			fp.readline() #consume one line

		self.query = query
		self.hit = hit
		self.query_start_index = query_start_index 
		self.hit_start_index = hit_start_index + self.index_offset
		self.query_end_index = query_end_index
		self.hit_end_index = hit_end_index + self.index_offset
		
		return 1

	
	
	def get_aligned_hit_sequence( self, query_length ) :
		'''
		returns hit aligned to query without gap.
		'''
		aln_fp = StringIO()

                #Pre alingment insertion
                pre_aln_string = '-'* self.query_start_index
                aln_fp.write( pre_aln_string )

                #alignment main part
                for a, b in zip( self.query, self.hit ) :
                        if a.isalpha() :
                                aln_fp.write(b)

                #post alignment insertion
                post_aln_length = query_length - aln_fp.tell()
                post_aln_string = '-'*post_aln_length
                aln_fp.write( post_aln_string )
                aln_fp.write('\n')
         
                return aln_fp.getvalue()

	def get_aligned_hit_sequence_with_query_gap( self, query_length, max_gap_info ) :
		'''
		returns hit aligned to query with gap.

		Note that the max_gap_info needs to be calculated 
		to know how many insertions should be made compared to the
		hit with longest insertions at each position.
		'''

		aln_fp = StringIO() 

		pre_aln_gap_length = self.query_start_index
		for i in range( self.query_start_index ) :
			if i in max_gap_info :
				pre_aln_gap_length += max_gap_info[i]
			
		pre_aln_string = '-'*pre_aln_gap_length
		aln_fp.write( pre_aln_string )

		#alignment main part
		query_index = self.query_start_index

		for a, b in zip( self.query, self.hit ) :
			if a.isalpha() :
				if query_index in max_gap_info :
					max_insertion = max_gap_info[query_index]
				else :
					max_insertion = 0

				#adjust gap
				if max_insertion and max_insertion > len(gapped_hit) :
					center = len( gapped_hit )/2
					diff = max_insertion-len(gapped_hit)
					gapped_hit = gapped_hit[:center] + '-'*diff + gapped_hit[center:]
				aln_fp.write( gapped_hit )
				aln_fp.write( b )

				#reset gapped_hit
				gapped_hit = ''
					
				query_index += 1
			else :
				gapped_hit += b

	
		post_aln_gap_length = query_length - self.query_end_index
		for i in range( self.query_end_index+1, query_length ) :
			if i in max_gap_info :
				post_aln_gap_length += max_gap_info[i]

		aln_fp.write( '-'*post_aln_gap_length )
		aln_fp.write( '\n' )
			
		
	def get_gapped_positions( self ) :
		'''
		returns a dictionary of query positions associated with non-zero gap lengths.
		Note that the length is gap before the position indices.

		E.g. the query and gap 

		---AAA--AA-A--
		AAAAA-AAAAAAAA

		will return dictionary of 
		{ 0:3, 3:2, 5:1, 6:2 }.

		This information is desgined to be used to reconstruct the
		gapped alignment based on pairwise alignment result of PSIBLAST.
		'''
		
		gap_info = {}
		query_index = self.query_start_index
		gap_count = 0
		for i in range( len( self.query ) ) :
			a = self.query[i]
			#b = self.hit[i]

			if a.isalpha() :
				#save the current gap count if
				#it is 1 or more.
				if gap_count > 0 :
					gap_info[query_index] = gap_count

				query_index += 1 #move to the next query position
				gap_count = 0 #reset the gap count!
			#gaps in the position
			#adds gap count!
			else  :
				gap_count += 1 
		else :
			#after the query...
			if gap_count > 0 :
				gap_info[query_index] = gap_count

		return gap_info
	
	def check_within( self, pos, query=True ) :
                '''
                Checks if the given position is within Query (or Hit) range.
                check_within( pos ) returns -1 if pos is before the range.
                0 if pos is within the range.
                1 if pos is after the range.
                '''
                start = self.query_start_index
                end = self.query_end_index

                if not query : 
                        start = self.hit_start_index
                        end = self.hit_end_index
         
                if pos <= start : 
                        return -1
                elif start < pos < end : 
                        return 0
                else : #elif end <= pos : same as this
                        return 1


	def is_overlapping( self, alignment=None, s2=None, e2=None, query=True ) :
		'''
		Check fore the query (if query==True) or hit side overlapping 
		between self and alignment or hsp.
		s2 and e2 can be specified for convenient checking.
		'''
		if query :
			s1, e1 = self.query_start_index, self.query_end_index
		else :
			s1, e1 = self.hit_start_index, self.hit_end_index

		if alignment != None:
			if query :
				s2, e2 = alignment.query_start_index, alignment.query_end_index
			else :
				s2, e2 = alignment.hit_start_index, alignment.hit_end_index

		if s2 == None or e2 == None :
			raise TypeError( "alignment or (s2,e2) should be specified!" )

			
		if s1 <= s2 < e1 or s2 <= s1 < e2 or s1 < e2 <= e1 or s2 < e1 <= e2 :
			return True
		else :
			return False



	def remove_overlapped_region( self, alignment, query=True ) :
		'''
		Removes the overlapping region between the input alignment and self.
		If query == True, the query indices are used to determine the overlap,
		otherwise, it uses the hit indices.

		Note that this function does not change self.
		But, it returns a list of modified (if overlapped region is removed)
		alignment objects.
		If no modification made, it will return Null list.
		'''
		if verbose :
			print("@@ remove_overlapped_region @@")
			print(str(self))
			print(str(alignment))

		if query :
			s2, e2 = alignment.query_start_index, alignment.query_end_index
		else :
			s2, e2 = alignment.hit_start_index, alignment.hit_end_index

		if verbose :
			if query :
				print("removing overlaps in queries...")
				s1, e1 = self.query_start_index, self.query_end_index
			else :
				print("removing overlaps in hits...")
				s1, e1 = self.hit_start_index, self.hit_end_index
			print('between (%s-%s) and (%s-%s).'%(s1, e1, s2, e2))

		#divide at the starting point
		pre_alignment, middle_alignment = self.divide( s2, query )
		if pre_alignment == None and middle_alignment == None :
			middle_alignment = self

		#dividing at the ending point
		if middle_alignment.check_within( e2 ) == 0 :
			#this check is done for the boundary condition.
			middle_alignment, post_alignment = middle_alignment.divide( e2, query )
		else :
			post_alignment = None

		return_list = []
		if pre_alignment :
			return_list.append( pre_alignment )
		if post_alignment :
			return_list.append( post_alignment )

		return return_list
		

		
	def divide( self, breaking_position, query=True ) :
		'''
		Divides the current alignment object and returns two
		new alignment objects at the given breaking position.
		Note that the breaking position is 0 a based index of the query sequence.

		If the breaking position is outside of the alignment boundary,
		returns None, None

		None, None will be returned for the matching number of variables. :)
		'''

		if self.check_within( breaking_position, query=query ) == 0 :
			pre_alignment = copy.deepcopy( self )
			try :
				if query : 
					pre_alignment.adjust_query_region( self.query_start_index, breaking_position )
				else : 
					pre_alignment.adjust_hit_region( self.hit_start_index, breaking_position )
			except ValueError:
				pre_alignment = None

			post_alignment =copy.deepcopy( self )
			try :
				if query :
					post_alignment.adjust_query_region( breaking_position, self.query_end_index )
				else  :
					post_alignment.adjust_hit_region( breaking_position, self.hit_end_index )
			except ValueError :
				post_alignment = None

			return pre_alignment, post_alignment

		return None, None



	def divide_multiple_points( self, breaking_positions ) :
		'''
		Divides the current alignment objects at the given list of points
		and returns a list of alignment objects..

		If no breaking positions were able to break the alignment,
		a null list will be returned.
		'''

		breaking_positions = list( breaking_positions )
		breaking_positions.sort() #needs to be sorted list!!

		divided_alignments = []
		target = self
		for breaking_position in breaking_positions :
			if not target :
				break

			pre, post = target.divide( breaking_position )
			if pre == None :
				target = post
				continue
			elif post == None :
				divided_alignments.append( pre )
				break
			else :
				divided_alignments.append( pre )
				target = post
		else :
			#picking up the last object!
			if target and target != self :
				divided_alignments.append( target )

		return divided_alignments


        def find_corresponding_query_start_index_with_real_index( self, hit_index ) :
                '''
                Find query start point.
                Since query string may have a gap for the exact point equivalent for hit_index
                the function will return right after the non-gap query index and the corresponding hit_index
                '''
                qi = self.query_start_index -1
                hi = self.hit_start_index -1
                qi_last = qi
                hi_last = hi
		real_index=0

                for i,(a,b) in enumerate(zip(self.query, self.hit)) :
                        querygap = 1
                        if a.isalpha()  :
                                qi += 1
                                querygap = 0

                        if b.isalpha() :
                                hi += 1

                        if hi >= hit_index and querygap == 0 :
                                qi_last = qi
                                hi_last = hi
				real_index = i
                                break
		else :
			#cannot find...
			qi_last += 1
			hi_last += 1
		#else :
			#raise ValueError( "hit_index %s was not found within the hit range %s-%s." %
				#(hi, self.hit_start_index, self.hit_end_index ) )

                return qi_last, hi_last, real_index



        def find_corresponding_hit_start_index_with_real_index( self, query_index ) :
                '''
                Find query start point.
                Since query string may have a gap for the exact point equivalent for hit_index
                the function will return right after the non-gap query index and the corresponding hit_index
                '''
                qi = self.query_start_index -1
                hi = self.hit_start_index -1
                qi_last = qi
                hi_last = hi
		real_index=0
                for i,(a,b) in enumerate(zip(self.query, self.hit)) :
                        hitgap = 1
                        if a.isalpha()  :
                                qi += 1 

                        if b.isalpha() :
                                hitgap = 0
                                hi += 1 

                        if qi >= query_index and hitgap == 0 :
                                qi_last = qi
                                hi_last = hi
				real_index = i
                                break   
		else :
			#cannot find
			qi_last += 1
			hi_last += 1
		
		#else :
			#raise ValueError( "query_index %s was not found within the query range %s-%s." %
				#(query_index, self.query_start_index, self.query_end_index ) )

                return qi_last, hi_last, real_index



        def find_corresponding_query_end_index_with_real_index( self, hit_index ) :
                '''
                Find query end point.
                Since hit string may have a gap for the exact point equivalent for query_index
                the function will return right before the non-gap query index and the corresponding hit_index
                '''
                qi = self.query_start_index -1
                hi = self.hit_start_index -1
                qi_last = qi
                hi_last = hi
		real_index = 0
                for i,(a,b) in enumerate(zip(self.query, self.hit)) :
                        querygap = 1
                        if a.isalpha()  :
                                qi += 1
				querygap=0

                        if b.isalpha() :
                                hi += 1

                        if hi <= hit_index and querygap == 0 :
                                qi_last = qi
                                hi_last = hi
				real_index = i

                        if hi >= hit_index :
                                break
                else :
                        return qi_last+1, hi_last+1, real_index+1

                return qi_last, hi_last, real_index
                #return qi_last+1, hi_last+1, real_index+1



        def find_corresponding_hit_end_index_with_real_index( self, query_index ) :
                '''
                Find hit end point.
                Since hit string may have a gap for the exact point equivalent for query_index
                the function will return right before the non-gap query index and the corresponding hit_index
                '''
                qi = self.query_start_index -1
                hi = self.hit_start_index -1
                qi_last = qi
                hi_last = hi
		real_index = 0
                for i,(a,b) in enumerate(zip(self.query, self.hit)) :
                        hitgap = 1
                        if a.isalpha()  :
                                qi += 1

                        if b.isalpha() :
                                hi += 1
                                hitgap = 0

                        if qi <= query_index and hitgap == 0 :
                                qi_last = qi
                                hi_last = hi
				real_index = i

                        if qi >= query_index :
                                break
                else :
                        return qi_last+1, hi_last+1, real_index+1

                return qi_last, hi_last, real_index
                #return qi_last+1, hi_last+1, real_index+1

		
	
	def adjust_query_region( self, new_query_start, new_query_end ) :
		'''
		Change current query and hit region by reducing to new query start and qeury end.
		new_query_start and new_query_end are 0 based indices for the query sequence.
		'''
		if verbose :
			print('adjust query region from %s-%s to %s-%s.'%(self.query_start_index, 
				self.query_end_index, new_query_start, new_query_end ))
			

		#check if the change can be done!
		if not self.is_overlapping( s2=new_query_start, e2=new_query_end ) :
			raise ValueError( "Query adjustment has failed due to invalid change from %s-%s to %s-%s." % 
				(self.query_start_index, self.query_end_index, new_query_start, new_query_end) )

		
		self.modified = 1
		new_qsi, new_hsi, si1 = self.find_corresponding_hit_start_index_with_real_index( new_query_start )
		new_qei, new_hei, ei1 = self.find_corresponding_hit_end_index_with_real_index( new_query_end )
		

		new_query = self.query[ si1:ei1 ]
		if not new_query :
			raise ValueError( "Query adjustment has resulted Null query sequence changing from %s-%s to %s-%s." % 
				(self.query_start_index, self.query_end_index, new_query_start, new_query_end) )

		new_query = self.query[si1:ei1]
		new_hit = self.hit[si1:ei1]

		if verbose :
			print('Old alignment:')
			print(str(self))

		self.query = new_query
		self.hit   = new_hit
		self.query_start_index, self.query_end_index = new_qsi, new_qei
		self.hit_start_index, self.hit_end_index = new_hsi, new_hei

		if verbose :
			print('New alignment:')
			print(str(self))

		return 1


	def adjust_hit_region( self, new_hit_start, new_hit_end ) :
		'''
		Change current query and hit region by reducing to new hit start and new hit end.
		new_hit_start and new_hit_end are 0 based indices for the hit sequence.
		'''
		if verbose :
			print('adjust hit region from %s-%s to %s-%s.'%(self.hit_start_index, 
				self.hit_end_index, new_hit_start, new_hit_end ))

		#check if the change can be done!
		if not self.is_overlapping( query=False, s2=new_hit_start, e2=new_hit_end ) :
			raise ValueError( "Hit adjustment has failed due to invalid change from %s-%s to %s-%s." % 
				(self.hit_start_index, self.hit_end_index, new_hit_start, new_hit_end) )

		
		self.modified = 1
		new_qsi, new_hsi, si1 = self.find_corresponding_query_start_index_with_real_index( new_hit_start )
		new_qei, new_hei, ei1 = self.find_corresponding_query_end_index_with_real_index( new_hit_end )
		

		new_query = self.query[ si1:ei1 ]
		new_hit = self.hit[si1:ei1]
		if not new_query :
			sys.stdout.flush()
			raise ValueError( "Query adjustment has resulted Null query sequence changing from %s-%s to %s-%s (%s-%s)." % 
				(self.query_start_index, self.query_end_index, new_qsi, new_qei, si1, ei1) )

		if verbose :
			print('Old alignment:')
			print(str(self))

		self.query = new_query
		self.hit   = new_hit
		self.query_start_index, self.query_end_index = new_qsi, new_qei
		self.hit_start_index, self.hit_end_index = new_hsi, new_hei

		if verbose :
			print('New alignment:')
			print(str(self))

		return 1

        def convert_position_to_real_index( self, pos, query=True ) : 
                '''
                Convert sequence position (for hit or query) into real index according to aligned_seq.                The aligned_seq should be hit string or query string.
                if conversion failes then this function returns length of aligned_seq
                '''
		if query :
			aligned_seq = self.query
			start_index = self.query_start_index
		else :
			aligned_seq = self.hit
			start_index = self.hit_start_index
			
                current_pos = start_index -1
                real_index = -1
                for i,a in enumerate(aligned_seq) :
                        if a.isalpha() :
                                current_pos += 1
         
                        if pos == current_pos :
                                real_index = i
                                break   
                else :
                        real_index = len( aligned_seq )
         
                return real_index

