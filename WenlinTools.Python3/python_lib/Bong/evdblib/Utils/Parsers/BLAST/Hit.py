from evdblib.Utils.Parsers.BLAST import *
from .Alignment import Alignment
from io import StringIO
from evdblib.Utils.Parsers.FASTA import parse_a3m, FASTA

verbose = 0

class Hit :
	"""
	This class can contains multiple HSPs for a hit.
	GI is used as a unique identifier for a hit.
	"""
	def __init__(self) :
		self.gi = 0
		self.protein_name = ""
		self.hsps = []
		self.length=0
		
		self.hit_start = hit_record_start #'>gi'
		self.hit_info_end = '          Length ='
		self.index_offset = 0 #for gir correction
		self.hsp_length_cutoff = 30 #length cutoff for hsp

	def has_data( self ) :
		'''
		returns True when Hit is successfully parsed and has some
		alignments in it.

		Otherwise, returns False.
		'''
		if self.hsps :
			return True
		else :
			return False

	def add_hit_as_hsp( self, hit ) :
		self.hsps.extend( hit.hsps )

	def write_hit_fasta( self, fp ) :
		for hsp in self :
			hsp.write_hit_fasta( fp, self.gi )

	def write_a3m_fasta( self, fp, queryseq, remove_gap=False, combine_hsps=False ) :
		'''
		Write A3M fasta strings into fp.
		Note that the remove_gap sets to be True automatically
		when combine_hsps is True.
		This limitation is due to save some computation rather than
		biologically motivated.
		'''

		if combine_hsps :
			remove_gap = True

		if verbose :
			print("#################writing a3m fasta##############")
			print(self.protein_name)

		#save fasta records into a temporary 
		#output buffer to combine
		outstream = StringIO()
		for hsp in self.hsps :
			hsp.write_a3m_fasta( outstream, self.protein_name, queryseq, remove_gap=remove_gap )

		if verbose :
			print("written fastas:")
			print(outstream.getvalue())

		#rewind
		outstream.seek(0)

		if combine_hsps and len(self.hsps) > 1 :
			#read fasta records in self.hsps
			fastas = parse_a3m( outstream )
			#prepare to reuse outstream
			outstream.seek(0)
			outstream.truncate()
			#combining
			#assuming remove_gap is True.
			combined_fasta = FASTA()
			combined_fasta.header = fastas[0].header
			for i in range(len(fastas[0].sequence)) :
				for j in range(len(self.hsps)) :
					if fastas[j].sequence[i].isalpha() :
						outstream.write(fastas[j].sequence[i])
						break
				else :
					outstream.write('-')
			combined_fasta.sequence = outsteam.getvalue()
		else :
			fp.write( outstream.getvalue() )
			
	def __len__(self) :
		return self.hsps.__len__()

	def __iter__( self ) :
		return self.hsps.__iter__()

	def __getitem__(self, i ) :
		return self.hsps[i]

	def __str__( self ) :
		s = 'gi: %s\nprotein: %s\n'%(self.gi, self.protein_name)
		return s+'\n'.join( [ str(hsp) for hsp in self ] )

	def get_aligned_hit_sequence_with_query_gaps( self, query_length, max_gap_info ) :
		'''
		returns hit sequence aligned to the query with query gap retained.
		The return sequence is similar to the blastpgp -m6 option.

		max_gap_info is a dictionary with the gap positions as keys and
		maximum length of query insertions as values.
		'''

		annotation = self.gi[:20]
		if not annotation :
			annotation = os.path.split( os.tempname )[-1]

		return ''.join( [ '%-20s %s'%(annotation, hsp.get_aligned_hit_sequence_with_query_gaps(query_length, max_gap_info )) for hsp in self.hsps ] )
		

	def get_aligned_hit_sequence( self, query_length, combine_hsps=False ) :
		'''
		returns hit sequence aligned to the query.
	
		if combine_hsps is True then the first non-gap amino acid
		will be returned for each position.
		
		'''
                sequences = []
                for hsp in self.hsps :
                        sequences.append( hsp.get_aligned_hit_sequence( query_length ) )

		if combine_hsps == False :
			annotation = self.gi[:20]
			if not annotation :
				annotation = os.path.split( os.tempname )[-1]

			return ''.join( [ '%-20s %s'%(annotation, seq) for seq in sequences ] )

                if len(sequences ) == 0 :
                        return ""
                elif len(sequences) == 1 :
                        return "%-20s %s" %( self.gi[:20], sequences[0] )

                str_fp = cStringIO.StringIO()
                for j in range( len( sequences[0] )-1 ) : #position j and it is adjust to get rid of the last new line
                        for i in sequences : #sequence i
                                if i[j].isalpha() :
                                        str_fp.write( i[j] )
                                        break
                        else :  
                                str_fp.write( '-' )
                str_fp.write( '\n' )

                return "%-20s %s" %( self.gi[:20], str_fp.getvalue() )


	def parse_a3m_fasta( self, hitfasta, queryfasta ) :
		hsp = Alignment( offset=self.index_offset )
		hsp.parse_a3m_fasta( hitfasta, queryfasta )
		self.hsps.append( hsp )

		self.gi = hitfasta.gi
		self.protein_name = hitfasta.get_header()
		self.length = hitfasta.parse_hit_length()

		return 0


	def parse( self, fp ) :
		'''
		Returns 0 for nomal end 
		Returns -2 for nomal psiblast/blast end 
		Returns -1 for errorneous parsing
		returns 1 for successful parsing
		'''
		hit_line = fp.readline()
		if verbose :
			print('first hit line check:', hit_line)

		while not hit_line.startswith(self.hit_start) :
			if hit_line.startswith( psiblast_end_of_result_marker ) : 
				if verbose :
					print("End of Blast Result reached!", hit_line) 
				fp.seek(-len(hit_line), 1 ) #rewind one line back for the parameters for blast run
				raise NoMoreIteration( "Alignment parsing reached at the end of BLAST result!" )
				return -2

			elif hit_line.startswith(psiblast_start) :
				if verbose :
					print("End of Blast Iteration reached!", hit_line) 
				#fp.seek(-len(hit_line), 1 ) #rewind one line back for the parameters for blast run
				raise NoMoreHit( "Alignment parsing reached at the end of Iteration")
				return -2

			hit_line = fp.readline()
			if verbose :
				print('subsequent hit line check:', hit_line)

			if not hit_line  :
				raise NoMoreIteration( 'Alignment parsing reached at the end of BLAST result!' )
				return -1

		if verbose :
			print('finally supposed to be correct hit line:', hit_line)

		#getting gi from a line
		#>gi|74188315|dbj|BAE25816.1| unnamed protein product [Mus musculus]
		try :
			gi = hit_line.split('|')[1] #representative hit gi is the first one
		except :
			gi = hit_line.split()[0]

		protein_name = hit_line[1:-1]
		#consume all other lines
		#stops at 
		#          Length = 566
		info_line = fp.readline()
		while info_line and not info_line.startswith(self.hit_info_end) :
			info_line = fp.readline()
		
		length = int(info_line.split()[-1])
		fp.readline() #consume one empty line

		#########################################
		#adjustment for gir: gi with range
		#gir consists of gi_start..end
		#e.g. 74188315_1..120
		index_offset = 0
		if gi.find('_') == -1 :
			pass
		else :
			gi, r = gi.split('_')
			try :
				index_offset = int(r.split('..')[0])
			except :
				gi = '_'.join( [gi,r] )
				pass
		#########################################

		self.gi = gi
		self.protein_name = protein_name
		self.length = length
		self.index_offset = index_offset

		hit_fp = cStringIO.StringIO()
		while True :
			line = fp.readline()
			if not line :
				break

			if line.startswith(self.hit_start) :
				fp.seek( -len(line), 1 )
				break
			hit_fp.write(line)

		hit_fp.seek(0) #rewinde hit_fp

		while 1 :
			try :
				hsp = Alignment( offset=self.index_offset )
				parse_flag = hsp.parse(hit_fp)
				if verbose :
					print('@ after the subsequent BlastHSP.parse(), parse_flag:', parse_flag)
			except NoMoreAlignment :
				if verbose :
					print('@ End of BlastHit.parse(), parse_flag:', parse_flag)
				return parse_flag
			else :
				self.hsps.append( hsp )

		if verbose :
			print('@ End of BlastHit.parse(), but should not be visited!:', parse_flag)
			
		'''
		if verbose :
			print '@ after the first BlastHSP.parse(), parse_flag:', parse_flag
			
		while parse_flag == 1 :
			self.hsps.append(hsp)
			hsp = BlastHSP( offset=self.index_offset )
			parse_flag = hsp.parse( hit_fp )
			if verbose :
				print '@ after the subsequent BlastHSP.parse(), parse_flag:', parse_flag
		if verbose :
			print '@ End of BlastHit.parse(), parse_flag:', parse_flag
	
		if parse_flag == 0 :
			return 1 #just HSP end!!
		'''

		
	
	def purge_overlapping_hsps( self, indices, pssm ) :
		'''
		Purges HSPs by selecting good regions in the HSPs
		'''

		if not indices :
			return

		if verbose :
			print('###', 'Purging Hit', self.gi, '###')
		
                for i, hsp in enumerate(self.hsps[:]) :
                        #select correct region among the segments
                        #first, divide query string into fragments using given indices
                        #second, calculate the score for each fragments.
                        new_hsps = hsp.divide_multiple_points( indices )
                        if verbose :
                                print('HSP #%d are divided into'%i, new_hsps)

                        if len(new_hsps)<=1 :
                                if verbose :
                                        print('purging cannot be done', 'new_hsps=', new_hsps)
                                continue

                        if verbose :
                                for j,h in enumerate( new_hsps ) :
                                        print('new_hsps:', j, h.query, h.hit)

                        pssm_scores = [pssm.cal_pssm_score(h) for h in new_hsps]
                        scores_hsps = list(zip(pssm_scores, new_hsps))
                        scores_hsps.sort()
                        if verbose :
                                print(scores_hsps)
                        self.hsps[i] = scores_hsps[-1][1]

                        if verbose :
                                print('score-hsp pairs:', scores_hsps)
                                print('##selected HSP:', self.hsps[i].query_start_index, self.hsps[i].query, self.hsps[i].query_end_index, '\n')

                #check the correct fragment mapping
                #by overlapping of regions

                for i, hsp in enumerate( self.hsps ) :
                        for hsp2 in self.hsps[i+1:] :
                                if hsp.is_overlapping( hsp2 ) :
                                        if verbose :
                                                print("Warning HSP selection has problem!")
                                        #removing all hsps
                                        self.hsps = []
                                        ########################
                                        #this hit instance should also be removed
                                        ########################
                else :  
                        if verbose :
                                print('HSP checking is done!')


                """
                Blocked out the removing small fragment part
                #remove small fragments
                #maybe not necessary
                for hsp in self.hsps[:] :
                        if len(hsps.hit) <= self.hsp_length_cutoff :
                                self.hsps.remove( hsp )
                """

	def remove_overlapped_hit_region( self, hit ) :
                """
                Remove all the common hit matching regions.
                """

                for hsp2 in hit.hsps :
                        for i, hsp in enumerate( self.hsps[:] ) : #to use slice to avoid problem in removing elements.
				if verbose :
                                        print("#### Hit", self.gi, "####")
                                        print("hit1:", hsp.hit_start_index, hsp.hit_end_index)
                                        print("hit2:", hsp2.hit_start_index, hsp2.hit_end_index)
	
				if not hsp.is_overlapping( hsp2, query=False ) :
					if verbose :
						print("\tskipping:\thit1 is not overlapping with hit2.", file=sys.stderr)
					continue
 
				#remove overlapping hit region.
				#query=False should be set!!
                                new_hsps = hsp.remove_overlapped_region( hsp2, query=False )

			
				if new_hsps :
					#if new_hsps has members,
					#the members are good without the overlapping region.
					#and the original hsp should be deleted.
                                        self.hsps.remove( hsp )
					self.hsps.extend( new_hsps )
				else :
					#if new_hsps do not have members,
					#it means the hsp do not overlap with hsp2 or
					#contained within the range.
					#Note that the overlapping check should be done
					#on the hit side.. needs the option (query=False).
					if hsp.is_overlapping( hsp2, query=False ) :
						#need to delete the hsp.
						self.hsps.remove( hsp )
					else :
						#do not need to do anything.
						pass


        def remove_overlapping_hsps( self ) :
                ''' 
                This is a testing function to solve the profile corruption.
                The idea is removing all possibly contaminated HSPs from the hit list.
                '''
                hsps_toberemoved = set()
                for hsp1 in self.hsps :
                        for hsp2 in self.hsps :
                                if hsp1 == hsp2 :
                                        continue
                                                
                                if hsp1.is_overlapping( hsp2 ) :
                                        hsps_toberemoved.add( hsp1 )
                                        break   

                for hsp in hsps_toberemoved :
                        self.hsps.remove( hsp )

