
from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentMethodRecord
from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentRecords, AlignmentAddingError
from evdblib.Utils import parse_profile_filename, rewind

debug = 0

class HHsearchNullResultError( Exception ) :
	pass

class HHsearchParseError( Exception ) :
	pass

class HHsearch :
	def __init__ ( self, queryid, result_fn=None, iteration=None, use_iteration_info=True, strict_id_match=False ) :
		'''
		HHsearch parses the HHsearch result file.
		Its major functionality is parsing the alignments and
		scores and return the GenericPairwiseAlignment object
		So that the alignments and scores can be easily manipulated.
		'''

		self.queryid = queryid
		self.result_fn = result_fn
		self.alignments = None
		self.iteration = iteration
		self.strict_id_match = strict_id_match # to control id match in query and target lines
		
		#if iteration information is not given,
		#yet still use_iteration_info is True,
		#try to parse the profile filename
		if self.iteration == None :
			result_dir, basename, iter, suffix = parse_profile_filename( result_fn )
			self.iteration = iter

		if self.result_fn :
			self.alignments = self.parse()

	def get_alignments( self ) :
		'''
		returns self.alignments
		'''

		if self.alignments :
			return self.alignments

		else :
			self.alignments = self.parse()
			if not self.alignments :
				raise HHsearchNullResultError( self.result_fn )
			

	def parse( self, check_integrity=1 ) :
		'''
		parse the whole result file 
		and returns alignments.
		'''

		fp = open( self.result_fn )
		line = fp.readline()
		if debug :
			print("first line:", line)
		if check_integrity :
			header, queryid = line.split()
			if header != 'Query' or queryid != self.queryid :
				raise HHsearchParseError( self.result_fn, "Query info is not right!" )

		alignments = PairwiseAlignmentRecords()
		while not line.startswith('No ') :
			line = fp.readline()
			if not line :
				raise HHsearchNullResultError( self.result_fn )
		else :
			if line.startswith( 'No ' ) :
				rewind( fp, len(line) )

		lines = []
		for line in fp :
			if line.startswith('No ') or line.startswith( 'Done!' ) :
				if lines :
					if debug :
						print("*"*50)
						print(''.join(lines))
						print("*"*50)

					methodrecord = self.parse_record( lines )
					try :
						alignments.add( methodrecord )
					except AlignmentAddingError :
						print("Record addition has a problem probably due to duplicated records.", file=sys.stderr)

					#End of record!
					if line.startswith( 'Done!' ) :
						break
					lines = []

			lines.append( line )
		else :
			raise HHsearchParseError( self.result_fn, "Result file format is not complete!" )

		return alignments
		
	

	def parse_record( self, content ) :
		'''
		parses a record in HHsearch result.
		'''

		aln1 = ''
		aln2 = ''
		start1 = 0
		start2 = 0

		#print content

		additional_scores = {}
		query_line_marker = 'Q ' + self.queryid
		for i, l in enumerate(content) :
			if l.startswith( 'No ' ) :
				hit_number = int(l.split()[-1])
				additional_scores['No'] = hit_number
				continue
		
			elif l.startswith('Probab=' ) :
				l = l.split()
				probability = float(l[0][7:])/100
				evalue = float( l[1][8:] )
				content = content[i+1:]
				
				for pair in l :
					scorename, score = pair.split('=')
					if scorename == 'Aligned_columns' :
						additional_scores[scorename] = int(score)
					elif scorename == 'Identities' :
						additional_scores[scorename] = int(score[:-1])
					else :
						additional_scores[scorename] = float(score)
				continue

			elif l == '\n' :
				continue

			elif l.startswith( '>' ) :
				hitid = l.split()[0][1:]

				if not hitid :
					raise HHsearchParseError( "Hit definition line cannot be understood!", hitid )
				hit_line_marker = 'T ' + hitid

			elif l.startswith( query_line_marker ) :
				l = l.split()
				if not aln1 :
					start1 = int(l[2]) -1 #correction for 0 based numbering
				aln1 += l[3]

			elif (not self.strict_id_match) and l.startswith( query_line_marker[:14] ) :
				l = l.split()
				if not aln1 :
					start1 = int(l[2]) -1 #correction for 0 based numbering
				aln1 += l[3]

			elif l.startswith( hit_line_marker ) :
				l = l.split()
				if not aln2 :
					start2 = int(l[2]) - 1 #correction for 0 based numbering
				aln2 += l[3]

			elif (not self.strict_id_match) and l.startswith( hit_line_marker[:14] ) :
				l = l.split()
				if not aln2 :
					start2 = int(l[2]) - 1 #correction for 0 based numbering
				aln2 += l[3]
			
		methodname = 'hhsearch'
		if self.iteration :
			methodname = methodname+'_'+str(self.iteration)

		if debug :
			print('##parse result##')
			print('queryid:', self.queryid)
			print('hitid:', hitid)
			print('al1:', aln1)
			print('al2:', aln2)
			print("additional scores:", additional_scores)

		methodrec = PairwiseAlignmentMethodRecord( 
			id1=self.queryid, id2=hitid,
			method_name = methodname, 
			alignment1 = aln1, alignment2 = aln2,
			start_position1 = start1, start_position2 = start2,
			additional_scores = additional_scores,
			raw_score = additional_scores['E-value'],
			norm_score = additional_scores['Probab'] / 100.0
		)
		if debug :
			print("final methodrec:")
			print(methodrec)

		#return ( evalue, probability, start1, start2, aln1, aln2 )
		return methodrec

######################################################
#Codes in the following part are legacy tools.
#They should be used only for reference purposes.
######################################################
# HHsearch parsing code
def parse_hhsearch( content ) :
        for i, l in enumerate(content) :
                if l[:7] == 'Probab=' :
                        l = l.split()
                        probability = float(l[0][7:])/100
                        evalue = float( l[1][8:] )
                        raw_score = float( l[2][6:] )
                        content = content[i+1:]
                        break

        length = len(content)
        if length%11 == 2 :
                pass
        else :
                print("Error!! The format is not right!! code1", file=sys.stderr)
                sys.stderr.writelines( content )
                sys.exit()

        nlines = length/11
        aln1 = ''
        aln2 = ''
        start1 = 0
        start2 = 0
        for i in range( nlines ) :
                #process query line
                l = content[i*11+3]
                l = l.split()
                if l[0] == 'Q' :
                        pass
                else :
                        print("Error!! The format is not right!! code2", file=sys.stderr)
                        sys.stderr.writelines( content )
                        sys.exit()

                if not aln1 :
                        start1 = int(l[2]) -1 #correction for 0 based numbering
                aln1 += l[3]
                #process subject line
                l = content[i*11+3+4]
                l = l.split()
                if l[0] == 'T' :
                        pass
                else :
                        print("Error!! The format is not right!! code3", file=sys.stderr)
                        sys.stderr.writelines( content )
                        sys.exit()

                if not aln2 :
                        start2 = int(l[2]) - 1 #correction for 0 based numbering
                aln2 += l[3]

        ######################
        #modification 08/04/09
        #to gather raw scores too.
        # return ( evalue, probability, start1, start2, aln1, aln2 )
        return ( evalue, probability, start1, start2, aln1, aln2, raw_score )


def run_hhsearch( id1_fn, id2_fn ) :
        hhsearch_cmd = "/usr1/HorAServer/local/hhsearch1.5_2/hhsearch -i %s -d %s -o stdout -b 1 -B 1 -z 1 -Z 1"
        pp = os.popen( hhsearch_cmd % (id1_fn, id2_fn ) )
        stdout = pp.readlines()
        #sys.stdout.writelines( stdout )
        #content = open( 'temp.hhr' ).readlines()
        return parse_hhsearch( stdout )

def run_hhsearch_old( id1_fn, id2_fn ) :
        #id2_fn_only = id2_fn.split('/')[-1]
        id2_fn_only = 'hhsearch_result' #id2_fn.split('/')[-1]
        hhsearch_cmd = "/usr1/HorAServer/local/hhsearch1.5_2/hhsearch -i %s -d %s -o %s.hhr -b 1 -B 1 -z 1 -Z 1 > /dev/null"
        os.system( hhsearch_cmd % (id1_fn, id2_fn, id2_fn_only ) )

        content = open( id2_fn_only + ".hhr" ).readlines()##read output from file!!

        return parse_hhsearch( content )


