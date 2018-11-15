import os, sys, io, traceback, re
from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentMethodRecord
verbose = 0

class TMalignParseError( Exception ) :
	pass

class TMalign :
	def __init__( self, alignment_string, cafn1, cafn2, id1=None, id2=None, strict=False, seqfn1=None, seqfn2=None, verbose=False ) :
		'''
		Main Parser class for TMalign structure alignment program result.

		Note that the first argument is alignment string not dccp file like in DaliLite.
		'''

		self.alignment_string = alignment_string
		self.cafn1 = cafn1
		self.cafn2 = cafn2
		self.id1 = id1
		self.id2 = id2
		self.seqfn1 = seqfn1
		self.seqfn2 = seqfn2
		self.strict = strict
		self.verbose = verbose

		self.alignment = PairwiseAlignmentMethodRecord(id1=self.id1, id2=self.id2, method_name='tmalign')

		if not self.id1 :
			self.id1 = parse_sequence_filename( self.cafn1)[1]

		if not self.id2 :
			self.id2 = parse_sequence_filename( self.cafn2)[1]
		
		self.fasta1 = None
		self.fasta2 = None
		if self.seqfn1 :
			self.fasta1 = FASTA( self.seqfn1 )
		if self.seqfn2 :
			self.fasta2 = FASTA( self.seqfn2 )

		if self.alignment_string and self.id1 and self.id2 :
			self.alignment = self.parse( io.StringIO(self.alignment_string), self.fasta1, self.fasta2, self.id1, self.id2 )

	def parse( self, fp, fasta1, fasta2, id1, id2 ) :
		'''
		parse a TMalign result and returns a PairwiseAlignmentMethodRecord object.
		'''
		alignment =  PairwiseAlignmentMethodRecord( id1=id1, id2=id2, method_name='tmalign' )
		
		for i in range( 24 ) :
			l = fp.readline()
			#case of no more record
			if not l and i == 0:
				#return None
				if self.strict :
					raise TMalignParseError( "Expected lines are missing.")
				else :
					return alignment
			if not l :
				if self.strict :
					raise TMalignParseError( "Expected lines are missing.")
				else :
					return alignment

			if i == 8 :
				child =  l[13:17]
			elif i == 9 :
				parent = l[13:17]
			elif i == 11 :
				tm_score =  float( l[43:50] )  #tm_align score
				rmsd = float( l[26:32] )  #RMSD
				aligned_length = int( l[15:19] )  # aligned elgnth
				identity = float( l[55:60] )  #identity
			elif i == 20 :
				length = len( l )
				c_alignment = l[:-1]
			elif i == 21 :
				equivalence = l[:-1]
				if len(l) == length :
					pass
				else :
					if self.verbose :
						print('WARNING: aligned length is not correct\n' + 'parsed_len:', length, 'real:', len(l[:-1]) ,'\n'+l[:-1], file=sys.stderr)
					if self.strict :
						raise TMalignParseError( "Aligned length is not correct.", id1, id2, length, len(l[:-1]) )
					else :
						return alignment
						
			elif i==22 :
				p_alignment = l[:-1]
				if len(l) == length :
					pass
				else :
					if self.verbose :
						print('WARNING: aligned length is not correct\n' + 'parsed_len:', length, 'real:', len(l[:-1]) ,'\n'+l[:-1], file=sys.stderr)
					if self.strict :
						raise TMalignParseError( "Aligned length is not correct.", id1, id2, length, len(l[:-1]) )
					else :
						return alignment

		#return (parent, child, tm_score, rmsd, aligned_length, identity, p_alignment, c_alignment )
		alignment.alignment1 = p_alignment
		alignment.alignment2 = c_alignment
		alignment.raw_score = tm_score * aligned_length
		alignment.normalized_score = tm_score

		return alignment




#legacy code for parsing TMalign alignment result.
def parse_tm_record( fp ) :
        for i in range( 18 ) :
                l = fp.readline()
                #case of no more record
                if not l and i == 0:
                        return None
                if not l :
                        return None, None
                if i == 8 :
                        child =  l[13:17]
                elif i == 9 :
                        parent = l[13:17]
                elif i == 11 :
                        tm_score =  float( l[43:50] )  #tm_align score
                        rmsd = float( l[26:32] )  #RMSD
                        aligned_length = int( l[15:19] )  # aligned elgnth
                        identity = float( l[55:60] )  #identity
                elif i == 14 :
                        length = len( l )
                        c_alignment = l[:-1]
                elif i == 15 :
                        equivalence = l[:-1]
                        if len(l) == length :
                                pass
                        else :
                                print('WARNING: aligned length is not correct\n' + 'parsed_len:', length, 'real:', len(l[:-1]) ,'\n'+l[:-1], file=sys.stderr)
                elif i==16 :
                        p_alignment = l[:-1]
                        if len(l) == length :
                                pass
                        else :
                                print('WARNING: aligned length is not correct\n' + 'parsed_len:', length, 'real:', len(l[:-1]) ,'\n'+l[:-1], file=sys.stderr)

        return (parent, child, tm_score, rmsd, aligned_length, identity, p_alignment, c_alignment )



def remove_xmasks_from_aln( seq, aln ) :
        new_aln = list( aln )
        #aalist = 'ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy'
        #notaalist = 'BUXZbuxz' 
        #xpos = [ (a == 'X' or a=='x') for a in seq ]
        #print xpos

        '''#debugging
        print len(seq), seq
        seq_from_aln = []
        for a in aln :
                if a.isalpha() :
                        acount += 1
                        seq_from_aln.append(a)
        print acount, aln
        print acout, ''.join( seq_from_aln )
        acount = -1
        '''

        acount = -1
        for i, a in enumerate(aln) :
                if a.isalpha() :
                        acount += 1
                        if seq[acount] in notaalist :
                                if ( a == 'A' ) :
                                        new_aln[i] = seq[acount].upper()
                                elif (a == 'a' ) :
                                        new_aln[i] = seq[acount].lower()
                                else :
                                        print("Error! this message should not be seen!", seq, aln, file=sys.stderr)
                                        sys.exit()

        return ''.join(new_aln)


