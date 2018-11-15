import os, sys, io, traceback, re
from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentMethodRecord
verbose = 0

class FASTParseError( Exception ) :
	pass

class FAST :
	def __init__( self, alignment_string, cafn1, cafn2, id1=None, id2=None, strict=False, seqfn1=None, seqfn2=None ) :
		'''
		Main Parser class for FAST structure alignment program result.

		Note that the first argument is alignment string not dccp file like in DaliLite.
		'''

		self.head_pattern = re.compile( 'FAST ALIGNMENT: (.{4}).\S+\s(.{4}).\S+' )
		self.score_pattern = re.compile( 'L=(\S+)\sSX=(\S+)\sSN=(\S+)\sL1=\S+\sL2=\S+\sRMSD=(\S+)' )
		self.alignment_pattern = re.compile( ' [12]:\t(\S*)\**' )

		self.alignment_string = alignment_string
		self.cafn1 = cafn1
		self.cafn2 = cafn2
		self.id1 = id1
		self.id2 = id2
		self.seqfn1 = seqfn1
		self.seqfn2 = seqfn2
		self.strict = strict
			
		self.alignment = PairwiseAlignmentMethodRecord(id1=self.id1, id2=self.id2, method_name='fast')

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
		parse a FAST result and returns a PairwiseAlignmentMethodRecord object.
		'''
		head_pattern = self.head_pattern
		score_pattern = self.score_pattern
		alignment_pattern = self.alignment_pattern
		alignment =  PairwiseAlignmentMethodRecord( id1=id1, id2=id2, method_name='fast' )
		strict = self.strict
		
		#first line
		l = fp.readline()
		#case of no more record
		if not l:
			if verbose :
				print("No results are found!", file=sys.stderr)
			return alignment

		#print l
		match = head_pattern.match( l )
		#checking for starting point
		if not match :
			if verbose :
				print('WARNING: wrong starting point!!', id1, id2, l, file=sys.stderr)

			while( not match ) :
				l = fp.readline()
				if not l :
					#print >>sys.stderr, 'End of file reached!'
					if self.strict :
						raise FASTParseError( "No results found", id1, id2 ) 
					else :
						return alignment

				if verbose :
					print('following lines:', l, id1, id2, file=sys.stderr)

				match = head_pattern.match( l )

		parent = match.group(1) #id1
		child = match.group(2) #id2

		if verbose :
			if parent != id1 :
				print("WARNING: ID1 is different from expected!", parent, id1, file=sys.stderr)
			if child != id2 :
				print("WARNING: ID2 is different from expected!", child, id2, file=sys.stderr)

		#second line
		l = fp.readline()
		match = score_pattern.match( l )
		aligned_length = int(match.group(1))
		fast_score = match.group(2)
		norm_score = match.group(3)
		rmsd = match.group(4)

		if not aligned_length :
			l1 = fp.readline()
			l2 = fp.readline()
			if l1 == l2 == '\n' :
				return alignment
			else:
				if verbose :
					print("Error! while parsing fast record. Aborted!", parent, child, file=sys.stderr)
				if strict :
					raise FASTParseError( "Incorrect format", id1, id2 ) 
				else :
					return alignment

		l = fp.readline()
		if l != '\n':
			if verbose :
				print("Error, unexpected pattern at line3:", parent, child, l, file=sys.stderr)
			if strict :
				raise FASTParseError( "Unexpected format at line3", id1, id2, l ) 
			else :
				return alignment
				

		#alignment process
		paseq = ''
		caseq = ''
		while ( True ) :
			l = fp.readline()
			if l == '\n' :
				break
			match = alignment_pattern.match( l )
			if not match :
				break
			paseq += match.group(1)

			l = fp.readline()
			match = alignment_pattern.match( l )
			caseq += match.group(1)

			l = fp.readline()
			if l == '\n' :
				pass
			else :
				if verbsoe :
					print("Error, unexpected pattern in alignment region!", parent, child, paseq, caseq, file=sys.stderr)
				if strict :
					raise FASTParseError( "Unexpected format in alignment region", id1, id2, paseq, caseq ) 
				else :
					return alignment


		#need to consume the trailing equvalence parts
		while ( True ) :
			l = fp.readline()
			if l == '\n' :
				l = fp.readline()
				if l == '\n' :
					break
				else :
					if verbose :
						print('Error, unexpected pattern at the end!', parent, child, l, file=sys.stderr)
					if strict :
						raise FASTParseError( "Unexpected format at the end", id1, id2, paseq, caseq ) 
					else :
						return alignment

		'''
		eq = ''
		for p,c in zip( paseq, caseq ) :
			if p.isalpha() and c.isalpha() :
				eq += ':'
			else :
				eq += ' '
		'''

		alignment.raw_score = fast_score
		alignment.normalized_score = norm_score
		alignment.alignment1 = paseq[:-1]
		alignment.alignment2 = caseq[:-1]
		return alignment



#legacy code for parsing FAST alignment result.
head_pattern = re.compile( 'FAST ALIGNMENT: (.{4}).\S+\s(.{4}).\S+' )
score_pattern = re.compile( 'L=(\S+)\sSX=(\S+)\sSN=(\S+)\sL1=\S+\sL2=\S+\sRMSD=(\S+)' )
alignment_pattern = re.compile( ' [12]:\t(\S*)\**' )
def parse_fast_record( fp ) :
        #first line
        l = fp.readline()
        #case of no more record
        if not l:
                return None
        #print l
        match = head_pattern.match( l )
        #checking for starting point
        if not match :
                print('WARNING: wrong starting point!!' ,l, file=sys.stderr)
                while( not match ) :
                        l = fp.readline()
                        if not l :
                                print('End of file reached!', file=sys.stderr)
                                return None
                        print('following lines:', l, file=sys.stderr)
                        match = head_pattern.match( l )

        parent = match.group(1)
        child = match.group(2)
        #second line
        l = fp.readline()
        match = score_pattern.match( l )
        aligned_length = int(match.group(1))
        fast_score = match.group(2)
        norm_score = match.group(3)
        rmsd = match.group(4)

        if not aligned_length :
                l1 = fp.readline()
                l2 = fp.readline()
                if l1 == l2 == '\n' :
                        return parent, child, fast_score, rmsd, aligned_length, norm_score, None, None
                else:
                        print("Error! while parsing fast record. Aborted!", parent, child, file=sys.stderr)

        l = fp.readline()
        if l != '\n':
                print("Error, unexpected pattern at line3:", parent, child, l, file=sys.stderr)

        #alignment process
        paseq = ''
        caseq = ''
        while ( True ) :
                l = fp.readline()
                if l == '\n' :
                        break
                match = alignment_pattern.match( l )
                if not match :
                        break
                paseq += match.group(1)

                l = fp.readline()
                match = alignment_pattern.match( l )
                caseq += match.group(1)

                l = fp.readline()
                if l == '\n' :
                        pass
                else :
                        print("Error, unexpected pattern in alignment region!", parent, child, paseq, caseq, file=sys.stderr)


        #need to consume the trailing equvalence parts
        while ( True ) :
                l = fp.readline()
                if l == '\n' :
                        l = fp.readline()
                        if l == '\n' :
                                break
                        else :
                                print('Error, unexpected pattern at the end!', parent, child, l, file=sys.stderr)

        '''
        eq = ''
        for p,c in zip( paseq, caseq ) :
                if p.isalpha() and c.isalpha() :
                        eq += ':'
                else :
                        eq += ' '
        '''
        return (parent, child, fast_score, rmsd, aligned_length, norm_score, paseq[:-1], caseq[:-1] )

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


