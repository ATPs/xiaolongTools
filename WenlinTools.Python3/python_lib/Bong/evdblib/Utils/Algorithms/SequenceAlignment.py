import sys, io

verbose = 0

def aascore( aa1, aa2, match=1, penaltyA=-10000, penaltyX=0, gap_penalty=-1 ) :
	'''return matching score for the pair of AAs.
	'''
	aa1 = aa1.upper() 
	aa2 = aa2.upper()
	if not aa1.isalpha() or not aa2.isalpha() :
		return gap_penalty
	if aa1 == aa2 :
		return 1
	elif aa1=='X' or aa2=='X' :
		return penaltyX 
	else :
		return penaltyA
	
def align_sequences( seq1, seq2 ) :
	'''This function globally align two sequences using simple identity matrix.
	identical residue are +1 and mismatch is very high penalty : -1000, except for X
	which is promiscuous and can be aligned to any letter without any penalty.
	
	This function is designed to detect indels without misaligning amino acids.
	'''
	penaltyA  = - len(seq1) - len(seq2) #penalty for X
	gap = '-'
	
	#score_mat = [ [0]*(len(seq2)+1) for i in xrange(len(seq1)+1) ]
	#trace_mat = [ [0]*(len(seq2)+1) for i in xrange(len(seq1)+1) ]
	score_mat, trace_mat = {}, {}
	for i in range( len(seq1)+1 ) :
		for j in range( len(seq2)+1 ) :
			score_mat[ i,j ] = 0
			trace_mat[ i,j ] = ''
	#can be arbitrarily set.
	#but open penalty should be severer than the extention panalty -1.
	gap_open_penalty = -5 
	#filling the matrix
	for i in range( len(seq1)+1 ) :
		if i > 0 :
			aa1 = seq1[i-1]
		for j in range( len(seq2)+1 ) :
			if j > 0 :
				aa2 = seq2[j-1]
			#case for the origine
			if i == 0 and j == 0 :
				trace_mat[i,j] = '00'
				score_mat[i,j] = 0
				continue
			#all insertion in y
			elif j == 0 :
				trace_mat[i,j] = '10'
				score_mat[i,j] = score_mat[i-1,j] + aascore( aa1, gap ) + gap_open_penalty
				continue
			#all insertion in x
			elif i == 0 :
				trace_mat[i,j] = '01'
				score_mat[i,j] = score_mat[i,j-1] + aascore( gap, aa2 ) + gap_open_penalty
				continue
			temp = []
			#diagonal
			trace = '11'
			score = score_mat[i-1,j-1] + aascore( aa1, aa2, penaltyA=penaltyA )
			temp.append( (score, trace) )
			#x
			trace = '01'
			score = score_mat[i,j-1] + aascore( aa1, gap )
			if trace_mat[i,j-1] != trace :
				score = score + gap_open_penalty
			temp.append( (score, trace) )
			#y
			trace = '10'
			score = score_mat[i-1,j] + aascore( gap, aa2 )
			if trace_mat[i-1,j] != trace :
				score = score + gap_open_penalty
			temp.append( (score, trace) )
			
			temp.sort()
			score_mat[i,j], trace_mat[i,j] = temp[-1] #the biggest score
	#print score_mat[-1,-1]
	final_score = score_mat[ len(seq1),len(seq2) ]
			
	#retrace from the end
	fp1 = io.StringIO()
	fp2 = io.StringIO()
	i = len(seq1)
	j = len(seq2)
	while i-1 >= 0 or j-1 >= 0 :
		aa1 = seq1[i-1]
		aa2 = seq2[j-1]
		if trace_mat[i,j] == '11' :
			fp1.write( aa1 )
			fp2.write( aa2 )
			i = i - 1
			j = j - 1
		elif trace_mat[i,j] == '10' :
			fp1.write( aa1 )
			fp2.write( gap )
			i = i-1 
			j = j
		elif trace_mat[i,j] == '01' :
			fp1.write( gap )
			fp2.write( aa2 )
			i = i 
			j = j-1
		elif trace_mat[i,j] == '00' :
			break 
			
	aln1 = [ a for a in fp1.getvalue() ]
	aln1.reverse() 
	aln1 = ''.join(aln1)
	aln2 = [ a for a in fp2.getvalue() ]
	aln2.reverse() 
	aln2 = ''.join(aln2)
	return aln1, aln2, final_score

def align_sequences_with_affine_gap( seq1, seq2, aascore=aascore ) :
        '''This function globally align two sequences using simple identity matrix.
        identical residue are +1 and mismatch is very high penalty : -1000, except for X
        which is promiscuous and can be aligned to any letter without any penalty.
        
        This function is designed to detect indels without misaligning amino acids.
        '''
        len1 = len(seq1) + 1
        len2 = len(seq2 ) + 1
        penaltyA  = - len1 - len2 #penalty for X
        gap = '-'
        #can be arbitrarily set.
        #but open penalty should be severer than the extention panalty -1.
        gap_open_penalty = -2
        lowest_number = ( penaltyA + gap_open_penalty )* len1*len2
        #score_mat = [ [0]*(len(seq2)+1) for i in xrange(len(seq1)+1) ]
        #trace_mat = [ [0]*(len(seq2)+1) for i in xrange(len(seq1)+1) ]
        M, T, Ix, Tx, Iy, Ty = {}, {}, {}, {}, {}, {}
        for i in range( len1 ) :
                for j in range( len2 ) :
                        #initialize all score matrices with very negative number
                        #so that the unused elements are not visited!!
                        M[ i,j ] = lowest_number
                        T[ i,j ] = ' '
                        Ix[i,j ] = lowest_number
                        Tx[i,j ] = ' '
                        Iy[i,j ] = lowest_number
                        Ty[i,j ] = ' '
        #initialize 0,0 points
        T[0,0] = '00'
        M[0,0] = 0
        Tx[0,0] = '00'
        Ix[0,0] = 0
        Ty[0,0] = '00'
        Iy[0,0] = 0
        #End gap panelty setting from the starting point!
        #if you don't want to penalize end gaps
        #make the gap_open_penalty and aascore(*, gap) to 0
        #initialize '0th' insertion column of 'X' or sequence1
        j=0
        for i in range( 1, len1 ) :
                aa1 = seq1[i-1]
                if i == 1 :
                        Tx[i,j] = '00'
                        Ix[i,j] = M[i-1,j] + gap_open_penalty
                else :
                        Tx[i,j] = 'Ix'
                        Ix[i,j] = Ix[i-1,j] + aascore( aa1, gap )
        #initialize '0th' insertion row of 'Y' or sequence2
        i = 0
        for j in range( 1, len2 ) :
                aa2 = seq2[j-1]
                if j == 1 :
                        Ty[i,j] = '00'
                        Iy[i,j] = M[i,j-1] + gap_open_penalty
                else :
                        Ty[i,j] = 'Iy'
                        Iy[i,j] = Iy[i,j-1] + aascore( gap, aa2 )

        #filling the matrix
        for i in range( 1, len1 ) :
                aa1 = seq1[i-1]
                for j in range( 1, len2 ) :
                        aa2 = seq2[j-1]
                        '''
                        If you don't want to penalize end gap penalty,
                        you need to remove the right end gap penalties.
                        This is little bit more difficult than remove left end gap penalties.
                        
                        #you need to implement the following concept.
                        if i or j reached end of sequence the rest of the sequence should be filled with no penalties.
                        '''
                        temp = []
                        #getting the M, T matrix values
                        #diagonal
                        trace = 'M'
                        score = M[i-1,j-1] + aascore( aa1, aa2, penaltyA=penaltyA )
                        temp.append( (score, trace) )
                        #x
                        trace = 'Ix'
                        score = Ix[i-1,j-1] + aascore( aa1, aa2, penaltyA=penaltyA )
                        temp.append( (score, trace) )
                        #y
                        trace = 'Iy'
                        score = Iy[i-1,j-1] + aascore( aa1, aa2, penaltyA=penaltyA )
                        temp.append( (score, trace) )
                        temp.sort()
                        M[i,j], T[i,j] = temp[-1] #the biggest score
                        temp = []
                        #Ix, Tx matrices
                        trace = 'M'
                        score = M[i-1,j] + gap_open_penalty
                        temp.append( (score, trace) )
                        #Ix, Tx matrices
                        trace = 'Ix'
                        score = Ix[i-1,j] + aascore( aa1, gap )
                        temp.append( (score, trace) )
                        trace = 'Iy'
                        score = Iy[i-1,j] + gap_open_penalty
                        temp.append( (score, trace) )
                        temp.sort()
                        Ix[i,j], Tx[i,j] = temp[-1] #the biggest score
                        temp = []
                        #Iy, Ty matrices
                        trace = 'M'
                        score = M[i,j-1] + gap_open_penalty
                        temp.append( (score, trace) )
                        #Ix, Tx matrices
                        trace = 'Ix'
                        score = Ix[i,j-1] + gap_open_penalty
                        temp.append( (score, trace) )
                        trace = 'Iy'
                        score = Iy[i,j-1] + aascore( gap, aa2 )
                        temp.append( (score, trace) )
                        temp.sort()
                        Iy[i,j], Ty[i,j] = temp[-1] #the biggest score
        temp = [ (M[ len1-1, len2-1 ], T,  'M'),
               (Ix[len1-1, len2-1 ],   Tx, 'Ix'),
               (Iy[len1-1, len2-1 ],   Ty, 'Iy') ]
        temp.sort()
	#print temp[-1][:2]
        final_score, trace_mat, status = temp[ -1 ]
        #retrace from the end
        fp1 = io.StringIO()
        fp2 = io.StringIO()
        i = len1-1
        j = len2-1
        while i >= 0 and j >= 0 :
                if i > 0  : aa1 = seq1[i-1]
                else : aa1 = ''
                if j > 0 : aa2 = seq2[j-1]
                else : aa2 = ''
                if verbose :
                        print('status:', status, 'i:', i, 'j:', j, aa1, aa2)
                if not (i,j) in trace_mat :
                        print("(i,j)=(%d,%d)"%(i,j), "not found in the matrix", file=sys.stderr)
                        break
                current_status = status
                status = trace_mat[i,j] #setting up next status
                #trace_mat and status should be delivered from the above code
                if current_status == 'M' :
                        fp1.write( aa1 )
                        fp2.write( aa2 )
                        i = i - 1
                        j = j - 1
                elif current_status == 'Iy' :
                        fp1.write( gap )
                        fp2.write( aa2 )
                        i = i
                        j = j-1
                elif current_status == 'Ix' :
                        fp1.write( aa1 )
                        fp2.write( gap )
                        i = i-1
                        j = j
                if status == 'M' : trace_mat = T
                elif status == 'Ix' : trace_mat = Tx
                elif status == 'Iy' : trace_mat = Ty
                elif status == '00' :
                        break
                else :
                        print("Error! This message should not be seen", i,j, status, file=sys.stderr)
        else :
                print("Error!", i,j, status, file=sys.stderr)
        aln1 = [ a for a in fp1.getvalue() ]
        aln1.reverse()
        aln1 = ''.join(aln1)
        aln2 = [ a for a in fp2.getvalue() ]
        aln2.reverse()
        aln2 = ''.join(aln2)
        return aln1, aln2, final_score

def align_sequences_with_affine_gap_position_specific_gap_opening_penalty( seq1, seq2, aascore=aascore, chain_continuous=[] ) :
        '''This function globally align two sequences using simple identity matrix.
        identical residue are +1 and mismatch is very high penalty : -1000, except for X
        which is promiscuous and can be aligned to any letter without any penalty.
        
        This function is designed to detect indels without misaligning amino acids.
	Note that the sequence2 gap opening penalty is now became 
	position dependent. This garantees to have correct insertions at chain breaking positions
	for PDB seqres mapping and missing residue positioning.

	chain_continuous contains list of boolean values of continuity between seq2[i], seq2[i+1] positions. 
	So by definition this continiutiy list should be 1 residue shorter than seq2.
	Seq2 is supposed to be "atom sequence" whereas seq1 is suppoed to be SEQRES sequence.
        '''
        len1 = len(seq1) + 1
        len2 = len(seq2 ) + 1
        penaltyA  = - len1 - len2 #penalty for X
        gap = '-'
        #can be arbitrarily set.
        #but open penalty should be severer than the extention panalty -1.
        gap_open_penalty = -2

	if chain_continuous :
		gapopen2 = []
		for b in chain_continuous :
			if b :
				gapopen2.append( gap_open_penalty * 2 )
			else :
				gapopen2.append( gap_open_penalty )
		gapopen2.append( gap_open_penalty )

		if len(gapopen2) != len(seq2 ) :
			if verbose :
				print("len1:", len(seq1), file=sys.stderr)
				print("len2:", len(seq2), file=sys.stderr)
				print("gapopen2", len(gapopen2), file=sys.stderr)
			raise IndexError( "The gapopen2 should have same length as sequence2." )

		if verbose :
			for pen, aa in zip( gapopen2, seq2 ) :
				print(pen, aa, file=sys.stderr)

	else :
		gapopen2 = [gap_open_penalty]*len(seq2)
		
        lowest_number = ( penaltyA + gap_open_penalty )* len1*len2
        #score_mat = [ [0]*(len(seq2)+1) for i in xrange(len(seq1)+1) ]
        #trace_mat = [ [0]*(len(seq2)+1) for i in xrange(len(seq1)+1) ]
        M, T, Ix, Tx, Iy, Ty = {}, {}, {}, {}, {}, {}
        for i in range( len1 ) :
                for j in range( len2 ) :
                        #initialize all score matrices with very negative number
                        #so that the unused elements are not visited!!
                        M[ i,j ] = lowest_number
                        T[ i,j ] = ' '
                        Ix[i,j ] = lowest_number
                        Tx[i,j ] = ' '
                        Iy[i,j ] = lowest_number
                        Ty[i,j ] = ' '
        #initialize 0,0 points
        T[0,0] = '00'
        M[0,0] = 0
        Tx[0,0] = '00'
        Ix[0,0] = 0
        Ty[0,0] = '00'
        Iy[0,0] = 0
        #End gap panelty setting from the starting point!
        #if you don't want to penalize end gaps
        #make the gap_open_penalty and aascore(*, gap) to 0
        #initialize '0th' insertion column of 'X' or sequence1
        j=0
        for i in range( 1, len1 ) :
                aa1 = seq1[i-1]
                if i == 1 :
                        Tx[i,j] = '00'
                        Ix[i,j] = M[i-1,j] + gap_open_penalty
                else :
                        Tx[i,j] = 'Ix'
                        Ix[i,j] = Ix[i-1,j] + aascore( aa1, gap )
        #initialize '0th' insertion row of 'Y' or sequence2
        i = 0
        for j in range( 1, len2 ) :
                aa2 = seq2[j-1]
                if j == 1 :
                        Ty[i,j] = '00'
                        Iy[i,j] = M[i,j-1] + gap_open_penalty
                else :
                        Ty[i,j] = 'Iy'
                        Iy[i,j] = Iy[i,j-1] + aascore( gap, aa2 )

        #filling the matrix
        for i in range( 1, len1 ) :
                aa1 = seq1[i-1]
                for j in range( 1, len2 ) :
                        aa2 = seq2[j-1]
                        '''
                        If you don't want to penalize end gap penalty,
                        you need to remove the right end gap penalties.
                        This is little bit more difficult than remove left end gap penalties.
                        
                        #you need to implement the following concept.
                        if i or j reached end of sequence the rest of the sequence should be filled with no penalties.
                        '''
                        temp = []
                        #getting the M, T matrix values
                        #diagonal
                        trace = 'M'
                        score = M[i-1,j-1] + aascore( aa1, aa2, penaltyA=penaltyA )
                        temp.append( (score, trace) )
                        #x
                        trace = 'Ix'
                        score = Ix[i-1,j-1] + aascore( aa1, aa2, penaltyA=penaltyA )
                        temp.append( (score, trace) )
                        #y
                        trace = 'Iy'
                        score = Iy[i-1,j-1] + aascore( aa1, aa2, penaltyA=penaltyA )
                        temp.append( (score, trace) )
                        temp.sort()
                        M[i,j], T[i,j] = temp[-1] #the biggest score

                        temp = []
                        #Ix, Tx matrices
                        trace = 'M'
                        score = M[i-1,j] + gapopen2[j-1]
                        temp.append( (score, trace) )
                        #Ix, Tx matrices
                        trace = 'Ix'
                        score = Ix[i-1,j] + aascore( aa1, gap )
                        temp.append( (score, trace) )
                        trace = 'Iy'
                        score = Iy[i-1,j] + gapopen2[j-1]
                        temp.append( (score, trace) )
                        temp.sort()
                        Ix[i,j], Tx[i,j] = temp[-1] #the biggest score

                        temp = []
                        #Iy, Ty matrices
                        trace = 'M'
                        score = M[i,j-1] + gap_open_penalty
                        temp.append( (score, trace) )
                        #Ix, Tx matrices
                        trace = 'Ix'
                        score = Ix[i,j-1] + gap_open_penalty
                        temp.append( (score, trace) )
                        trace = 'Iy'
                        score = Iy[i,j-1] + aascore( gap, aa2 )
                        temp.append( (score, trace) )
                        temp.sort()
                        Iy[i,j], Ty[i,j] = temp[-1] #the biggest score


        temp = [ (M[ len1-1, len2-1 ], T,  'M'),
               (Ix[len1-1, len2-1 ],   Tx, 'Ix'),
               (Iy[len1-1, len2-1 ],   Ty, 'Iy') ]
        temp.sort()
	#print temp[-1][:2]
        final_score, trace_mat, status = temp[ -1 ]
        #retrace from the end
        fp1 = io.StringIO()
        fp2 = io.StringIO()
        i = len1-1
        j = len2-1
        while i >= 0 and j >= 0 :
                if i > 0  : aa1 = seq1[i-1]
                else : aa1 = ''
                if j > 0 : aa2 = seq2[j-1]
                else : aa2 = ''
                if verbose :
                        print('status:', status, 'i:', i, 'j:', j, aa1, aa2)
                if not (i,j) in trace_mat :
                        print("(i,j)=(%d,%d)"%(i,j), "not found in the matrix", file=sys.stderr)
                        break
                current_status = status
                status = trace_mat[i,j] #setting up next status
                #trace_mat and status should be delivered from the above code
                if current_status == 'M' :
                        fp1.write( aa1 )
                        fp2.write( aa2 )
                        i = i - 1
                        j = j - 1
                elif current_status == 'Iy' :
                        fp1.write( gap )
                        fp2.write( aa2 )
                        i = i
                        j = j-1
                elif current_status == 'Ix' :
                        fp1.write( aa1 )
                        fp2.write( gap )
                        i = i-1
                        j = j
                if status == 'M' : trace_mat = T
                elif status == 'Ix' : trace_mat = Tx
                elif status == 'Iy' : trace_mat = Ty
                elif status == '00' :
                        break
                else :
                        print("Error! This message should not be seen", i,j, status, file=sys.stderr)
        else :
                print("Error!", i,j, status, file=sys.stderr)
        aln1 = [ a for a in fp1.getvalue() ]
        aln1.reverse()
        aln1 = ''.join(aln1)
        aln2 = [ a for a in fp2.getvalue() ]
        aln2.reverse()
        aln2 = ''.join(aln2)
        return aln1, aln2, final_score


def print_matrix( m, len1, len2 ) :
        for i in range( len1 ):
                for j in range( len2 ) :
                        print('%5s'%str(m[i,j]), end=' ')
                print()
        print()

