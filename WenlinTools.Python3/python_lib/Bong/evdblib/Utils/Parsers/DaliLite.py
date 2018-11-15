
import os, sys, io, traceback

from .DaliLiteDAT import DaliLiteDAT
from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentMethodRecord

verbose = 0

class DaliLite :
	def __init__( self, dccp_fn='', dat_fn1='', dat_fn2='', query_id='', hit_id='', strict=False ) :
		'''
		Main parser class for DaliLite structural alignment
		result file, dccp.

		Note that to parse and build a DaliLite alignment correctly
		we need 
		'''

		self.dccp, self.dat1, self.dat2 = None, None, None
		self.query_id = query_id
		self.hit_id = hit_id
		self.strict = strict

		self.alignment = PairwiseAlignmentMethodRecord( id1=query_id, id2=hit_id, method_name='dali' )

		if dat_fn1 :
			self.dat1 = DaliLiteDAT( dat_fn1 )
			if not self.query_id :
				self.query_id = self.dat1.identifier


		if dat_fn2 :
			self.dat2 = DaliLiteDAT( dat_fn2 )
			if not self.hit_id :
				self.hit_id = self.dat2.identifier


		if dccp_fn and os.path.exists( dccp_fn ) :
			pass #good!
		else :
			return
			

		if self.hit_id and dccp_fn :
			self.dccp = DCCP( dccp_fn, self.hit_id, strict=self.strict )

		if verbose  :
			print("dccp_fn")
			print(dccp_fn)
			print("dccp")
			print(self.dccp.header_records)
			print(self.dccp.alignment_records)


		if self.dccp and self.dat1 and self.dat2 :
			self.alignment = self.parse(self.dccp, self.dat1, self.dat2, self.query_id, self.hit_id )

	def parse( self, dccp, dat1, dat2, query_id, hit_id ) :
		dali_alignment = self.dccp.build_best_daliz_alignment( self.dat1.sequence, self.dat2.sequence )
		if dali_alignment :
			dalirec = PairwiseAlignmentMethodRecord(
				method_name='dali',
				id1=query_id,
				id2=hit_id,
				alignment1 = dali_alignment[4],
				alignment2 = dali_alignment[5],
				raw_score= dali_alignment[0],
				norm_score=dali_alignment[1]
			)
		else :
			dalirec = PairwiseAlignmentMethodRecord(
				method_name= 'dali',
				id1=query_id,
				id2=hit_id )

		#check if the dalirec contains an alignment
		if dalirec.alignment1 :
			#checking sequence!
			if self.dat1.sequence == dalirec.seq1 and self.dat2.sequence == dalirec.seq2 :
				pass
			else :
				print("Error: Sequence does not match", dalirec.id1, dalirec.id2, file=sys.stderr)
				print('dat1_seq:', self.dat1.sequence, file=sys.stderr)
				print('aln1_seq:', dalirec.seq1, file=sys.stderr)
				print('dat1_seq:', self.dat1.sequence, file=sys.stderr)
				print('aln1_seq:', dalirec.seq1, file=sys.stderr)
				raise DCCPParsingError( "Parsed Alignment and Input Sequence does not match!" )

		return dalirec
		

class DCCPParsingError( Exception ) :
	pass

class DCCP :
	'''	
	DCCP is a raw DaliLite output parser
	and container for sub-optimal dali alignments.
	'''
	def __init__( self, dccp_fn='', code2='', strict=True ) :
		self.dccp_fn = dccp_fn
		self.code2 = code2
		self.alignment_records = None #after correct parsing this will be  a list
		self.header_records = None
		self.strict = strict

		if not code2 :
			raise DCCPParsingError( "code2 or hit identifier should be given!" )

		if dccp_fn :
			fp = open( dccp_fn )
			self.header_records, self.alignment_records = self.parse_dccp_records( fp )
			fp.close()

		if verbose :
			print("DCCP.__init__")
			print("self.code2:", self.code2)
			print("self.header_records", self.header_records)

	def __str__( self ) :
		fp = io.StringIO()
		for header, alignment in zip( self.header_records, self.alignment_records ) :
			print(header, file=fp)
			print(alignment, file=fp)

		return fp.getvalue()


	def parse_dccp_records( self, fp ) :
		'''
		reads dccp file descriptor fp
		and break down all dccp records into a list.

		Note that the method returns
		list of dccp string records.
		'''
		records = []
		headers = []
		str_fp = None
		for l in fp :
			if verbose :
				print(l)

			if l.startswith( ' DCCP' ) :
				if str_fp :
					records.append( str_fp.getvalue() )
					str_fp.close()
		
				str_fp = io.StringIO()

				header = self.parse_dccp_header( l )
				headers.append( header ) 

			elif l.startswith( ' alignment' ) :
				#nothing to be done!
				continue
			else :
				#should be alignment line!
				str_fp.write( l )
		else :
			if str_fp : 
				records.append( str_fp.getvalue() )

		return headers, records


	def parse_dccp_header( self, header ) :
		'''
		returns a dccp header.
		A header is a tuple of Z-score, raw_score, rmsd, switch_flag, child_id.
		The switch_flag indicates 
		if the alignmet record should be inverted or not. 
		The child_id is for error checking.
		'''
		#dccp_index += 1
		raw_score = float( header[9:18] )
		Z_score = float( header[27:34] )
		rmsd = float( header[19:23] )

		if  header[69:74] == self.code2 : #'mol2' :
			child_id = header[69:74]
			switch_flag = 1
		elif header[75:80] == self.code2 :
			child_id = header[75:80]
			switch_flag = 0
		else :
			raise DCCPParsingError( "Switching identification has problem!", header, self.code2 )
	

		return ( Z_score, raw_score, rmsd, switch_flag, child_id )


	def parse_alignment_line( self, line ) :
		'''
		Parses the alingment line.
	
		Simple string.split() should not be used
		for parsing lines in DCCP file.
		As like PDB format it is formatted with fixed string
		lengths.
		'''
		lengthrec = 10 #the alignment definition '????  ????' 4+2+4

		nrec = len( line ) / lengthrec #
		if len(line) != lengthrec*nrec :
			print("Error the alignment line has some strange bug!!", file=sys.stderr)
			print(line, file=sys.stderr)
			raise DCCPParsingError( line, "Alignment line parsing has problem!")

		parsed = []
		for i in range( nrec ) :
			temp = line[i*10:(i+1)*10]
			if verbose :
				print(line)
				print(temp)
			l1, l2 = temp.split()
			parsed.append( int(l1) )
			parsed.append( int(l2) )

		return parsed


	def parse_alignment_string( self, alignment_string, swap_flag ) :
		
		alignments = []
		for l in alignment_string.split('\n') :
			alignments.append( self.parse_alignment_line( l ) )

		parent, child = [], []
		if not swap_flag :
			for i in range( len(alignments) / 2 ) :
				for j in range( len( alignments[i] ) / 2 ) :
					parent.append( [alignments[i][j*2],alignments[i][j*2+1] ] )
			for i in range( len(alignments) / 2, len(alignments) ) :
				for j in range( len( alignments[i] ) / 2 ) :
					child.append( [alignments[i][j*2], alignments[i][j*2+1]] )
		else :
			for i in range( len(alignments) / 2, len( alignments ) ) :
				for j in range( len( alignments[i] ) / 2 ) :
					parent.append( [alignments[i][j*2],alignments[i][j*2+1] ] )
			for i in range( len(alignments) / 2 ) :
				for j in range( len( alignments[i] ) / 2 ) :
					child.append( [alignments[i][j*2], alignments[i][j*2+1]] )

                ##new checking!#
                #getting rid of the trailing zeros!
                #trailing zeros should be removed!
                while [0,0] in parent :
                        parent.remove( [0,0] )


                while [0,0] in child :
                        child.remove( [0,0] )

                #checking monotonic increase in the alignment length!
                for (b1,e1),(b2,e2) in zip(parent[:-1],parent[1:]) :
                        if b1 <= e1 < b2 <= e2 :
                                pass
                        else :
                                #print >>sys.stderr, "Error, non-monotonic increase in query residue indices! while parsing dccp file!"   
				raise DCCPParsingError( (b1, e1, b2, e2,), "Non-monotonic increase in query alignment!" )

                for (b1,e1),(b2,e2) in zip(child[:-1],child[1:]) :
                        if b1 <= e1 < b2 <= e2 :
                                pass
                        else :
                                #print >>sys.stderr, "Error, non-monotonic increase in query residue indices! while parsing dccp file!"   
				raise DCCPParsingError( (b1, e1, b2, e2,), "Non-monotonic increase in hit alignment!" )

                if len( parent ) == len( child ) :
                        pass
                else :
                        #print >>sys.stderr, "Error, while parsing dccp file!"
			raise DCCPParsingError( "Alignment lengths does not match up!" )


		return parent, child
		

	def build_best_daliz_alignment( self, seq1, seq2 ) :
		'''
		returns an alignment record of
		the select dccp record with best DaliLite Z-score
		(for ties, select higher raw-score).
		'''
		for header, alignment in reversed( sorted( zip(self.header_records,self.alignment_records) ) ) :
			try :
				aln1, aln2 = self.build_alignment( alignment, seq1, seq2, header[3] )
				return ( header[1], header[0], 0, 0, aln1, aln2 )

			except DCCPParsingError as e :
				print(repr(e), file=sys.stderr)
				continue

			except Exception as e:
				print(repr(e), file=sys.stderr)
				if self.strict :
					raise
				continue
		return None

	def build_alignment( self, alignment_string, seq1, seq2, swap_flag ) :
		def1, def2 = self.parse_alignment_string( alignment_string, swap_flag )

		aln_fp1 = io.StringIO()
		aln_fp2 = io.StringIO()
		
		eqcount1, eqcount2 = 0, 0 #counting equivalent positions from definition
		prep, prec = 0, 0
		for (ps,pe),(cs,ce) in zip( def1, def2 ):
			eqcount1 += -ps+pe+1 #accounting purposes
			eqcount2 += -cs+ce+1 #accounting purposes
			#adding gaps
			if (ps-cs-prep+prec) > 0 :
				aln_fp2.write( '-'*(ps-cs-prep+prec) )
			elif (cs-ps+prep-prec) > 0 :
				aln_fp1.write( '-'*(cs-ps+prep-prec) )
			#adding non-equivalent parts
			seq1_noneq = seq1[prep:ps-1].lower()
			seq2_noneq = seq2[prec:cs-1].lower()
			aln_fp1.write( seq1_noneq )
			aln_fp2.write( seq2_noneq )

			#equivalent parts
			seq1_eq = seq1[ps-1:pe].upper()
			seq2_eq = seq2[cs-1:ce].upper()
			aln_fp1.write( seq1_eq )
			aln_fp2.write( seq2_eq )

			if seq1_eq and seq2_eq :
				pass
			else :
				raise DCCPParsingError( "Equivalent Block building has problem!", (ps, pe), (cs, ce) )

			#saving position
			prep = pe
			prec = ce

		##
		##Adding last part beyond definition
		##

		ps = len(seq1)+1
		pe = len(seq1)+1
		cs = len(seq2)+1
		ce = len(seq2)+1
		#adding gaps
		if (ps-cs-prep+prec) > 0 :
			aln_fp2.write( '-'*(ps-cs-prep+prec) )
		elif (cs-ps+prep-prec) > 0 :
			aln_fp1.write( '-'*(cs-ps+prep-prec) )
		#adding non-equivalent parts
		aln_fp1.write( seq1[prep:ps-1].lower() )
		aln_fp2.write( seq2[prec:cs-1].lower() )
		#equivalent parts
		aln_fp1.write( seq1[ps-1:pe].upper() )
		aln_fp2.write( seq2[cs-1:ce].upper() )
		#saving position
		prep = pe
		prec = ce

		aln1 = aln_fp1.getvalue()
		aln2 = aln_fp2.getvalue()
		if len( aln1 ) != len( aln2 ) :
			#print >>sys.stderr, "Error! alignment lengths are not equal.", len( alignment[0]), len( alignment[1])
			raise DCCPParsingError( "Alignment lengths are not equal", len(aln1), len(aln2) )

		eq1count = 0
		eq2count = 0
		for a, b in zip(aln1, aln2 ) :
			if a.isupper() :
				eq1count += 1
			if b.isupper() :
				eq2count += 1

		if eq1count == eq2count == eqcount1 == eqcount2 :
			pass
		else :
			raise DCCPParsingError( "Equivalent position mapping is not right!", eqcount1, eqcount2, eq1count, eq2count )
				

		return aln1, aln2



##############################
#legacy codes
##############################
#for referencing purposes only!
#bug found!!!
def parse_dccp( fp, code2 ) :
        dccp_array = [] #temporary array to save dccp information
        raw_score_array = []
        switch_flag = 0

        dccp_index = 0  #index for max_z ---- 1 based!!! careful!!!
        max_index = 0

        l = fp.readline()
        while l :
                temp = []
                if l[1:5] == 'DCCP':
                        dccp_index += 1
                        raw_score = float( l[9:18] )
                        Z_score = float( l[27:34] )
                        rmsd = float( l[19:23] )
                        raw_score_array.append( raw_score )
                        if not max_index :
                                max_index = 1
                        elif Z_score > dccp_array[max_index-1][2] :
                                max_index = dccp_index

                        #previously
                        #rmsd < rmsd_array[max_index-1] :
                        #was the second condition

                        elif Z_score == dccp_array[max_index-1][2] and raw_score > raw_score_array[max_index-1]  :
                                max_index = dccp_index

                        if  l[69:73] == code2[:4] : #'mol2' :
                                child_id = l[69:74]
                                switch_flag = 1
                        else :  
                                child_id = l[75:80]
                                switch_flag = 0

                        temp.append( child_id )
                        temp.append( raw_score )
                        temp.append( Z_score )

                #read alignment
                l = fp.readline()
                l = fp.readline()
                parent = []
                child = []
                alignments = []
                while not l[1:5] == 'DCCP' :
                        #Bug!!! fixed length should be used
                        #l = l.split()
                        l = parse_alignment_line( l ) #fixed!!

                        for i in range( len(l) ) :
                                l[i] = int( l[i] )
                        alignments.append( l )
                        l = fp.readline()
                        #end of loop when the file is read
                        if not l :
                                break

                if not switch_flag :
                        for i in range( len(alignments) / 2 ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        parent.append( [alignments[i][j*2],alignments[i][j*2+1] ] )
                        for i in range( len(alignments) / 2, len(alignments) ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        child.append( [alignments[i][j*2], alignments[i][j*2+1]] )
                else :
                        for i in range( len(alignments) / 2, len( alignments ) ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        parent.append( [alignments[i][j*2],alignments[i][j*2+1] ] )
                        for i in range( len(alignments) / 2 ) :
                                for j in range( len( alignments[i] ) / 2 ) :
                                        child.append( [alignments[i][j*2], alignments[i][j*2+1]] )

                if len( parent ) == len( child ) :
                        pass
                else :
                        print("Error, while parsing dccp file!", file=sys.stderr)

                temp.append( parent )
                temp.append( child )
                dccp_array.append( temp )

        if not max_index :
                return None

        #debug
        #print dccp_array[max_index-1]
        return dccp_array[max_index-1]

        #returns only maximum Z score alignment information

############################################
# parse alignment definition in dccp file
# simple l.split() can have serious bug
# since the alignment definition is written 
# fixed length 
############################################
def parse_alignment_line( line ) :
        lengthrec = 10 # the alignment definition '????  ????' 4+2+4
        if line[-1] == '\n' :
                line = line[:-1]

        nrec = len( line ) / lengthrec #
        if len(line) != lengthrec*nrec :
                print("Error the alignment line has some strange bug!!", file=sys.stderr)
                print(line, file=sys.stderr)
                sys.exit()

        parsed = []
        for i in range( nrec ) :
                temp = line[i*10:(i+1)*10]
                l1, l2 = temp.split()
                parsed.append( l1 )
                parsed.append( l2 )

        return parsed

###########################################################
# Alignment building function
###########################################################
def build_alignment( seq1, seq2, def1, def2 ) :
        alignment = ['','']

        prep = 0
        prec = 0
        for (ps,pe),(cs,ce) in zip( def1, def2 ):
                #adding gaps
                if (ps-cs-prep+prec) > 0 :
                        alignment[1] += '-'*(ps-cs-prep+prec)
                elif (cs-ps+prep-prec) > 0 :
                        alignment[0] += '-'*(cs-ps+prep-prec)
                #adding non-equivalent parts
                alignment[0] += seq1[prep:ps-1].lower()
                alignment[1] += seq2[prec:cs-1].lower()
                #equivalent parts
                alignment[0] += seq1[ps-1:pe].upper()
                alignment[1] += seq2[cs-1:ce].upper()
                #saving position
                prep = pe
                prec = ce

        ##
        ##Adding last part beyond definition
        ##

        ps = len(seq1)+1
        pe = len(seq1)+1
        cs = len(seq2)+1
        ce = len(seq2)+1
        #adding gaps
        if (ps-cs-prep+prec) > 0 :
                alignment[1] += '-'*(ps-cs-prep+prec)
        elif (cs-ps+prep-prec) > 0 :
                alignment[0] += '-'*(cs-ps+prep-prec)
        #adding non-equivalent parts
        alignment[0] += seq1[prep:ps-1].lower()
        alignment[1] += seq2[prec:cs-1].lower()
        #equivalent parts
        alignment[0] += seq1[ps-1:pe].upper()
        alignment[1] += seq2[cs-1:ce].upper()
        #saving position
        prep = pe
        prec = ce

        if len( alignment[0] ) != len( alignment[1] ) :
                print("Error! alignment lengths are not equal.", len( alignment[0]), len( alignment[1]), file=sys.stderr)
                print(alignment[0], file=sys.stderr)
                print(alignment[1], file=sys.stderr)

        return alignment

