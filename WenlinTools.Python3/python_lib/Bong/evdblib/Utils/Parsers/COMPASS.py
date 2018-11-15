
import sys
from math import log

from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentMethodRecord
from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentRecords, AlignmentAddingError
from evdblib.Utils import parse_profile_filename

verbose = 0
#special constant for tiny value
#when compass e-value hits zero!
tinyfloat = 2.2250738585072014e-308

class COMPASSParseError( Exception ) :
	pass

class COMPASS :
	def __init__( self, queryid=None, result_fn=None, iteration=None, use_iteration_info=True ) :
                '''
                COMPASS parses the COMPASS result file.
                Its major functionality is parsing the alignments and
                scores and return the GenericPairwiseAlignment object
                So that the alignments and scores can be easily manipulated.
                '''
                self.queryid = queryid
                self.result_fn = result_fn
                self.alignments = None
		self.iteration = iteration

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
                                raise COMPASSNullResultError( self.result_fn )

	def parse( self ) :
                '''
                parse the whole result file 
                and returns alignments.
                '''

		fp = open( self.result_fn )

		alignments = PairwiseAlignmentRecords()

		lines = []
		for line in fp :
			if line.startswith( 'Ali1:' ) :
				if lines :
					methodrecord = self.parse_record( lines ) 
					try :
						alignments.add( methodrecord )
					except AlignmentAddingError :
						print('WARNING: Method addition has a problem probably due to duplicated records.', file=sys.stderr)

					lines = []
			lines.append( line )
		else :
			methodrecord = self.parse_record( lines )
			try :
				alignments.add( methodrecord )
			except AlignmentAddingError :
				print("WARNING: Method addition has a problem probably due to duplicated records.", file=sys.stderr)

		return alignments


	def parse_record( self, content ) :
		'''
                parses a record in HHsearch result.
                '''

                aln1 = ''
                aln2 = ''
                start1 = 0
                start2 = 0

                additional_scores = {}

		if self.queryid :
			query_line_marker = self.queryid
			queryid = self.queryid
		else :
			query_line_marker = content[-4].split()[0]
			queryid = query_line_marker
		
		hitid = content[-2].split()[0]
		hit_line_marker = hitid

		if verbose :
			print("".join( content ))

		alignment_line_number = 0
		for i, l in enumerate(content) :
			if l.startswith( 'length1=' ) or l.startswith( 'Nseqs1=' ) :
				l = l.split()

				n1, s1 = l[1].split('=')
				additional_scores[n1] = float(s1)

				n2, s2 = l[3].split('=')
				additional_scores[n2] = float(s2)
				
			elif l.startswith( 'Smith-Waterman score' ) :
				l = l.split()
				raw_score = int(l[3])
				try :
					norm_score = -log(float( l[6] ))
				except ValueError :
					#chatching null results
					norm_score = -log(10000)

					methodname = 'compass'
					if self.iteration :
						methodname = methodname+'_'+str(self.iteration)
					methodrec = PairwiseAlignmentMethodRecord( 
						id1=queryid, id2=hitid,
						method_name = methodname,
						additional_scores = additional_scores,
						raw_score = raw_score,
						norm_score = norm_score)
					return methodrec

				except OverflowError :
					norm_score = -log(tinyfloat) #special case..

				alignment_line_number = i+1 #to offset the first null line
				break

		if not alignment_line_number :
			raise COMPASSParseError( "Error, Unexpected COMPASS out Format!" )

		for i, line in enumerate(content[alignment_line_number:]) :
			if verbose :
				print("Parsing line #", i, ":", line[:-1])
	
			if line == '\n' :
				continue
			elif (i%5)== 1 :
				l = line.split()
				if not aln1 :
					start1 = int(l[1]) -1 #0 based index!
					
				aln1 += l[-1]

				if not line.startswith( query_line_marker ) :
					print("WARNING: Query line marker not matching!", query_line_marker, l, file=sys.stderr)

			elif (i%5) == 2 :
				pass

			elif (i%5) == 3 :
				l = line.split()
				if not aln2 :
					start2 = int(l[1]) -1 #0 based index!
				aln2 += l[-1]

				if not line.startswith( hit_line_marker ) :
					print("WARNING: Hit line marker not matching!", hit_line_marker, l, file=sys.stderr)

			else :
				print("WARNING: You should not see this output!", file=sys.stderr)
				print(line[:-1], file=sys.stderr)

                methodname = 'compass'
                if self.iteration :
                        methodname = methodname+'_'+str(self.iteration)

                methodrec = PairwiseAlignmentMethodRecord( 
			id1=queryid, id2=hitid,
                        method_name = methodname,
                        alignment1 = aln1, alignment2 = aln2,
                        start_position1 = start1, start_position2 = start2,
                        additional_scores = additional_scores,
                        raw_score = raw_score,
                        norm_score = norm_score)

                #return ( evalue, probability, start1, start2, aln1, aln2 )
                return methodrec



#############################################################################
#COMPASS legacy codes
#The following codes should be used coding references and examples.
#Cannot be used as part of evdblib...
#############################################################################
#maximum value is treated as lowest score!! bug is fixed!
def parse_compass_result( content ) :
        new_content = []
        smith_waterman_score = 0
        negloge = 0

        in_flag = 0
        for l in content :
                if l[:23] == "Smith-Waterman score = " :
                        in_flag = 1
                        score_line = l.split()
                        try :
                                smith_waterman_score = float(score_line[3])
                        except :
                                print('WARNING! possible problem in compass raw score parsing', score_line, file=sys.stderr)
                                print(content, file=sys.stderr)
                                smith_waterman_score = 0.0

                        try :
                                evalue = float( score_line[6] )
                        except :
                                print("WARNING! possible problem in compass e-value parsing", score_line, file=sys.stderr)
                                print(content, file=sys.stderr)
                                negloge = -350.0
                                continue
                        
                        if evalue != 0.0 :
                                negloge = -math.log( evalue )
                        else :
                                negloge = 350
                                print('WARNING! e-value is zero!', score_line, file=sys.stderr)
                                continue

                        continue 

                if in_flag and (l[:11] == "Parameters:" or l[:6] == 'Ali1: ') :
                        in_flag = 0
                        break

                if in_flag and l != '\n' and l[0] != ' ' :
                        new_content.append( l )

        l = new_content[0].split()
        aln1 = l[2]
        start1 = int( l[1] ) -1 #correction for 0 based numbering

        l = new_content[1].split()
        aln2 = l[2]
        start2 = int( l[1] ) -1 #correction for 0 based numbering

        for i in range( 1, len(new_content)/2 ) :
                l = new_content[i*2].split()
                aln1 += l[1]

                l = new_content[i*2+1].split()
                aln2 += l[1]
        return ( smith_waterman_score, negloge, start1, start2, aln1, aln2 )


def get_compass_profile_info( fn ) :
        cmd = '/usr1/HorAServer/local/compass/compass_241_db1Xdb2'

        pp = os.popen( cmd + ' -i ' + fn + ' -j ' + fn )

        #I will parse the following two lines
        #length1=99      filtered_length1=98     length2=99      filtered_length2=98
        #Nseqs1=1        Neff1=9.929     Nseqs2=1        Neff2=9.929

        content = pp.readlines()
        for l in content :
                if l[:8] == 'length1=' :
                        l = l.split()
                        fl1 = int(l[1].split('=')[-1])  #filtered length1
                        fl2 = int(l[3].split('=')[-1])  #filtered length2
                        if fl1 == fl2 :
                                pass
                        else :
                                print("Filtered length1 and filtered length2 are different!", file=sys.stderr)
                                sys.stderr.writelines( content )
                                sys.exit()

                if l[:7] == 'Nseqs1=' :
                        l = l.split()
                        neff1 = float( l[1].split('=')[-1] ) #Neff1
                        neff2 = float( l[3].split('=')[-1] )

                        if neff1 == neff2 :
                                pass
                        else :
                                print("Neff1 and Neff2 are different!", file=sys.stderr)
                                sys.stderr, writelines( content )
                                sys.exit()
                        return fl1, neff1

def calculate_compass_evalue( alilen_mat1, dblen, n_eff1, n_eff2, score_final ) :

        lam_p1 = [0.12961, 0.00587, -0.00054];
        lam_p2 = [1.7977, 0.01421];
        k_p1 = [-0.01033];
        k_p2 = [10.46486, 0.05784, -0.02049];
        #lambda_est, K_est; #/* estimated lambda_g and K_g values */

        nefterm = 0.5*(n_eff1 + n_eff2);
        lenterm = 0.5*(1.0/dblen + 1.0/alilen_mat1);

        lamest1 = lam_p1[0] + lam_p1[1]*nefterm + lam_p1[2]*nefterm*nefterm;
        lamest2 = lam_p2[0] + lam_p2[1]*nefterm;
        lambda_est = lamest1 + lamest2*math.sqrt(lenterm);

        kest1 = k_p1[0];
        kest2 = k_p2[0] + k_p2[1]*nefterm + k_p2[2]*nefterm*nefterm;
        K_est = kest1 + kest2*lenterm;
        if(K_est<0) :
                K_est=0.002;

        Evalue = K_est*alilen_mat1*dblen*math.exp(-1.0*score_final*lambda_est);

        return Evalue


