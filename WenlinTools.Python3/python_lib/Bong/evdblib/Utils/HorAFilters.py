#!/usr/bin/env python

#this is part of Hua's negative filter
#supposedly read in three pairs of alignments
#and return maximum of consensus number between three equivalences
#divided by shorter length of two domains.
def consensus_number( eq1, eq2, eq3, length1, length2 ) :
        result = []
        for i, eqa in enumerate((eq1, eq2, eq3)) : 
                for j, eqb in enumerate( (eq1, eq2, eq3) ) :
                        if i < j : #this setting only runs upper triangle between eq's
                                result.append( compare_two_alignments( eqa, eqb ) ) 

        #print result

        return max( result )*1.0 / min( length1, length2 )
#however there's some question in this consensus calculation,
#because the consensus reuquires at least two of three methods should have some structural alignment.
#sometimes it is not quite right.. I guess "comment by bong-hyun"

#two arguments are equivalent maps between two alignments
#from the same sequences
#this function returns # of equivalent positions between two alignments
def compare_two_alignments( eq1, eq2 ) : 
        if eq1 and eq2 :
                pass
        else :
                return 0

        i = 0; #starting index of eq1
        j = 0 ; #starting index of eq2
        count = 0; #variable for counting same
        while i < len(eq1[0]) and j < len(eq2[0]) :
                if eq1[0][i] == eq2[0][j] : #checking parent side
                        if eq1[1][i] == eq2[1][j] :
                                count += 1
                                #print i, eq1[0][i], j, eq2[0][j]

                        i += 1;
                        j += 1;

                elif eq1[0][i] < eq2[0][j] :
                        i += 1
                elif eq1[0][i] > eq2[0][j] :
                        j += 1

                else :
                        print("This should not be shown", file=sys.stderr)
                        print("in comapare_two_alignments in hora_util.py", file=sys.stderr)
                        sys.exit(-1)
        return count


#This function is implementation of Hua's alignment post-processing 
def equivalent_post_process( eq1, eq2, sec1, sec2 ) :
        neq1 = []
        neq2 = []
        for peq, ceq in zip( eq1, eq2 ) :
                if (sec1[peq] == sec2[ceq] ) or (sec1[peq] == 'L' or sec2[ceq] == 'L') :
                        neq1.append( peq )
                        neq2.append( ceq )

        return neq1, neq2


#Given equivalency mapping and coordinates
#contact positions are counted
#  default sequence cutoff is 10  (two residues should be farther than 10 residues)
# default structural cutoff is 14 angstrom (two CA atoms should be closer than 14 angstrom)
def contact_number( eq1, eq2, ca1, ca2, seq_cutoff=10, str_cutoff=14 ) :
        count = 0
        for i in range( len(eq1) ) :
                for j in range( i+1, len(eq1) ) :
                        d1ij = dist( ca1[eq1[i]], ca1[eq1[j]] )
                        d2ij = dist( ca2[eq2[i]], ca2[eq2[j]] )
                        sd1ij = eq1[j]-eq1[i]
                        sd2ij = eq2[j]-eq2[i]
                        if ( sd1ij >= seq_cutoff and sd2ij >= seq_cutoff ) :
                                if d1ij <= str_cutoff and d2ij <= str_cutoff :
                                        count += 1

        return count

#distance calculation
import math
def dist( a, b ) :
        return math.sqrt( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) )


###############################################
###############################################
## Negative filter value calculation
##
## The following code is directly copied and translated from
## Hua's code
###############################################
###############################################

#  Filter scop pairs by two lines and one ellipse
#  line1: contact = s * consensus + t
#  line2: contact = m * consensus + n
#  ellipse: contact = 1 + (r2/r1) * [(r1+consensus-1)(r1-consensus+1)]^0.5    Note: ellipse is centered at (1,1)
#  Suppose that a pair has (consensus, contact)=(x,y), 
#  we plug x in the three equations and calculate the threshold values for contact y1, y2, y3. The final threshold y0 = max(y1, y2, y3).
#  If y > y0, the pair passes the filter.
#  r1>0 and r2>0 are the x and y radii of the ellipse;
#  (p1, p2) are (x, y) coordinates of the point at which the two lines intersect
#  b is the y-intercept of the first line
#  c is the x-intercept of the second line
#  s = (p2-b)/p1; t = b; m = p2/(p1-c); n = -c*p2/(p1-c)
#  Usage: filter.pl datafile contactcolumn r1 r2 p1 p2 b c

def negative_filter_value( consensus, contact, r1=55, r2=180, p1=30, p2=30, b=420, c=100 ) :
        #Hua's input
        #($datafile, $contactcolumn, $r1, $r2, $p1, $p2, $b, $c)=@ARGV;
        #function argument default values are from Hua's suggestion. 
        #seems not written in paper. Probably needs to be written in suplement materials

        s=(p2-b)*1.0/p1;
        t=b; m=p2*1.0/(p1-c);
        n=-1.0*c*p2/(p1-c);

        x = 100.0 * consensus + 1;
        y = contact + 1.0;
        y1 = s*x + t;
        y2 = m*x + n;

        if x <= 1.0 + r1:  # otherwise, the value in square root is negative
                y3 = 1.0 + (r2*1.0/r1)*((r1+x-1.0)*(r1-x+1.0))**0.5;
        else :
                y3 = 0.0;


        if y > y1 and y > y2 and y > y3 :
                return 1
        else :
                return 0


def get_hhsearch_probability( hitrec  ) :
	'''
	initial version of hhsearch probability extractor.
	from hitrec, PairwiseAlignmentHirRecord object.
	'''
	for i in range( 8, 0, -1 ) :
		hhrec = hitrec[ 'hhsearch_%d'%i ]
		if hhrec :
			return float(hhrec.get_normalized_score() )


def get_highest_hhsearch_probability( hitrec ) :
	'''
	Tries to extract highest HHsearch probability
	in the hitrec record.
	If no HHsearch records are found,
	return 0.0
	'''
	probs = []
	for methodname in hitrec :
		if methodname.startswith( 'hhsearch' ) :
			hhrec = hitrec[methodname]
			probs.append( float(hhrec.get_normalized_score()) )
	return max(probs)


def positive_filter( hitrec, hhsearch_threshold=0.9 ) :
	'''	
	returns HoraServer positive filter conclusion, hhsearch probability
	The conclusion is simply the hhsearch_probability is higher
	than the 
	'''
	hhsearch_probability = get_highest_hhsearch_probability( hitrec )
	if hhsearch_probability >= hhsearch_threshold :
		positive_filter_value = 1
	else :  
		positive_filter_value = 0

	return positive_filter_value, hhsearch_probability


def negative_filter( hitrec, ssd1, coord1, ssd2, coord2 ) :
	'''
	returns HorAServer negative filter conclusion, consensus value, contact value.
	'''

	dalirec = hitrec['dali']
	fastrec = hitrec['fast']
	tmalignrec = hitrec['tmalign']

	length1 = None
	length2 = None
	if dalirec and dalirec.alignment1 :
		dalieq = dalirec.get_equivalent_map()
		selected_eq_map = dalieq
		length1 = dalirec.get_sequence_length_in_alignment1()
		length2 = dalirec.get_sequence_length_in_alignment2()
	else :
		dalieq = None

	if fastrec and fastrec.alignment1 :
		fasteq = fastrec.get_equivalent_map()
		if not dalieq :
			selected_eq_map = fasteq
			length1 = fastrec.get_sequence_length_in_alignment1()
			length2 = fastrec.get_sequence_length_in_alignment2()
	else :
		fasteq = None

	if tmalignrec and tmalignrec.alignment1 :
		tmaligneq = tmalignrec.get_equivalent_map()
		if not dalieq and not fasteq :
			selected_eq_map = tmaligneq
			length1 = tmalignrec.get_sequence_length_in_alignment1()
			length2 = tmalignrec.get_sequence_length_in_alignment2()
	else :
		tmaligneq = None

	consensus = consensus_number( dalieq, fasteq, tmaligneq, length1, length2 )
	processed_eq1, processed_eq2 = equivalent_post_process( selected_eq_map[0], selected_eq_map[1], ssd1, ssd2 )

	contact = contact_number( processed_eq1, processed_eq2, coord1, coord2 )
	negative_value = negative_filter_value( consensus, contact )

	return negative_value, consensus, contact


if __name__ == '__main__' :
	'''
	This main function is used as an working example of
	HorAServer Negative & Positive filters.
	'''
	import sys, os
	sys.path.append( os.path.abspath(os.path.join(  __file__, '..', '..','..' )) )
	#print sys.path
	
	from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentRecords
	from evdblib.Utils.Parsers.SimpleParsers import parse_palsse_output, parse_coordinates_from_ca_only_pdb

	parecord = PairwiseAlignmentRecords()
	parecord.parse( sys.argv[1] )

	def get_secondary_structure( id, 
		ldir="/local_scratch/bhk/pssd/", 
		gdir='/home/bhk/projects/final_scop_update/results/pssd/' 
		) :

		if os.path.exists(ldir + '%s.ssd' % id) :
			return parse_palsse_output( ldir + '%s.ssd'%id )
		else :
			return parse_palsse_output( gdir + '%s.ssd'%id )

	def get_coordinates( 
		id, 
		ldir='/local_scratch/bhk/ca_structure/', 
		gdir='/home/bhk/projects/final_scop_update/structure_dir/', suffix='.pdb' 
		) :
		'''
		get coordinates from local scratch dir or global structure dir
		'''

		if os.path.exists( ldir ) :
			dir = ldir
		else :
			dir = gdir
		
		fn = os.path.join( dir, id ) + suffix 
		return parse_coordinates_from_ca_only_pdb( fn )
	
	for queryid in parecord :
		qrec = parecord[queryid]
		id1 = queryid
		ssd1 = get_secondary_structure( id1 )
		coord1 = get_coordinates( id1 )

		for hitid in qrec :
			hitrec = qrec[hitid]
			id2 = hitid
			ssd2 = get_secondary_structure( id2 )
			coord2 = get_coordinates( id2 )
			
			negative_filter_result, consensus, contact = negative_filter( hitrec, ssd1, coord1, ssd2, coord2 )
        		print("Negative_filter: %d consensus: %f contact: %d # %s %s" % (negative_filter_result, consensus, contact, id1, id2))
			
			##############################
			#positive filter
			##############################
			positive_filter_result, hhsearch_probability = positive_filter( hitrec )
			print("Positive_filter: %d hhsearch_probability: %f # %s %s\n" % (positive_filter_result, hhsearch_probability, id1, id2))

