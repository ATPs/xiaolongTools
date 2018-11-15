'''
This module contains codes to compare different alignments
and calculate some simple properties.

Current implementation contains,
1. simple overlap test.
2. advance overlap test using aligned regions.

One of the important usage of this module is 
detecting multidomain proteins.
'''


########################################################
#Main functions to compare overlaps in this module
########################################################


def calculate_simple_range_overlap( range1, range2 ) :
	'''
	Compares two ranges, range1 and range2, and returns the size of overlapped region.
	Both range1 and range2 are vectors of two integer values of ranges like in python slices.

	'''
	i, j = range1
	m, n = range2


	if n < i :
		return 0
	elif m > j :
		return 0

	#overlap start
	a = max(i,m)
	b = min(j,n)

	overlap_size = b-a
	if overlap_size < 0 :
		raise NegativeOverlapValueError( range1, range2, overlap_size )

	return overlap_size

def calculate_advanced_range_overlap_using_aligned_regions( aligned_pos1, aligned_pos2 ) :
	'''
	Compares two alignments using two lists of aligned positions, aligned_pos1 and aligned_pos2
	and returns number of commonly aligned positions in the two aligned position lists.
	
	Both aligned_pos1 and aligned_pos2 are lists of integer values.
	Each integer value means the index value of the aligned positions.
	'''

	s1 = set(aligned_pos1)
	s2 = set(aligned_pos2)
	
	return len( s1&s2 )



####################################################
#Easy Interface to the GenericAlignmentParser Objects.
####################################################


def get_query_overlap_value_between_two_methodrecs( methodrec1, methodrec2 ) :

	'''
	Using two PairwiseAlignmentMethodRecords,
	this function first run calculate_simple_range_overlap
	and if the range_overlap is more than half of the smaller one,
	it runs advance comparison to use aligned positions.
	'''

	range1 = ( methodrec1.start_position1, methodrec1.start_position1+methodrec1.get_sequence_length_in_alignment1() )
	range2 = ( methodrec2.start_position1, methodrec2.start_position1+methodrec2.get_sequence_length_in_alignment1() )

	length1 = range1[1]-range1[0]
	length2 = range2[1]-range2[0]

	#not a good measure
	#minlength = min(length1,length2)*1.0

	overlap1 = calculate_simple_range_overlap( range1, range2 )
	if overlap1 == 0 :
		return 0.0


	eq1 = methodrec1.get_equivalent_map()[0]
	eq2 = methodrec2.get_equivalent_map()[0]
	overlap2 = calculate_advanced_range_overlap_using_aligned_regions( eq1, eq2 )

	minlength = min(len(eq1),len(eq2))*1.0
	
	return overlap2/minlength


def get_query_overlap_value_between_two_aligned_postion_sets( aligned_pos1, aligned_pos2 ) :
	'''
	Two aligned positions will be compared and return overlap ratio
	scaled to the smaller aligned positions.
	'''
	minlength = min( len(aligned_pos1), len(aligned_pos2) ) * 1.0
	overlap = len(aligned_pos1 & aligned_pos2)

	return overlap/minlength


class MultidomainDetector :
	'''
	Detector of a multidomain protein.

	Here the main idea is to detect the region covered by two different
	SCOP superfamilies do not overlap, then the query domain might be a 
	multidomain protein.
	'''
	def __init__( self, alignment_record, query_id, alignment_type, homology_threshold, overlap_threshold, scop_sf_mapping ) :
		self.arec = alignment_record
		self.quid = query_id
		self.atyp = alignment_type
		self.hcut = homology_threshold
		self.ocut = overlap_threshold
		self.anno = scop_sf_mapping

		self.qusf = scop_sf_mapping[self.quid]
		self.aligned_position_mapping = {}
		self.found_non_overlapping_sfs = []

	def detect( self, type=None, lengths=None, structure_coverage_threshold=0.5 ) :
		if type == None :
			return self.detect_with_scop_groups()

		elif type == 'structure' :
			return self.detect_for_structure( structure_coverage_threshold, lengths )

	def detect_for_structure( self, structure_coverage_threshold, lengths ) :
		
		alignments = self.arec
		check_candidates = set([])
		found_unit_domains = {}
		domain_information = self.anno
		methodname = self.atyp
		score_threshold = self.hcut
		overlap_threshold = self.ocut

		for queryid in alignments :
			if not queryid in domain_information :
				continue

			queryrec = alignments[queryid]
			for hitid in queryrec :
				if not hitid in domain_information :
					continue

				hitrec = queryrec[hitid]
				methodrec = hitrec[methodname]

				if not methodrec.alignment1 :
					continue

				hit_length = lengths[methodrec.id2]
				aligned_length = len(methodrec.get_equivalent_map()[0])
				hit_coverage = aligned_length * 1.0 / hit_length

				if methodrec.normalized_score >= score_threshold and hit_coverage >= structure_coverage_threshold :
					#debugging
					#print hitid, methodrec.normalized_score
					check_candidates.add( methodrec )

		if not check_candidates :
			return

		query_id = alignments.get_query_ids()[0]
		query_annot = domain_information[query_id]

		check_candidates = list(check_candidates)
		for i, methodrec1 in enumerate(check_candidates) :

			for methodrec2 in check_candidates[i+1:] :
				annot1 = domain_information[ methodrec1.id2 ]
				annot2 = domain_information[ methodrec2.id2 ]

				#supposedly not in the same superfamily.
				if annot1 != annot2 :
					overlap_value = get_query_overlap_value_between_two_methodrecs( methodrec1, methodrec2 )
					if overlap_value < overlap_threshold :
						#print query_id, query_annot, 'region1:', methodrec1.id2, annot1, 'region2:', methodrec2.id2, annot2

						key = (min( annot1, annot2), max(annot1, annot2))
						if key in found_unit_domains :
							found_unit_domains[key].append((min(methodrec1.normalized_score, methodrec2.normalized_score), methodrec1, methodrec2 ))
						else :  
							found_unit_domains[key] = [(min(methodrec1.normalized_score, methodrec2.normalized_score), methodrec1, methodrec2 )]

		representatives = []
		#select best for each pair
		for k,v in found_unit_domains.items() :
			v.sort()
			representatives.append( v[-1] )
		representatives.sort()
		representatives.reverse()

		for s, mrec1, mrec2 in representatives :
			print(query_id, query_annot, ":", s, mrec1.id2, domain_information[mrec1.id2], mrec2.id2, domain_information[mrec2.id2])


	def detect_with_scop_groups( self ) :
		'''
		The default detector is using this relationship to the scop groups.
		'''
		queryrec = self.arec[self.quid]

		for hitid in queryrec :
			hitrec = queryrec[hitid]
			methodrec = hitrec[self.atyp]

			if not methodrec :
				continue
	
			if not methodrec.alignment1 :
				continue

			huid = hitrec.id2
			if not huid in self.anno :
				continue

			husf = self.anno[huid]

			"""
			if husf == self.qusf :
				continue
			"""

			if methodrec.normalized_score < self.hcut :
				continue

			aligned_positions = methodrec.get_equivalent_map()[0]

			if not husf in self.aligned_position_mapping :
				self.aligned_position_mapping[husf] = set()
				
			self.aligned_position_mapping[husf].update( aligned_positions )
				
		sf_list = list(self.aligned_position_mapping)
		for i, sfid1 in enumerate(sf_list) :
			for sfid2 in sf_list[i+1:] :
				overlap = get_query_overlap_value_between_two_aligned_postion_sets( self.aligned_position_mapping[sfid1], self.aligned_position_mapping[sfid2] )
				if overlap < self.ocut :
					self.found_non_overlapping_sfs.append( (overlap, sfid1, sfid2) )
					#print overlap
					#print self.aligned_position_mapping[sfid1]
					#print self.aligned_position_mapping[sfid2]
	
	def detect_with_scop_groups_and_coverage( self, hit_coverage_threshold, lengths ) :
		'''
		The default detector is using this relationship to the scop groups.
		'''
		queryrec = self.arec[self.quid]

		for hitid in queryrec :
			hitrec = queryrec[hitid]
			methodrec = hitrec[self.atyp]

			if not methodrec :
				continue
			if not methodrec.alignment1 :
				continue

			huid = hitrec.id2
			if not huid in self.anno :
				continue

			husf = self.anno[huid]

			"""
			if husf == self.qusf :
				continue
			"""

			if methodrec.normalized_score < self.hcut :
				continue

			aligned_positions = methodrec.get_equivalent_map()[0]
			hit_coverage = len(aligned_positions)*1.0/lengths[huid]

			if hit_coverage < hit_coverage_threshold :
				continue

			if not husf in self.aligned_position_mapping :
				self.aligned_position_mapping[husf] = set()
				
			self.aligned_position_mapping[husf].update( aligned_positions )
				
		sf_list = list(self.aligned_position_mapping)
		for i, sfid1 in enumerate(sf_list) :
			for sfid2 in sf_list[i+1:] :
				overlap = get_query_overlap_value_between_two_aligned_postion_sets( self.aligned_position_mapping[sfid1], self.aligned_position_mapping[sfid2] )
				if overlap < self.ocut :
					self.found_non_overlapping_sfs.append( (overlap, sfid1, sfid2) )
					#print overlap
					#print self.aligned_position_mapping[sfid1]
					#print self.aligned_position_mapping[sfid2]
						
	def print_detected_domains( self ) :
		if self.found_non_overlapping_sfs :
			print(self.quid, self.qusf, end=' ') 
			for (o, i, j) in self.found_non_overlapping_sfs :
				print(":", o, i, j, end=' ') 
			print()


class NegativeOverlapValueError( Exception ) :
	'''
	Raised when the overlap calculation is negative!
	'''
	pass
