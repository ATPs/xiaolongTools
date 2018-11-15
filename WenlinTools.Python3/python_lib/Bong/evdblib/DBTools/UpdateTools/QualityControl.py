'''
This module contains codes to check the quality of DB 
or results, like profiles generated, alignments generated, etc.
'''
import os, sys, glob
from evdblib.DBTools import Settings
from evdblib.Utils import build_sequence_filename, parse_profile_filename

verbose = 0

def check_generated_profile( dominfo ) :
	'''
	returns number of iterations when all iteration results
	meet integrity check criteria.
	
	Currently, simply check
	the existence of the all iteration files have the same 
	number of iterations.
	A more rigorous checking should be implemented
	e.g. checking the contents of each profile type.
	'''

	if verbose :
		print('checking profile integrity...', dominfo['uniqueid'])

	if dominfo.get( 'profile_integrity' ) :
		return dominfo[ 'profile_integrity' ]

	domid = dominfo['uniqueid']
	dompath = dominfo['domain_path']

	if not os.path.exists( dompath ) :
		return 0

	basename = os.path.join( dompath, domid )
		
	all_profiles = glob.glob( os.path.join( basename + ".*.*" ) )
	pnp = [ fn for fn in all_profiles if fn.endswith( '.pnp' ) and parse_profile_filename(fn)[2] ]
	a3m = [ fn for fn in all_profiles if fn.endswith( '.a3m' ) and parse_profile_filename(fn)[2] ]
	hhm = [ fn for fn in all_profiles if fn.endswith( '.hhm' ) and parse_profile_filename(fn)[2] ]
	cnp = [ fn for fn in all_profiles if fn.endswith( '.cnp' ) and parse_profile_filename(fn)[2] ]
	cnplen = [ fn for fn in all_profiles if fn.endswith( '.cnp.len' ) and parse_profile_filename(fn)[2] ]

	#checking single item
	#for null case to be the True.
	if not pnp :
		if verbose :
			print("No result file is found!")

		return 0

	if len(pnp) == len(a3m) == len(hhm) == len(cnp) == len(cnplen) :
		#need to continue to test a little more on the maximum
		#iteration!!
		pass 
	else :
		return 0
	try :
		max_pnpi = max( [ int(parse_profile_filename(fn)[-2]) for fn in pnp ] )
		max_a3mi = max( [ int(parse_profile_filename(fn)[-2]) for fn in a3m ] )
		max_hhmi = max( [ int(parse_profile_filename(fn)[-2]) for fn in hhm ] )
		max_cnpi = max( [ int(parse_profile_filename(fn)[-2]) for fn in cnp ] )
		max_cnpleni = max( [ int(parse_profile_filename(fn)[-2]) for fn in cnplen ] )

		if max_pnpi == max_a3mi == max_hhmi == max_cnpi == max_cnpleni == len(pnp) :
			return max_pnpi
		else :
			return 0

	except ValueError :
		if verbose :
			print("The filename parsing routine has a problem!")
			print(pnp)
			print(a3m)
			print(hhm)
			print(cnp)
			print(cnplen)
		return 0


#checked generated profile is a little bit of misnomer.
#I added check profile integrity proxy function.
check_profile_integrity = check_generated_profile
	

def check_generated_profiles( domain_informations ) :
	'''
	Check each domain information dictionary 
	and saves the resulting information into the dictionary

	the checking result will be saved in keys, 
	"profile_integrity".

	Also modifies the dominfo['progress'] will be set to 1 
	for the reruning jobs.
	
	And if the profile has good results,
	"profile_integrity" will be the value of number of iterations.
	If not, profile_integrity will be set to 0.
	'''

	errored_domain_informations = {}
	for domid, dominfo in domain_informations.items() :

		#pass if the profile is already good!
		if dominfo.get( 'profile_integrity' ) :
			continue
		else :
			profile_integrity = check_generated_profile( dominfo )
			dominfo['profile_integrity'] = profile_integrity
			dominfo['progress'] = 1  #prefile rebuilding!
			
			if not profile_integrity :
				errored_domain_informations[domid] = dominfo

	return errored_domain_informations

from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentRecords


def check_profile_alignment( dominfo ) :
	'''
	returns the integrity of profile alignment integrity
	'''
	return dominfo.get('profile_alignment_integrity')

check_profile_alignment_integrity = check_profile_alignment

def _check_profile_alignment( dominfo, filtered_ids ) :

	if dominfo.get('profile_alignment_integrity') :
		return dominfo['profile_alignment_integrity']

	domid = dominfo['uniqueid']
	ext = Settings.get( "alignment_suffix" )
	aln_fn = build_sequence_filename( dominfo['domain_path'] , domid, ext )

	if os.path.exists( aln_fn ) :
		alignments = PairwiseAlignmentRecords()
		alignments.parse(aln_fn )

		for alignment_method in Settings.get( "profile_comparison_methods" ).split() :
			method_name = alignment_method.lower() + "_1"
			if not alignments.count( domid, filtered_ids, method_name ) == len(filtered_ids) :
				return 0
		else :
			return 1
	return 0


def check_profile_alignments( domain_informations ) :
	'''
	Check each domain information dictionary
	and saves the checking informaiton into the dictionary.

	It checks only first iteration alignments are available 
	in the alignments.
	And currently simple flag of good or no good is saved!

	Better resolution "like what method is bad" and 
	'''

	filtered_ids = set()
	for domid, dominfo in domain_informations.items() :
		if check_generated_profile( dominfo ) :
			filtered_ids.add( domid )

	errored_domain_informations = {}
	for domid in filtered_ids :
		dominfo = domain_informations[domid]
		if _check_profile_alignment( dominfo, filtered_ids ) :
			dominfo['profile_alignment_integrity'] = 1
		else :
			dominfo[ 'profile_alignment_integrity' ] = 0
			errored_domain_informations[domid] = dominfo

	return errored_domain_informations

def check_structure_alignments( domain_informations ) :
	'''
	Check each domain information dictionary
	and saves the checking information into the dictionary.
	'''
	#needs to be implemented!
	return []
