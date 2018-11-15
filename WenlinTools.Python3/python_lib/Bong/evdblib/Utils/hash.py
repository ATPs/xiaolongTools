#!/usr/bin/env python
import sys, random

#This script will generate add-hoc hash values for the id when the 

def ad_hoc_numeric_hash( hashable_string, number_of_letters=2 ) :
	''' 
	This function simply use the python hash function to generate initial hash values
	And then shrink the value into desired digits.
	
	Default is 2 digits. 
	
	Note that the number of letters cannot be more than 4 digits, due to the design issues
	the maximal hash value cannot be more than ~32000 and this in turn limit the hash digit size
	be in between 2 to 4.
	'''
	return ('%0' + str(number_of_letters) + 'd') % ( hashable_string.__hash__()%( 10**int(number_of_letters) ) )

def random_numbers( number_of_letters=2 ) :
	return ( '%0' + str(number_of_letters) + 'd') % int((random.random() * 10**int(number_of_letters)))

import string
def pdb_style_hash( hashable_string, position=1 ) :
	'''
	This function simply returns a hash value of length 2.
	To make each directory containing the files not too big, PDB had simple idea of
	taking two middle letters from its ID.
	It mimick this PDB style hash value and returns two letters from the position in lower cases.

	If the ID is not long to get the hash values from the given position,
	randomly selected two alphanumeric characters will be returns
	'''

	s = str(hashable_string)
	hash = s[position:position+2]

	if hash and len(hash) == 2 :
		return hash.lower()
	else :
		selections = string.lowercase + string.digits
		return random.choice(selections) + random.choice( selections )
	

if __name__ == '__main__' :
	id_fn = '/home/bhk/projects/hua_update/scop40.1.75'
	fp = open( id_fn )
	id_list = fp.read().split()
	
	number_of_letters = 2
	
	hash_table = {}
	for id in id_list :
		#hash = ad_hoc_numeric_hash( id, 2 )
		hash = pdb_style_hash( id, 2 )
		print(id, hash)
