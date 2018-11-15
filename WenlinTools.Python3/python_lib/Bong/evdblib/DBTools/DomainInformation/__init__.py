'''
This subpackage contains modules to manage information about proteins 
or domains in the database.
'''
#debug = 1
verbose=1

import os, shelve
from evdblib.DBTools import Settings

domain_info_db = Settings.get('domain_info_db')
classification_info_db = Settings.get('classification_info_db')
classification_levels = int( Settings.get('classification_levels' ) )

#########################################
#Using local database
#########################################
if Settings.get( 'use_local_db' ) :
	cwd = os.getcwd()
	domain_info_db = os.path.join( cwd, os.path.basename( domain_info_db ) )
	classification_info_db = os.path.join( cwd, os.path.basename( classification_info_db ) )

##################################k
#temporary blocking!!!
##################################k
#class DomInfoDB :
	#def __init__( self, domain_info_db ) :
	#	self.dominf = shelve.open( domain_info_db )

	'''
	def __del__( self ) :
		if hasattr( self, 'dominfo' ) and self.dominf :
			self.dominf.close()
	'''
	
#class ClaInfoDB :
	#def __init__( self, classification_info_db ) :
	#	self.clainf = shelve.open( classification_info_db )

	'''
	def __del__( self ) :
		if hasattr( self, 'clainf' ) and self.clainf :
			self.clainf.close()
	'''
	

dominfodb = None #DomInfoDB( domain_info_db ) #shelve.open( domain_info_db )
#domcla = shelve.open( domain_classification_db )
#classification is included into the domain information data

#contains sets of classifications
#csets, psets
clainfodb = None #ClaInfoDB( classification_info_db ) #shelve.open( classification_info_db )

def generate_intermediate_dbs( domain_informations ) :
	#Settings.get( 'intermediate_result_dir' )
	intermediate_domain_info_db = Settings.get( 'intermediate_domain_info_db' )
	intermediate_classification_info_db = Settings.get( 'intermediate_classification_info_db' )

	#Need to think if the intermediate db is needed!!


def get_record( uniqid, domaininfo_db=None ) :
	'''
	returns the records of the unique ID.

	domaininfo_db can be used for extracting info from 
	non-default db.
	If no record was found for the unique ID, 
	this function will return None.
	'''
	if domaininfo_db == None :
		global dominfodb
		if not dominfodb :
			dominfodb = shelve.open( domain_info_db )
		domaininfo_db = dominfodb #dominfodb.dominf

	inf = domaininfo_db.get( uniqid ) #dictionary of keys "uniqueid, name, comment, exptype, resolution, range"
	#cla = domcla.get( uniqid ) #tuple of classification words
	
	return inf

def get_all_records( domaininfo_db = None ) :
	'''
	returns list of all records in the database
	'''
	if  domaininfo_db == None :
		global dominfodb
		if not dominfodb :
			dominfodb = shelve.open( domain_info_db )
		domaininfo_db = dominfodb#.dominf
	else :
		pass

	#temporarily return []
	#need to fix the nast None-Type calling error 
	#at the end of program run.
	#that probably came from two DB files..

	return list(domaininfo_db.values())

def get_domaininfo_db( domaininfo_db = None ) :
	'''
	returns list of all records in the database
	'''
	if  domaininfo_db == None :
		global dominfodb
		if not dominfodb :
			dominfodb = shelve.open( domain_info_db )
		domaininfo_db = dominfodb#.dominf
	else :
		pass

	#temporarily return []
	#need to fix the nast None-Type calling error 
	#at the end of program run.
	#that probably came from two DB files..

	return domaininfo_db


def is_domain_information_unique( domain_informations ) :
	'''
	try to check if the uniqueness of the uniqid.
	If given domain_information dictionary (new information adding) 
	contains keys that is overlapping with the current uniqueID,
	this will return False
	'''
	for uniqid in domain_informations.keys() :
		if get_record(uniqid) :
			#overlapping uniqid key found!
			return False

	return True

def add_records( dominfos ) :
	for domid in dominfos :
		if not set_record( domid, dominfos[domid] ) :
			raise Exception( "Error in setting record" )

def del_records( domids, strict=False ) :
	for domid in domids :
		if not del_record( domid ) :
			if strict :
				raise Exception( "Error in deleting record" )
			else :
				print(domid, "deletion failed!")

def set_record( inf, dominfo, domaininfo_db=None, classinfo_db=None) :
	'''
	returns True if the domain information setting has done successfully.
	Otherwise, returns False.
	'''
	#check if the non-default db's usage
	if domaininfo_db :
		dominf = domaininfo_db
	else :
		global dominfodb
		if not dominfodb :
			dominfodb = shelve.open( domain_info_db )
			
		dominf = dominfodb#.dominf

	if classinfo_db :
		clainf = classinfo_db
	else :
		global clainfodb
		if not clainfodb :
			clainfodb = shelve.open( classification_info_db )
		clainf = clainfodb#.clainf

	#check the uniqueness of the path first!!
	if set_classification_path( dominfo, classinfo_db=clainf ) :
		pass
	else :
		return False
	
	dominf[inf] = dominfo

	return True

def del_record( inf, domaininfo_db=None, classinfo_db=None ) :
	'''
	returns True if the domain information has the domid and has successfully removed.
	'''
	#check if the non-default db's usage
	if domaininfo_db :
		dominf = domaininfo_db
	else :
		global dominfodb
		if not dominfodb :
			dominfodb = shelve.open( domain_info_db )
			
		dominf = dominfodb#.dominf

	if classinfo_db :
		clainf = classinfo_db
	else :
		global clainfodb
		if not clainfodb :
			clainfodb = shelve.open( classification_info_db )
		clainf = clainfodb#.clainf

	if inf in dominf :
		del dominf[inf]
		return True
	else :
		return False



def finalize_dbs(*args) :
	'''
	Closes all shelve db's so that the update will be saved into the disk.

	Note the the shelve db's should be closed to properly save the changes in the db.

	Optionally the list of shelve db's given in the argument,
	can be all properly closed!
	'''

	if args :
		for db in args :
			db.close()
	else :
		global dominfodb, clainfodb

		if dominfodb :
			dominfodb.close()
			dominfodb = None
		if clainfodb :
			clainfodb.close()
			clainfodb = None


from evdblib.Utils import string2pathname
def is_classification_path_unique( dominfo, classinfo_db=None ) :
	'''
	returns True when the domain classification defined in the dominfo 
	is also converted to a unique pathname. 
	Otherwise, returns False.

	This check is necessary since the classification names will be used and
	directory names but all special characters should be removed.

	classinfo_db can be set to test with the custom made classification information 
	db.
	'''
	clainf_c2p_key_temp = '%dc2p' 
	clainf_p2c_key_temp = '%dp2c'

	uniqueid = dominfo['uniqueid']

	#if classinfo_db is given 
	#the check is done on the custom db.
	if classinfo_db == None :
		global clainfodb
		if not clainfodb :
			clainfodb = shelve.open( classification_info_db )
		classinfo_db = clainfodb

	for level in range( 1, 1+classification_levels ) :
		c2pkey = clainf_c2p_key_temp % level
		p2ckey = clainf_p2c_key_temp % level

		#classification is tuple

		#classification is tuple
		claval = dominfo['classification'][level-1] 
		clapath = string2pathname( claval )

		#compare to the previously saved paths
		if c2pkey in classinfo_db :
			c2p = classinfo_db[ c2pkey ] 
		else :
			c2p = {}

		if p2ckey in classinfo_db :
			p2c = classinfo_db[ p2ckey ]
		else :
			p2c = {}

		#two conditions that are OK
		if claval not in c2p and clapath not in p2c :
			continue
		elif claval in c2p and clapath == c2p[claval] :
			continue
		else :
			return False
			
	return True #at the end of check.. it is good.
		
def set_classification_path( dominfo, classinfo_db=None ) :
	'''
	updates the domain information and classification record of
	the uniqid.

	classinfo_db can be set to add to the custom made classification information 
	db. Note that this function will not set the path and classification mapping
	when the pathname and classification do not have one to one correspondence.

	If the update happens successfully, this function will return True,
	otherwise, it will return False.
	'''
	clainf_c2p_key_temp = '%dc2p' 
	clainf_p2c_key_temp = '%dp2c'

	uniqueid = dominfo['uniqueid']

	#if classinfo_db is given 
	#the check is done on the custom db.
	if classinfo_db == None :
		pass
	else :
		clainf = classinfo_db

	if is_classification_path_unique( dominfo, classinfo_db=classinfo_db ) :
		pass
	else :
		return False

	for level in range( 1, 1+classification_levels ) :
		c2pkey = clainf_c2p_key_temp % level
		p2ckey = clainf_p2c_key_temp % level

		#classification is tuple
		claval = dominfo['classification'][level-1] 
		clapath = string2pathname( claval )

		#compare to the previously saved paths
		if c2pkey in clainf :
			c2p = clainf[ c2pkey ] 
		else :
			c2p = {}

		if p2ckey in clainf :
			p2c = clainf[ p2ckey ]
		else :
			p2c = {}

		c2p[claval] = clapath
		p2c[clapath] = claval

		clainf[ c2pkey ] = c2p
		clainf[ p2ckey ] = p2c

	return True

from evdblib.Utils.hash import pdb_style_hash, random_numbers
def read_information_files( domain_information_file, domain_classification_file ) :
	'''
	reads in the domain_information_file and domain_classification_file and format the db record.
	'''
	domain_informations = {}
	number_of_info_columns = 8
	if os.path.exists( domain_information_file ) :
		fp = open( domain_information_file )
		for l in fp.readlines() :
			line = l.split()
			for i in range( len(line), number_of_info_columns  ) :
				line.append( '' )

			uniqueid, proteinfile, range, verserion, name, comment, exptype, resolution = line
			if uniqueid not in domain_informations :
				domain_informations[ uniqueid ] = {"uniqueid":uniqueid, "proteinfile":proteinfile, "range":range, "name":name, "usercomment":comment, "exptype":exptype, "resolution":resolution}
			else :
				raise DomainInformationReadingError( "Unique ID is defined redundantly. %s"%uniqueid ) 
	else :
		raise DomainInformationReadingError( "No Domain information file found! %s "% domain_information_file ) 

	if verbose :
		print(domain_information_file, os.path.exists( domain_information_file ))
		print(domain_classification_file, os.path.exists( domain_classification_file )) 
		

	if domain_classification_file :
		if not os.path.exists( domain_classification_file ) :
			raise DomainInformationReadingError( "No domain classification file found! %s"%domain_classificaiton_file )
		fp = open( domain_classification_file )
		for l in fp :
			line = l.strip().split('\t')

			if verbose :
				print(line)

			uniqid = line[0]
			if uniqid in domain_informations :
				if len(line) == classification_levels + 1:
					domain_informations[uniqid]['classification'] = line[1:]

					if verbose :
						print(domain_informations[uniqid])

				else :
					raise DomainInformationReadingError( "Classification column number is not right!"+ str(line) )
				
	else :
		for uniqid, dominfo in domain_informations.items() :
			if data_type == 'structure' and classification_levels == 1 and (len(uniqid)==4 or len(uniqid)==5):
				dominfo['classification'] = (pdb_style_hash(uniqid),)
			else :
				dominfo['classification'] = tuple([ random_numbers(len(uniqid)) for i in range(classification_levels) ])

	#final check for the classification assignment.
	for dominfo in domain_informations.values():
		if not "classification" in dominfo :
			raise DomainInformationReadingError( "classification information is not assinged. " + dominfo['uniqueid'] )

	return domain_informations

class DomainInformationReadingError( Exception ) :
	pass
