'''
PDB package can parse a PDB file
and perform various checking on them.

The package is developed with an emphasis 
on extracting protein sequence from PDB file.
'''
import io

#####################################
#amino acid name conversion tables
amino_acids = {'GLY':'G', 'PRO':'P', 'ALA':'A', 
	'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 
	'CYS':'C', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 
	'HIS':'H', 'LYS':'K', 'ARG':'R', 'GLN':'Q', 
	'ASN':'N', 'GLU':'E', 'ASP':'D', 'SER':'S', 
	'THR':'T'}

amino_acids_rev = { 'G':'GLY', 'P':'PRO', 'A':'ALA', 
	'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET', 
	'C':'CYS', 'F':'PHE', 'Y':'TYR', 'W':'TRP', 
	'H':'HIS', 'K':'LYS', 'R':'ARG', 'Q':'GLN', 
	'N':'ASN', 'E':'GLU', 'D':'ASP', 'S':'SER', 
	'T':'THR', 'X': 'UNK' }
#####################################
#verbose standard output setting
verbose = 0
#####################################


def parse( fn='', fp=None ) :
	'''
	returns a PDB class instance from PDB file given
	in the filename or file pointer.
	'''
	from .PDB import PDB
	return PDB( fn=fn, fp=fp )

##################################
#Package-wise utility function
##################################
def build_residue_id( chainid, resnum, insertion_code ) :
	'''
	returns a unique identifier of a residue in a PDB file.
	
	This function is designed to 
	return a meaningful unique identifier of a residue
	that can be a hash key.
	'''

	return (chainid, int(resnum), insertion_code)

def build_one_bigger_residue_id( residue_id ) :
	'''
	returns a unique identifier of a residue in a PDB file.

	This function is designed to be used internal purpose.
	'''
	chainid, resnum, insertion_code = residue_id
	return ( chainid, resnum+1, insertion_code )

def resnames2sequence( resnames ) :
	'''
	converts a list of three letter residue names
	as appeared in PDB file,
	into a string of one letter amino aicd codes.

	Note that the non-standard amino acid names will 
	be converted into X
	'''
	fp = io.StringIO()
	for resname in resnames :
		if resname in amino_acids :
			fp.write( amino_acids[resname] )
		else :  
			fp.write( 'X' )
	seq = fp.getvalue()
	fp.close()
	return seq


######################################
#defining Errors in Parsing PDB.
#in main PDB module level
######################################
from .. import ParseError 
class FormatError( ParseError ) :
	'''
	This error means that the PDB format is not curable.
	'''
	pass #shallow wrapper for the Error cases.


class NullError( ParseError ) :
	'''
	This error means that the PDB file or pointer does not 
	have any content in it.
	'''
	pass


class DetailedFormatError( ParseError ) :
	'''
	This error means that the PDB format has some problem
	but maybe curable using fallback PDB parser 
	which do not use advanced features.
	'''
	pass


class EOFErorr( ParseError ) :
	'''  
        This error occurs when the file pointer is trying to read
	at the end of the file.
        '''
	pass


class PDBInfoRecordFormatError( ParseError ) :
	'''
	This error occurs when parsing for the class PDBInfo has a problem.

	Note that the exact record where the error occured will be delivered through
	the message routine.
	'''
	pass


class ModelFormatError( ParseError ) :
	'''
	This error occurs when the Model has problem in Format.
	Or problem in attempting to parse the model part initially.

	???
	Not sure if this is correct place to put this exception
	or this is necessary.. let's see.
	'''
	pass

class MissingResidueFormatError( ParseError ) :
	'''
	This error occurs when the SequenceInfo failed to parse missing residues.
	'''
	pass

class AtomFormatError( ParseError ) :
        '''
        This error means that the Parsing in ATOM record
        of PDB file has problem.

        Most likely due to the formatting issues in the line.
        '''
        pass

class NoAtomLineErorr( ParseError ) :
        '''
        This error occurs when the line to be parsed by Atom class
        is not ATOM or HETATM or TER or ANISOU line.
        '''
        pass

class SEQRESMappingError( ParseError ) :
	'''
	This error occurs when the Mapping of SEQRES record
	to the ATOM records sequence.
	'''
	pass

###############################################
###############################################
# Internet Utility functions
###############################################
###############################################
# Downloading PDB files using pdb id
# from PDB website 
# and Download NCBI genebank file
import urllib.request, urllib.parse, urllib.error, gzip
def download_pdb_file( pdbid, save_dir = '', pdb_fn = '', check_local=1 ) :
	'''gets pdb file from the PDB website and return pdb_fn.
	if pdb_fn is not given, use the pdb_fn <pdbid>.pdb.
	Caution, if the file is already available, the file will be overwritten.
	
	If check_local is True, check the current working directory for the pdb file.
	'''
	if check_local and os.path.exists( pdbid + '.pdb' ) :
		return pdbid + '.pdb'
	
	#default pdb_fn is <pdbid>.pdb
	if not pdb_fn :
		if save_dir :
			pdb_fn = save_dir + '/' + pdbid + '.pdb'
		else :
			pdb_fn = pdbid + '.pdb'
	pdbid = pdbid.lower()
	
	#getting file from the PDB website
	urltemplate = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' 
	pdb_url = urltemplate % pdbid
	'''
	#proxy does not needed for urllib.urlretrieve function.. ??Strange??
	proxy_string = os.getenv( 'http_proxy' )
	if proxy_string == None :
		proxies = None
	else :
	#checking for proxy setting
		proxies = {'http': proxy_string }
	print proxies
	'''
	try : 
		downloaded_filename, http_message = urllib.request.urlretrieve( pdb_url)#, proxies=proxies )
	except IOError:
		print("WARNING! A problem occurred while downloading..", pdbid, chainid, file=sys.stderr)
		return ''
	
	try:
		gfp = gzip.GzipFile( downloaded_filename )
		#print gfp.read()
		pdb_content = gfp.read()
		gfp.close()
		pdb_fp = open( pdb_fn, 'w' )
		pdb_fp.write( pdb_content )
		pdb_fp.close()
	except :
		print("WARNING! A problem occured while uncompressing the downloaded file.", pdbid, chainid, file=sys.stderr)
		return ''
	return pdb_fn

def download_genbank_pdbchain_record( pdbid, chainid ) :
	''' returns the downloaded text as a string.
	Note that the chainid can be either upper or lower!!
	Interestingly, the NCBI genbank specify the chainid for lowercase chainids by 
	putting two characters.
	'''
	#########
	#building accession code!
	pdbid = pdbid.upper()
	if chainid.islower() :
		chainid = chainid.upper() + chainid.upper() #double the chainid
	accession = pdbid + '_' + chainid
	#retrieving genbank file using Entrez
	efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=%s&rettype=gb&email=kim@chop.swmed.edu'
	
	try :
		genbank_fn, http_message = urllib.request.urlretrieve( efetch % accession )
		fp = open(genbank_fn )
		return fp.read()
	except :
		print("WARNING! A problem occurred while downloading the file.", pdbid, chainid, file=sys.stderr)
		return ''


