'''
This subpackage contains classes and functions for parsing (PSI)BLAST results.

The main module BLAST contains MSA objects (whole alignments in a BLAST 
or alignments in a iteration of a PSI-BLAST).

MSA in turn contains hits.
Then a hit contains alingments of the query and the hit.

Note that the driver scripts that runs BLAST or PSI-BLAST
are seperated into seperate class BLASTRunner or its derived classes.
'''
############################
#verbose output for debugging
verbose = 0
############################

###########################
# Parsing constants
###########################
psiblast_end_of_result_marker = "  Database:"
psiblast_end_of_result_iteration = "Searching..................................................done"
psiblast_start = 'Results from round'
psiblast_start_m6 = 'QUERY '
hit_record_start = '>'
hsp_start = ' Score ='
###########################
import sys, os, io
from ..FASTA import FASTA
from evdblib.Utils import is_eof

#probably function parse should be defined here to parse a BLAST
#psiblast result.

class NoMoreIteration (Exception) :
	'''
	When no more Iteration is avaiable in the stream.
	'''
	pass

class NoMoreHit( Exception ) :
	'''
	When no more hits are available in the current iteration.
	'''
	pass

class NoMoreAlignment( Exception ) :
	'''
	When no more alignments are available in the current hit.
	'''
	pass

class IterationFormatError (Exception) :
	'''
	When the BLAST or PSIBLAST iteration has a problem!

	This Error frequently need to be caught because
	the PSI-BLAST ouptut might has some problem but
	the subsequent iteration it might be OK. 
	'''
	pass
	
class FormatError (Exception) :
	'''
	When the alignment has a problem.
	'''
	pass

class NullBLASTResultError (Exception) :
	'''
	No BLAST result is in File
	'''
	pass

from .MSA import MSA
class BLAST :
	'''
	The main class in the BLAST package.

	Parses a BLAST or PSI-BLAST result
	and provide a list of functions to manipulate them.

	Note:
	This class currently parses only a text output 
	of a legacy BLAST or PSI-BLAST (output of -m default option).
	'''
	def __init__( self, query_fn='',  blast_result_fn='', query=None, blast_result_fp=None, null_result_error=False ) :
		#each PSI-BLAST iteration will be save as each MSA.
		self.MSAs = []	 

		self.query_fn = query_fn
		self.query = query #in case for FASTA object is passed.
		self.null_result_error = null_result_error

		#if query file is given,
		#use informaiton in the file.
		if self.query_fn :
			self.query = FASTA( query_fn )

		self.fn = blast_result_fn
		self.fp = blast_result_fp
		self.fp_position = 0

		#if blast result file is given,
		#use the blast result file.
		if blast_result_fn :
			self.fp = open( blast_result_fn )
		elif self.fp :
			#record the initial point for later usage.
			self.fp_position = self.fp.tell()

		#no input stream available
		if not self.fp :
			if verbose :
				print("MSA: No input stream is available!")
			return

		self.parse()

	def __len__(self ) :
		return self.MSAs.__len__()

	def __getitem__( self, i ) :
		return self.MSAs[i]

	def __iter__( self ) :
		return self.MSAs.__iter__()

	def parse(self) :
		'''	
		Main parse funtion of the BLAST output.
		'''
		line = self.fp.readline()
		while line and not line.startswith( hit_record_start ) :
			line = self.fp.readline()
		else :
			if not line :
				if self.null_result_error :
					raise NullBLASTResultError( 'No BLAST result found!' )
				else :
					#build Null MSA object!
					msa = MSA(self.query, self.fp )
					self.MSAs.append(msa)
					return
					
			#rewinde to the point of line start
			self.fp.seek( -len(line), 1 )

		
		while not is_eof(self.fp) :
			try :
				#note that it is very crucial that 
				#the statement might stop within.. :)
				msa = MSA(self.query, self.fp ) 
				msa.parse()
			except NoMoreIteration :
				if msa :
					self.MSAs.append(msa)
				break
			except FormatError as e :
				#This Error occurs when the output of a certain
				#iteration in the middle has some problem!
				print("Error in parsing an iteration", e, file=sys.stderr)
				pass

			if len(msa) :
				self.MSAs.append( msa )


