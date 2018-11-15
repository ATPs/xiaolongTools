
import os, sys, tempfile, shutil
from subprocess import Popen, PIPE

from evdblib.DBTools import Settings
from evdblib.Utils import parse_sequence_filename
from evdblib.Utils.Parsers.PairwiseAlignmentGenericParser import PairwiseAlignmentRecords

class PairAligner :
	'''
	Pairwise aligner base class.
	'''
	def __init__( self, inputfile1=None, inputfile2=None, cmd=None, parser=None, outputfile=None, cwd=None, verbose=False, 
			identifier1=None, identifier2=None ) :

		self.inputfile1 = inputfile1 #input file1 or query file
		self.inputfile2 = inputfile2 #input file2 or hit file

		self.cmd = cmd
		self.parser = parser
		self.outputfile = outputfile

		self.alignment_raw_result = None
		self.alignment_result = None
		self.verbose = verbose
		self.cwd = cwd

		self.identifier1 = identifier1
		self.identifier2 = identifier2

		if self.outputfile and os.path.exists( self.outputfile ) :
			raise Exception( self.outputfile + " already exists!" )

	def get_command_line( self ) :
		raise NotYetImplementedError( "get_command_line in PairAligner." )

	def run( self ) :
		self.command_line = self.get_command_line()

		if self.outputfile :
			p = Popen( self.command_line, stdout=PIPE, stderr=PIPE, cwd=self.cwd )
			stdout, stderr = p.communicate()

			if self.verbose :
				sys.stderr.write( stderr )
				sys.stdout.write( stdout )
			
			if os.path.exists( self.outputfile ) :
				fp = open( self.outputfile )
				self.alignment_raw_result = fp.read()
				fp.close()

		else :
			p = Popen( self.command_line, stdout=PIPE, stderr=PIPE )
			self.alignment_raw_result = p.communicate()[0]
			self.returncode = p.wait()
			
		if self.parser :
			self.alignment_result = self.parser.parse( self.alignment_raw_result )
		else :
			self.alignment_result = self.alignment_raw_result

	def get_alignment( self ) :
		if self.parser :
			return self.parser.get_alignment()
		raise Exception( "No parser is defined." )

from evdblib.Utils.Parsers.DaliLiteDAT import DaliLiteDAT

class DaliLiteAligner( PairAligner ) :
	def __init__( self, 
		cmd=None,
		inputfile1=None, inputfile2=None, 
		identifier1=None, identifier2=None,
		fakeid1='1domA', fakeid2='2domA',
		outputfile=None, 
		parser=None, verbose=False ) :
		'''
		Alinger class for DaliLite pairwise program.
		
		Currently Parser does not work!
		If verbose is True, output from DaliLite will be
		printed out!
		Use symlink True when the DAT files are on
		fast access directories.
		For remote files, turn off symlink.
		This will make the local copy of .dat files for DaliLite.
		'''
		if cmd == None :
			cmd = Settings.get( 'dalilite' )

		self.temp_dir = tempfile.mkdtemp()

		if inputfile1 == None :
			raise AlignerError( "No input file 1 is given." )

		if inputfile2 == None :
			raise AlignerError( "No input file 2 is given." )

		if outputfile == None :
			self.outputfile = os.path.join( self.temp_dir, os.path.basename(inputfile1).replace('.dat','.dccp') )

		if identifier1 == None :
			identifier1 = parse_sequence_filename(inputfile1)[1]
		if identifier2 == None :
			identifier2 = parse_sequence_filename(inputfile2)[1]

		if os.path.exists( self.outputfile ) :
			raise AlignerError( outputfile, "Outputfile already exists!" )


		self.identifier1 = identifier1
		self.identifier2 = identifier2

		self.fakeid1 = fakeid1
		self.fakeid2 = fakeid2

		#preparing 
		self.dat_dir = os.path.join( self.temp_dir, 'DAT' )
		os.mkdir( self.dat_dir )
		self.dat_file1 = os.path.join( self.dat_dir, self.fakeid1 + '.dat' )
		self.dat_file2 = os.path.join( self.dat_dir, self.fakeid2 + '.dat' )

		dat1 = DaliLiteDAT( inputfile1 )
		dat1.convert_identifier( output=self.dat_file1, output_identifier=self.fakeid1 )
		dat2 = DaliLiteDAT( inputfile2 )
		dat2.convert_identifier( output=self.dat_file2, output_identifier=self.fakeid2 )

		self.temp_output = os.path.join( self.temp_dir, self.fakeid1+".dccp" )

		PairAligner.__init__( self, 
			cmd=cmd, cwd=self.temp_dir,
			inputfile1=inputfile1, 
			inputfile2=inputfile2, 
			identifier1=identifier1,
			identifier2=identifier2,
			outputfile=self.temp_output,
			parser=None,
			verbose=verbose )


	def __del__( self ) :
		if hasattr( self, 'temp_dir') and self.temp_dir :
			shutil.rmtree( self.temp_dir )


	def get_command_line( self ) :
		return [ self.cmd, '-align', self.fakeid1, self.fakeid2 ]


	def run( self ) :
		PairAligner.run( self )
		if os.path.exists( self.temp_output ) :
			if self.temp_output != self.outputfile :
				shutil.copy( self.temp_output, self.outputfile )



class FASTAligner( PairAligner ) :
	def get_command_line( self ) :
		#need to check self.cmd here to avoid writing 
		#another version of __init__
		#probably not a good practice.

		if self.cmd == None :
			self.cmd = Settings.get( "fast" )
			
		return [ self.cmd, self.inputfile1, self.inputfile2 ]

	def run( self ) :
		self.command_line = self.get_command_line()

		p = Popen( self.command_line, stdout=PIPE, stderr=PIPE, cwd=self.cwd )
		stdout, stderr = p.communicate()

		self.alignment_raw_result = stdout
		self.returncode = p.wait()

		if self.verbose :
			print(stdout)
			print(stderr)
			
		if self.parser :
			self.alignment_result = self.parser.parse( self.alignment_raw_result, self.inputfile1, self.inputfile2, )
		else :
			self.alignment_result = self.alignment_raw_result

		if self.outputfile :
			fp = open( self.outputfile, 'w' )
			fp.write( self.alignment_result )
			fp.close()

class TMalignAligner( FASTAligner ) :
	def get_command_line( self ) :
		#need to check self.cmd here to avoid writing 
		#another version of __init__
		#probably not a good practice.

		if self.cmd == None :
			self.cmd = Settings.get( "tmalign" )
			
		return [ self.cmd, self.inputfile2, self.inputfile1 ]

import copy
from evdblib.Utils.Parsers.DaliLite import DaliLite
from evdblib.Utils.Parsers.FAST import FAST
from evdblib.Utils.Parsers.TMalign import TMalign

def run_dalilite( *args ) :
	'''
	Easy interface to run DaliLite.
	returns a PairwiseAlignmentMethodRecord object.

	input parameter is a iterator containing 
	1. inputfile1 (dat file of query)
	2. inputfile2 (dat file of hit)
	3. query id (in case ID in dat1 file should be changed)
	4. hit id (similarly to change hit id in dat2).
	'''
	inputfile1, inputfile2, query_id, hit_id = args 

	#print inputfile1, inputfile2, query_id, hit_id

	dalirunner = DaliLiteAligner( inputfile1=inputfile1, inputfile2=inputfile2 )
	dalirunner.run()
	#the query_id and hit_id should be fake ids.
	#daliparser = DaliLite( dalirunner.outputfile, inputfile1, inputfile2, query_id, hit_id )
	daliparser = DaliLite( dalirunner.outputfile, inputfile1, inputfile2, query_id=dalirunner.fakeid1, hit_id=dalirunner.fakeid2 )
	print('aligned', args, file=sys.stderr)
	daliparser.alignment.id1 = query_id
	daliparser.alignment.id2 = hit_id
	

	return copy.deepcopy( daliparser.alignment )


def run_fast( *args ) :
	'''
	Easy interface to run FAST.
	returns a PairwiseAlignmentMethodRecord object.

	input parameter is a iterator containing 
	1. inputfile1 (CA only PDB file of query)
	2. inputfile2 (CA only PDB file of hit)
	3. query id
	4. hit id
	'''
	inputfile1, inputfile2, query_id, hit_id = args 

	runner = FASTAligner( inputfile1=inputfile1, inputfile2=inputfile2 )
	runner.run()
	parser = FAST( runner.alignment_raw_result, inputfile1, inputfile2, query_id, hit_id )
	print('aligned', args, file=sys.stderr)

	return copy.deepcopy( parser.alignment )


def run_tmalign( *args ) :
	'''
	Easy interface to run TMalign.
	returns a PairwiseAlignmentMethodRecord object.

	input parameter is a iterator containing 
	1. inputfile1 (CA only PDB file of query)
	2. inputfile2 (CA only PDB file of hit)
	3. query id
	4. hit id
	'''
	inputfile1, inputfile2, query_id, hit_id = args 

	runner = TMalignAligner( inputfile1=inputfile1, inputfile2=inputfile2 )
	runner.run()
	parser = TMalign( runner.alignment_raw_result, inputfile1, inputfile2, query_id, hit_id )
	print('aligned', args, file=sys.stderr)

	return copy.deepcopy( parser.alignment )

class DBAligner :
	def __init__( self, cmd=None, inputfile=None, outputfile=None, dbfile=None, listfile=None, pairaligner=None, parser=None, echo=True ) :
		'''
		Database aligner base class.
		Aligns one query against the whole database defined in db file or
		list file.

		outputfile is a saving path of the raw alignment results,
		not for saving final output!
		'''
		self.inputfile = os.path.abspath( inputfile )
		self.dbfile = os.path.abspath( dbfile )

		if listfile :
			self.listfile = os.path.abspath( listfile )
		else :
			self.listfile = listfile

		self.cmd = os.path.abspath( cmd )
		self.pairaligner = pairaligner
		self.parser = parser
		self.outputfile = outputfile
		self.echo = echo

	def get_command_line( self ) :
		raise NotYetImplementedError( "get_command_list in DBAligner." )

	def run( self ) :
		if self.dbfile :
			command_line = self.get_command_line()
			if self.echo :
				print(''.join( command_line ))

			p = Popen( command_line, stdout=PIPE, stderr=PIPE )

			stdout, stderr = p.communicate()
			print(stdout)
			print(stderr, file=sys.stderr)

			if self.parser :
				self.ailgnment_results = self.parser.parse( self.alignment_raw_results )
			
		elif self.listfile :
			listfp = open( listfile )
			inputfile2_list = listfp.read().split()
			listfp.close()

			for inputfile2 in inputfile2_list :
				pairaligner = self.pairaligner( cmd=self.cmd, inputfile1=inputfile, inputfile2=inputfile2, parser=self.parser )
				pairaligner.run()
				self.alignment_results.append( pairaligner.alignment_result )

from evdblib.Utils.ThreadingTools import Pool

class DaliLiteDBAligner :
	def __init__( self, inputfile1_list=None, inputfile2_list=None, outputfile=None, query_id_list=None, hit_id_list=None, echo=True, nthreads=0, pairwise_run=run_dalilite ) :
		'''
		Runs DaliLite and parses the results between 
		inputfile (input1 or query) and dat files in inputfil2_list 
		(aka hit dat files).
		'''
		self.outputfile = outputfile
		self.echo = echo
		self.nthreads = nthreads
		self.inputfile1_list = inputfile1_list 
		self.inputfile2_list = inputfile2_list
		self.query_id_list = query_id_list
		self.hit_id_list = hit_id_list
		self.pairwise_run = pairwise_run

		if not self.inputfile1_list :
			raise AlignerError( "List of inputfile1 should be provided!" )
		
		if not self.inputfile2_list :
			raise AlignerError( "List of inputfile2 should be provided!" )

		if not self.query_id_list :
			self.query_id_list = [None]*len(self.inputfile1_list)

		if not self.hit_id_list :
			self.hit_id_list = [None]*len(self.inputfile2_list)

		self.alignments = PairwiseAlignmentRecords()

	def run( self ) :
		if self.nthreads == 1 or not self.nthreads :
			for inputfile1, inputfile2, query_id, hit_id in zip(self.inputfile1_list, self.inputfile2_list, self.query_id_list, self.hit_id_list) :
				methodrec = self.pairwise_run( inputfile1, inputfile2, query_id, hit_id )
				self.alignments.add( methodrec )
		else :
			
			pool = Pool( self.nthreads )

			apply_result = pool.map_async( 
				self.pairwise_run, 
				list(zip(self.inputfile1_list, 
					self.inputfile2_list, 
					self.query_id_list, 
					self.hit_id_list
				)) 
			)

			method_records = apply_result.get()
			for methodrec in method_records  :
				self.alignments.add( methodrec )

			pool.terminate()
			

class FASTDBAligner (DaliLiteDBAligner) :
	def __init__( self, inputfile1_list=None, inputfile2_list=None, outputfile=None, query_id_list=None, hit_id_list=None, echo=True, nthreads=0, pairwise_run=run_fast ) :
		'''
		Runs FAST and parses the results between 
		inputfile (input1 or query) and hit list files in inputfil2_list 
		(aka hit CA only files).
		'''
		DaliLiteDBAligner.__init__( self, inputfile1_list, inputfile2_list, outputfile, query_id_list, hit_id_list, echo, nthreads, pairwise_run )


class TMalignDBAligner (DaliLiteDBAligner) :
	def __init__( self, inputfile1_list=None, inputfile2_list=None, outputfile=None, query_id_list=None, hit_id_list=None, echo=True, nthreads=0, pairwise_run=run_tmalign ) :
		'''
		Runs TMalign and parses the results between 
		inputfile (input1 or query) and hit list files in inputfil2_list 
		(aka hit CA only files).
		'''
		DaliLiteDBAligner.__init__( self, inputfile1_list, inputfile2_list, outputfile, query_id_list, hit_id_list, echo, nthreads, pairwise_run )





class HHsearchDBAligner( DBAligner ) :
	def __init__( self, cmd=None, inputfile=None, outputfile=None, dbfile=None, cpu=1, db_size=None ) :
		if cmd == None :
			cmd = Settings.get( 'hhsearch_cmd' )
		
		if inputfile == None :
			raise AlingerError( "No input file is given." )

		if outputfile == None :
			self.outputfile_fp = tempfile.NamedTemporaryFile()
			outputfile = self.outputfile_fp.name

		if dbfile == None :
			raise AlignerError( "No DB file is given." )

		self.db_size = db_size

		DBAligner.__init__( self, cmd=cmd, inputfile=inputfile, outputfile=outputfile, dbfile=dbfile )
		#setting number of CPU's can be used.
		self.cpu = cpu

	def get_command_line( self ) :
		command_line = '%s -i %s -o %s -cpu %s -alt 1 -d %s -realign'%(self.cmd, 
			self.inputfile, self.outputfile, self.cpu, self.dbfile)

		if self.db_size :
			command_line += ' -Z %s -z %s -B %s -b %s' %( self.db_size, self.db_size, self.db_size, self.db_size )
		return command_line.split()
			

class COMPASSDBAligner( DBAligner ) :
	def __init__( self, cmd=None, inputfile=None, outputfile=None, dbfile=None ) :
		if cmd == None :
			cmd = Settings.get( 'compass_cmd' )

		if inputfile == None :
			raise AlingerError( "No input file is given." )

		if outputfile == None :
			self.outputfile_fp = tempfile.NamedTemporaryFile()
			outputfile = self.outputfile_fp.name

		if dbfile == None :
			raise AlignerError( "No DB file is given." )

		DBAligner.__init__( self, cmd=cmd, inputfile=inputfile, outputfile=outputfile, dbfile=dbfile )
		#setting number of CPU's can be used.

	def get_command_line( self ) :
		command_line = '%s -i %s -o %s -j %s'%( self.cmd, 
			self.inputfile, self.outputfile, self.dbfile )
		return command_line.split()
	

class AlignerError( Exception ) :
	pass
