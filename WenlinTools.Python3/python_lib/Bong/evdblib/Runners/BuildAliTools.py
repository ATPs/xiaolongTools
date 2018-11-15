'''
This module contains tools for building profiles. 
'''
import os, sys, tempfile, shutil, glob, pickle
from subprocess import Popen, PIPE

from evdblib.DBTools import Settings 
from evdblib.Runners import Runner
from evdblib.Runners.BLASTRunners import PSIBLASTRunner
from evdblib.Utils import parse_sequence_filename, parse_profile_filename, build_profile_filename
from evdblib.Utils.Parsers.Range import SequenceRange

verbose = 0

class BuildAli :
	'''
	Build alignments using buildali.pl written by Soding group in Germany.

	Note that this method is not able to run background.
	Also delecate control of temporary directory is not required since
	buildali program uses temporary directory by itself. 
	'''
	def __init__( self, cmd=None, dominfo=None, 
			save_dir=None, sequence=None, range=None, max_iterations=None, selected_iterations=None, 
			max_rerun=3, stdout=PIPE, stderr=PIPE, timeout=0 ) :
		if cmd == None :
			cmd = Settings.get( "buildali" )

		self.cmd = cmd
	
		#getting information out of dominfo
		if dominfo :
			self.save_dir = dominfo['domain_path' ]
			self.sequence = dominfo['profile_sequence_file']
			self.range = dominfo['profile_sequence_range']
		elif save_dir and sequence :
			self.save_dir = save_dir
			if not os.path.exists( self.save_dir ) :
				os.makedirs( self.save_dir )
			self.sequence = sequence
		else :
			raise ProfileBuildingError( "DomainInformation should be given." )

		if not os.path.exists( self.sequence ) :
			raise ProfileBuildingError( "Sequence file %s cannot be found!"%self.sequence )

		#getting profile iteration information from setitng
		if max_iterations == None :
			max_iterations = int(Settings.get( 'max_iterations' ))
		self.max_iterations = max_iterations
		
		if selected_iterations == None:
			selected_iterations = [ int(i) for i in Settings.get( 'selected_iterations' ).split() ]
		self.selected_iterations  = selected_iterations
			

		if not self.save_dir :
			raise ProfileBuildingError( "Saving directory should be set." )
	
		self.range=range
		if self.range :
			self.prepare_fasta(range=self.range)
		
		self.max_rerun = max_rerun

		#getting command line arguments and setting timedrunner
		self.command_line = self.get_command_line()
		#Runner.TimedRunner.__init__( self, self.command_line, **kwargs )
		self.timeout = timeout
		self.stdout = stdout
		self.stderr = stderr

	def _initiate_runner( self ) :
		'''
		returns runner class.
		'''
		return Runner.TimedRunner( 
				self.command_line, 
				timeout=self.timeout, 
				stdin=None, 
				stdout=self.stdout, 
				stderr=self.stderr 
			)

	def prepare_fasta( self, range=None ) :
		'''
		build temporary fasta input file by extracting regions defined in the range.
		'''

		if range == None :
			range = self.range
		
		if not range :
			return 

		sequencerange = SequenceRange()

		from evdblib.Utils.Parsers.FASTA import FASTA
		fasta = FASTA(self.sequence)
		extracted_fasta_string = fasta.extract_fasta_string( range )

		fp = open( os.path.join(self.save_dir, os.path.basename(self.sequence)), 'w' )
		#fp = tempfile.NamedTemporaryFile( suffix='.fa')
		fp.write( extracted_fasta_string )
		fp.flush()
		self.temporary_sequence_fp = fp
		self.sequence = fp.name


	def run( self, echo=True ) :

		if echo :
			print(self.command_line)

		self.runner = self._initiate_runner()
		self.runner.run()
		stdout, stderr = self.runner.communicate()
		returncode = self.runner.poll()

		print(stdout)
		print(stderr, file=sys.stderr)

		#rerunning until get it.
		rerun_counter = 0
		while returncode != 0:

			if rerun_counter > self.max_rerun :
				break

			self.runner = self._initiate_runner()
			self.runner.run()
			stdout, stderror = self.runner.communicate()
			returncode = self.runner.poll()
			rerun_counter += 1 

		self.set_result_files()

		return returncode


	def get_iteration_file( self, i ) :
		'''
		returns filename that contains the given iteration.
		If the filename with iteration is found, it will return None.
		'''
		for a3m_file in self.a3m_files :
			iteration = parse_profile_filename( a3m_file )[-2]
			if int(i) == int(iteration) :
				return a3m_file
		return None
		

	def set_result_files( self ) :
		'''
		returns list of a3m files in the save_dir.
		'''
		basename = parse_sequence_filename( self.sequence )[1]
		if not basename :
			return
		self.a3m_files = glob.glob( os.path.join( self.save_dir, basename + '*.a3m' ) )
		

	def get_command_line( self ) :
		return '%s -n %d -cn %s %s '%(self.cmd, int(self.max_iterations ), self.sequence, self.save_dir  )


class SimpleRunner :
	'''
	Abstract class for running prorams
	where complex controls are not needed!
	'''
	def __init__( self, command_line ) :
		self.command_line

	def run( self ) :
		p = Popen( self.command_line, stderr=PIPE, stdout=PIPE )
		stdout, stderr = p.communicate()
		self.returncode = p.poll()
		print(stdout)
		print(stderr, file=sys.stderr)


class ReformatMSA (SimpleRunner):
	'''
	Reformats MSA majorily in "A3M" to "PSI".
	This is a very shallow wrapper of Soding's Reformat.pl script.

	WARNING: This is a heftly written class that does not check
	rigorous error checkings.
	'''
	def __init__ ( self, cmd=None, input_file=None, output_file=None, input_type=None, output_type=None, remove_query_gap=True ) :

		if cmd == None :
			cmd = Settings.get( "reformat" )

		self.input_file = input_file
		self.output_file = output_file
		self.input_type = input_type
		self.output_type = output_type

		self.cmd = cmd

		self.command_line = self.get_command_line()
		self.remove_query_gap = remove_query_gap

		
	def get_command_line( self ) :

		if self.input_type == None :
			self.input_type = 'a3m'

		if self.output_type == None :
			self.output_type = 'psi'
			
		if self.output_file == None :
			dir, base_filename, iteration, ext = parse_profile_filename( self.input_file )
			self.output_file = build_profile_filename(dir, base_filename, iteration, '.'+self.output_type )
		command_line = [self.cmd, self.input_type, self.output_type, self.input_file, self.output_file]
		return command_line
			

	def run( self ) :
		'''
		Refomat from a3m to psi and then remove query gaps.
		'''
		p = Popen( self.command_line, stdout=PIPE, stderr=PIPE )
		stdout, stderr = p.communicate()
		print(stdout)
		print(stderr, file=sys.stderr)
		self.returncode = p.poll()

		if self.remove_query_gap :
			fp = tempfile.NamedTemporaryFile()
			temp_name = fp.name
			remove_query_gaps( self.output_file, temp_name )
			shutil.copy( temp_name, self.output_file )
			fp.close()


class COMPASSDBBuilder (SimpleRunner) :
	'''
	Class for building numerical COMPASS database.
	'''
	def __init__ ( self, cmd=None, input_file=None, output_file=None ) :

		if cmd == None :
			cmd = Settings.get( "mk_compass_db" )

		if input_file == None :
			raise TypeError ( "Input_file should be given." )

		self.input_file = input_file
		self.output_file = output_file
		#self.input_type = input_type
		#self.output_type = output_type
		self.cmd = cmd

		self.command_line = self.get_command_line()
		
	def get_command_line( self  ) :
		if self.output_file == None :
			dir, base_filename, iteration, ext = parse_profile_filename( self.input_file )

			if verbose :
				print("COMPASS Builder file analysis")
				print('input_file', self.input_file)
				print("dir:", dir)
				print("base_filename:", base_filename)
				print("ieration:", iteration)
				print("ext:", ext)
				
			self.output_file = build_profile_filename(dir, base_filename, iteration, '.cnp' )

		self.temporary_output_fp = tempfile.NamedTemporaryFile()
		self.temporary_output_file = self.temporary_output_fp.name
		list_fn = self.prepare_list_file()
		command_line = [self.cmd, '-i', list_fn, '-o', self.temporary_output_file]
		return command_line
			

	def prepare_list_file( self ) :
		'''
		returns temporary list filename.
		'''
		tempfp = tempfile.NamedTemporaryFile()
		tempfp.write( self.input_file )
		tempfp.flush()
		self.list_fp = tempfp

		return tempfp.name


	def run( self, echo=True ) :
		if echo :
			print(' '.join( self.command_line ))

		p = Popen( self.command_line, stdout=PIPE, stderr=PIPE )
		stdout, stderr = p.communicate()
		print(stdout)
		print(stderr, file=sys.stderr)

		try :
			shutil.copy( self.temporary_output_file, self.output_file )
			shutil.copy( self.temporary_output_file+'.len', self.output_file+'.len' )
		except :
			print("COMPASS DB file generation failed!", self.output, file=sys.stderr)
		
		self.returncode = p.poll()

		self.list_fp.close()
		self.temporary_output_fp.close()

		return self.returncode


class HHMBuilder (SimpleRunner) :
	'''
	Class for building HHsearch Hidden Markov Model.
	'''
	def __init__ ( self, cmd=None, input_file=None, output_file=None, calibration_db=None, calibrate=True, calibration_cmd=None ) :

		if cmd == None :
			cmd = Settings.get( "hhmake" )

		if input_file == None :
			raise TypeError ( "Input_file should be given." )

		if calibration_cmd == None :
			calibration_cmd = Settings.get( "hhsearch_cmd" )

		if calibration_db == None :
			calibration_db = Settings.get( "hhm_cal_db" )

		if output_file == None :
			dir, base_filename, iteration, ext = parse_profile_filename( input_file )
			output_file = build_profile_filename(dir, base_filename, iteration, '.hhm' )

		self.cmd = cmd
		self.calibration_cmd = calibration_cmd
		self.calibration_db = calibration_db
	
		self.calibrate = calibrate

		self.input_file = input_file
		self.output_file = output_file

		#the following part is added due to hhsearch bug
		#of cannot handle long input file name handling!
		self.tmpinput = tempfile.NamedTemporaryFile()
		self.tmpinputname = self.tmpinput.name
		try :
			shutil.copy( self.input_file, self.tmpinputname )
		except IOError :
			self.tmpinputname = self.input_file

		self.tmpoutput = tempfile.NamedTemporaryFile()
		self.tmpoutputname = self.tmpoutput.name
		#need to be copied after the execution!

		self.command_lines = self.get_command_lines()
	
	def get_command_lines( self ) :

		################################
		#to fix long filename bug
		################################
		#command_line = [self.cmd, '-i', self.input_file, '-o', self.output_file]
		command_line = [self.cmd, '-i', self.tmpinputname, '-o', self.tmpoutputname]
		command_lines = [command_line ]

		if self.calibrate :
			self.tmpfp = tempfile.NamedTemporaryFile()
			#command_lines.append(  [ self.calibration_cmd, '-cal', '-d', self.calibration_db, '-i', self.output_file, '-o', self.tmpfp.name ] )
			command_lines.append(  [ self.calibration_cmd, '-cal', '-d', self.calibration_db, '-i', self.tmpoutputname, '-o', self.tmpfp.name ] )

		return command_lines


	def run( self, echo=True ) :

		for command_line in self.command_lines :
			if echo :
				print(' '.join( command_line ))
			p = Popen( command_line, stdout=PIPE, stderr=PIPE )
			stdout, stderr = p.communicate()

			print(stdout)
			print(stderr, file=sys.stderr)

		try :
			shutil.copy( self.tmpoutputname, self.output_file )
		except :
			print("Failed to BuildHHM!", file=sys.stderr)

		self.tmpfp.close()
		self.tmpoutput.close()

		return p.poll()


import evdblib.ProfileScoreTools as profile_score_module
class NumericalProfileBuilder( SimpleRunner ) :
	'''
	Builds python numerical profile file.
	'''
	def __init__( self, input_file=None, output_file=None ) :
		if input_file == None :
			raise TypeError( "No Input file is given." )
		
		if output_file == None :
			dir, base_filename, iteration, ext = parse_profile_filename( input_file ) 
			output_file = build_profile_filename( dir, base_filename, iteration, ".pnp" )
		
		self.input_file = input_file
		self.output_file = output_file


	def run( self ) :
		msa_seqs = [ l.split()[-1] for l in open( self.input_file ).readlines() ]

                neffs, Qs, pssm = profile_score_module.get_neffs_Qs_and_pssm( msa_seqs )

                #save_fp = open( aln[:-4] + '.pnp', 'w' )
		save_fp = open( self.output_file, 'w' )
                pickle.dump( neffs, save_fp, -1 )
                pickle.dump( Qs, save_fp, -1 )
                pickle.dump( pssm, save_fp, -1 )
                save_fp.close()

class PSIPREDRunnerError( Exception ) :
	pass

class PSIPREDRunner :
	def __init__( self, msa=None, query=None, #main input
		number_of_processors=None, database=None, #BLAST Setting if MSA should be built
		echo = True, #echoing ouptut option
		#commands
		formatdb=None, 
		makemat=None, 
		psipred=None, 
		psipass2=None,
		#data directory containing weight file for psipred
		psipred_data_dir=None ) :
		'''
		Predict secondary structure using PSIPRED.

		if multiple sequence alignment or msa (MSA object) is given, 
		the msa is used for prediction.

		if query is given, blastpgp will be used
		for building multiple sequence alignment and
		predict the secondary structure.

		Mode 1 of this class is basically copied from
		runpsipred script.
		'''

		if formatdb == None :
			formatdb = Settings.get( 'formatdb' ) 
		if makemat == None :
			makemat = Settings.get( 'makemat' )
		if psipred == None :
			psipred = Settings.get( 'psipred' )
		if psipass2 == None :
			psipass2 = Settings.get( 'psipass2' )
		if psipred_data_dir == None :
			psipred_data_dir = Settings.get( 'psipred_data_dir' )

		
		print("PSIPRED is running...")

		self.temp_dir = None #tempfile.tempname()
		self.number_of_processors = number_of_processors


		#mode 2 stuff.
		self.database = database
		self.query = query

		#Mode 1. use the given MSA to predict 2nd Structures
		if msa != None :

			input_string = str(msa.query)
			self.temp_dir = tempfile.mkdtemp()

			alignment_input = os.path.join( self.temp_dir, "query.aln" )
			fp = open( alignment_input, 'w' )
			msa.build_psiblast_alignment_input( fp )
			fp.close()

			dummy_db = os.path.join( self.temp_dir, 'query.seq' )
			msa.query.save( dummy_db )
			os.system( formatdb + ' -i ' + dummy_db )

			checkpoint = os.path.join( self.temp_dir, 'query.chk' )
			output = '/dev/null'

			if verbose :
				print('temp_dir:', self.temp_dir)
				print('input_string:', input_string)
				print('dummy_db', dummy_db)
				print('checkpoint', checkpoint)
				print('alignment_input', alignment_input)
				print(open( alignment_input).read())

			runner = PSIBLASTRunner( input_string=input_string, 
				max_iterations=1,
				output = output,
				input_alignment = alignment_input, 
				database = dummy_db, 
				number_of_processors = number_of_processors,
				checkpoint = checkpoint )
			runner.run()

			basename = os.path.join( self.temp_dir , 'query' )
				
			fp = open( basename+".pn", 'w' )
			print("query.chk", file=fp)
			fp.close()

			fp = open( basename+".sn", 'w' )
			print("query.seq", file=fp)
			fp.close()

			if verbose :
				print("basename:", basename)
				print(basename+'.pn')
				print(open(basename+'.pn').read())
				print(basename+'.sn')
				print(open(basename+'.sn').read())

			os.system(makemat + " -P " + basename);
			weight_file = os.path.join( psipred_data_dir, "weights.dat" )
			os.system( psipred + " %(basename)s.mtx %(weight_file)s %(weight_file)s2 %(weight_file)s3 %(weight_file)s4 > %(basename)s.ss" % locals() )

			os.system("%(psipass2)s %(psipred_data_dir)s/weights_p2.dat 1 0.98 1.09 %(basename)s.ss2 %(basename)s.ss > %(basename)s.horiz" % locals() );

			self.output = basename + '.horiz' #important output !!!

		#Mode 2. build MSA and then predict 2nd structures
		elif query != None :
			raise NotYetImplementedError( 'Building Query Mode has not yet been implemented!' )
	
		else :
			raise PSIPREDRunnerError( 'MSA or query FASTA should be given!' )


	def __del__( self ) :
		if self.temp_dir :
			shutil.rmtree( self.temp_dir )

def add_psipred( fasta_file, ss, conf ) :
	'''
	Adds PSIPRED prediction results into the a3m or other fasta file
	formatted for hhmake to build HHM file.
	'''
	fp = open( fasta_file, 'a' )
	print('>ss_pred', file=fp) 
	print(ss, file=fp)
	print('>ss_conf', file=fp)
	print(conf, file=fp)
	fp.close()


def remove_query_gaps( read_fn, save_fn ) :
	'''
	This function removes gaps in the query sequence
	in the simple MSA format (as in the PSIBLAST Input option B).

	WARNING: This function relies in the model that the input file has
	32 character length MSA.
	'''
	content = open( read_fn ).readlines()
        label_length = 32
        #get columns to be deleted
        columns = [1]*len(content[0])
        for i, c in enumerate( content[0] ) :
                if i >= label_length and c == '-' : 
                        columns[ i ] = 0

        new_content = []
        for l in content :
                new_line = ''
                for i, c in enumerate(l) :
                        if columns[i] :
                                new_line += c
                new_content.append( new_line )

        open( save_fn, 'w' ).writelines( new_content )

class ProfileBuildingError( Exception ) :
	pass

if __name__ == '__main__' :
	import sys
	buildali= BuildAli( cmd='/home/kim/local/hhsearch1.5_2/buildali_casp8.pl', sequence = sys.argv[2], max_iterations=int(sys.argv[1]), save_dir=sys.argv[3], selected_iterations=[] )

	buildali.run()
