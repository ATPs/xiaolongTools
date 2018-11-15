import os, sys, tempfile, shutil, glob, pickle
from subprocess import Popen, PIPE

from evdblib.DBTools import Settings
from evdblib.Runners import Runner
from evdblib.Utils import parse_sequence_filename, parse_profile_filename, build_profile_filename
from evdblib.Utils.Parsers.Range import SequenceRange
from evdblib.Utils.Parsers.FASTA import FASTA

verbose = 0
preserve_temp_dir = 0

class BuildAliRunnerError( Exception ) :
	pass

class BuildAliRunner :
        '''
        Build alignments using buildali.pl written by Soding group in Germany.

        Note that this method is not able to run background.
        Also delecate control of temporary directory is not required since
        buildali program uses temporary directory by itself. 
        '''
        def __init__( self, cmd=None, 
                        save_dir=None, sequence=None, range=None, 
			max_iterations=1, save_all_iteration_results=False, msa_input_fn=None,
                        max_rerun=3, stdout=PIPE, stderr=PIPE, timeout=0, number_of_processors=1 ) :
		'''
		sequence: fasta file for the query sequence
		range: 
		'''

		self.temp_dir = None #initialized at the first
		if save_dir == None :
			self.save_dir = tempfile.mkdtemp()
			self.temp_dir = self.save_dir
		else :
			self.save_dir = save_dir
			
                if cmd == None :
                        cmd = Settings.get( "buildali" )

                self.cmd = cmd
		self.number_of_processors = number_of_processors

		self.sequence = sequence
		self.msa_input_fn = msa_input_fn
		if msa_input_fn :
			self.sequence = msa_input_fn
                self.range = range

		print(self.sequence, self.msa_input_fn)

		if self.sequence or self.msa_input_fn :
			pass
		else :
                        raise BuildAliRunnerError( "Query sequence or MSA is necessary to build profile." )
	
                if (self.sequence and os.path.exists( self.sequence )) or (self.msa_input_fn and os.path.exists(self.msa_input_fn) ):
			pass
		else :
                        raise ProfileBuildingError( "Query sequence file %s cannot be found!"%self.sequence )
			
                self.max_rerun = max_rerun
		self.save_all_iteration_results = save_all_iteration_results

		#buildali options
                #getting profile iteration information from setitng
                self.max_iterations = max_iterations
		self.msa_input_fn = msa_input_fn

                #getting command line arguments and setting timedrunner
                self.command_line = self.get_command_line()
                #Runner.TimedRunner.__init__( self, self.command_line, **kwargs )
                self.timeout = timeout
                self.stdout = stdout
                self.stderr = stderr


	def __del__( self ) :
		if (not preserve_temp_dir) and self.temp_dir :
			shutil.rmtree( self.temp_dir )

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

                sequencerange = SequenceRange(range)

                fasta = FASTA(self.sequence)
                extracted_fasta_string = fasta.extract_fasta( sequencerange )

                fp = tempfile.NamedTemporaryFile()
                fp.write( extracted_fasta_string )
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


	def get_final_iteration_file( self ) :
		if self.a3m_files :
			return self.a3m_files[-1]
		else :
			return None


        def set_result_files( self ) :
                '''
                returns list of a3m files in the save_dir.
                '''
                basename = parse_profile_filename( self.sequence )[1]
                if not basename :
                        return
		
		a3m_files = glob.glob( os.path.join( self.save_dir, basename + '*.a3m' ) )
		new_a3m_files = []
		for a3m in a3m_files :
			dir, basename2, iteration, ext = parse_profile_filename( a3m )
			newa3m = build_profile_filename( dir, basename, iteration, ext )
			shutil.move( a3m, newa3m )
			new_a3m_files.append( newa3m )

                self.a3m_files = new_a3m_files
		self.a3m_files.sort()


        def get_command_line( self ) :

		if self.save_all_iteration_results :
			cn = '-cn'
		else :
			cn = ''

		if self.number_of_processors > 1 :
			proc = '-cpu %s'%self.number_of_processors
		else :
			proc = ''

		inputfn=self.sequence

                return '%s -n %d %s %s %s %s '%(self.cmd, int(self.max_iterations), proc, cn, inputfn, self.save_dir  )

class IterativeBuildAliRunnerError( Exception ) :
	pass

from evdblib.Utils.Parsers.BLAST.MSA import MSA

class IterativeBuildAliRunner :
        def __init__( self, sequence=None, range=None, save_dir=None, max_iterations=8, msa_input_fn=None, **kwargs ) :
		'''
		sequence: buildali.pl query sequence fasta file.
		range: sequence range string
		save_dir: save directory, if None same directory of sequence file.
		max_iterations: 
		msa_input_fn: jump start a3m file, but does not implemented yet.
		'''
		self.temp_dir = tempfile.mkdtemp()

		if sequence == None :
			raise IterativeBuildAliRunnerError( "No input sequence is given!" )

		self.output_base = os.path.basename(sequence) 
		self.output_base = parse_sequence_filename(self.output_base)[1]
		sequence = os.path.abspath( sequence )

		#setting up saving directory
		self.save_dir = save_dir
		if not self.save_dir :
			self.save_dir = os.getcwd()

		if not os.path.exists( self.save_dir ) :
			os.makedirs( self.save_dir )

		self.sequence=sequence
		self.range = range
		if self.range :
			self.sequencerange = SequenceRange()
			self.sequencerange.parse(range)
		else :
			self.sequencerange = None

		self.buildali_arguments	 = kwargs
		self.max_iterations = max_iterations

		self.current_runner = None
		self.current_iteartion = 0
		self.current_msa_output = None

		self.a3m_files = []
		self.aln_files = []

		#################################
		#compatibility for IterativeBLAST
		self.input = sequence
		self.input_dir, self.input_base, self.input_ext = parse_sequence_filename( self.input )
		#convert into FASTA
		self.inputfasta = FASTA( self.input )
		self.inputfasta.sequencerange = self.sequencerange
		self.processed_inputfasta = self.inputfasta.extract_fastas()
		#process input fasta
		if self.sequencerange :
			self.processed_input_string = self.inputfasta.extract_fasta_string()
			self.processed_input = os.path.join( self.temp_dir, self.input_base+self.input_ext )
			fasta = FASTA()
			fasta.read_from_string( self.processed_input_string)
			fasta.save( self.processed_input )
		else :
			self.processed_input = self.input


	def __del__( self ) :
		if not preserve_temp_dir :
			shutil.rmtree( self.temp_dir )

	def get_msa( self, iteration ) :
		'''
		return MSA for the iteration
		'''
		if iteration == -1 :
			keys = list(self.msas.keys())
			keys.sort()
			return self.msas[ keys[-1] ]
		return self.msas.get( iteration )

	def get_a3m_file( self, iteration ) :
		'''
		return alignment file for the given iteration number.
		If the iteration file is not available, returns None

		Note that this function retuns a3m file in the temporary directory!
		'''
		try :
			return self.a3m_files[iteration-1]
		except IndexError :
			return None
	

	def initial_iteration( self ) :
		buildalirunner = BuildAliRunner( sequence=self.processed_input, 
			max_iterations=1,
			**self.buildali_arguments )
		buildalirunner.run()

		self.current_iteration = 1
		self.current_runner = buildalirunner
		
		a3mfile = buildalirunner.get_final_iteration_file()

		msa = MSA( blast_result_fp=open(a3mfile), buildalia3m=True )
		msa.parse_a3m_fasta()

		self.current_msa_output = os.path.join( self.temp_dir, '%s.1.a3m'%self.output_base )
		if os.path.exists( self.current_msa_output ) :
			raise Exception( self.current_msa_output + " already exists!" )

		msa.build_a3m_fasta( open(self.current_msa_output, 'w') )
		self.a3m_files.append( self.current_msa_output )

		return msa

	def next_iteration( self ) :
		buildalirunner = BuildAliRunner( msa_input_fn=self.current_msa_output, 
			max_iterations=1,
			**self.buildali_arguments )
		buildalirunner.run()

		self.current_iteration += 1
		self.previous_runner = self.current_runner
		self.current_runner = buildalirunner
		
		a3mfile = buildalirunner.get_final_iteration_file()
	
		msa = MSA( blast_result_fp=open(a3mfile), buildalia3m=True )
		msa.parse_a3m_fasta()

		self.current_msa_output = os.path.join( self.temp_dir, '%s.%d.a3m'%(self.output_base,self.current_iteration) )
		if os.path.exists( self.current_msa_output ) :
			raise Exception( self.current_msa_output + " already exists!" )

		msa.build_a3m_fasta( open(self.current_msa_output, 'w') )
		self.a3m_files.append( self.current_msa_output )

		return msa

	def check_convergence( self, previous, current ) :
		'''
		returns true if the previous msa has more elements than the current msa.
		'''
		if len(previous) >= len(current) :
			self.converged = True
			return True

		else :
			self.converged = False
			return False


	def cleanup( self ) :
		'''
		puts necessary data files into save_dir
		'''
		for i, fn in enumerate(self.a3m_files) :
			shutil.copy( fn, self.save_dir )
			self.a3m_files[i] = os.path.join( self.save_dir, os.path.basename(fn) )
		
        def run( self, echo=True ) :
                '''
                Main function to initiate the iterative blast runs.

                Most behaviors are defined as member functions.
                1. check_convergence()
                2. initial_iteration()
                3. next_iteration()
                4. cleanup()
                '''
                self.msas = {}

                if echo :
                        print("Iteration 1 starting...")

                #run initial blast
                current_msa = self.initial_iteration()
                self.msas[1] = current_msa

                for i in range( 2, self.max_iterations+1 ) :
                        if echo :
                                print("\nIteration", i, "starting...")

                        previous_msa = current_msa
                        current_msa = self.next_iteration()

                        ##################
                        #check convergence
                        ##################
                        #need to implemented
                        if self.check_convergence(previous_msa, current_msa) :
                                if echo :
                                        print("\nConverged!\n")
                                break
                        self.msas[i] = current_msa

                self.cleanup()

from evdblib.Runners.BLASTRunners import build_neighboring_msas, ScoreMat, BackblastPurger

class IterativePurgedOverlappingHSPsBuildAliRunner( IterativeBuildAliRunner ) :
	def __init__( self, sequence=None, range=None, neighboring_fastas=None, neighboring_msas=None,
			use_backblast_purging=True, use_overlapping_purging=True, max_iterations=8,
			number_of_processors=1, **kwargs ) :
		'''
                Removes the erroneously extended fragments first by removing
                backblasting, and then overlapping region.
                Also Known as HangOut.
                For detils see HangOut paper published in Bioinformatics.

                If neighboring_fastas (FASTA files fro neighboring sequences) 
                or neightboring_msas (Multiple sequence alignments for neighboring sequences) are given,
                those informations are used.

                If not, neighboring information is extracted from the sequencerange.
		Note that this class is very similar to that of original HangOut method
		but the profiles are built by buildali.pl method in HHsearch package.
                '''
		IterativeBuildAliRunner.__init__( self, 
			sequence=sequence, range=range,
			number_of_processors=number_of_processors, 
			**kwargs )
		self.kwargs = kwargs
		self.number_of_processors = number_of_processors

		self.use_backblast_purging = use_backblast_purging
		self.use_overlapping_purging = use_overlapping_purging
		self.inserted_positions = []
		self.neighboring_msas = []
		
		if not self.sequencerange :
			raise Exception( "Sequence range should be given!" )
		else :
			self.processed_inputfasta = self.inputfasta.extract_fastas()

			processed_sequencerange = self.processed_inputfasta.sequencerange
			self.inserted_positions = processed_sequencerange.get_breaking_positions()
			self.inserted_positions = [i-1 for i in self.inserted_positions ]

		if not self.inserted_positions :
			print("No insertion positions found! Overlapping purging cannot be done!", file=sys.stderr)
			self.use_overlapping_purging = False

		#build neighboring msas
		if self.use_backblast_purging :
			if verbose : print("Inferring neighboring sequence for backblasting")
			self.neighboring_msas.extend( build_neighboring_msas( input_fasta=self.inputfasta ) )

			if neighboring_fastas :
				self.neighboring_msas.extend( build_neibhboring_msas( neighboring_fastas=neighboring_fastas ))
			if neighboring_msas :
				self.neighboring_msas.extend( neighboring_msas )
				
			if not self.neighboring_msas :
				print("No neighboring sequence found! Backblast purging cannot be done!", file=sys.stderr)
				self.use_backblast_purging = False

			if verbose :
				print("Backblast purging preparation is done!")
				print("# of neighboring msas:", len(self.neighboring_msas))


	def initial_iteration( self, echo=True ) :
		'''
		Run initial buildali and prepare output
		'''
		runner = BuildAliRunner( sequence=self.processed_input, 
			number_of_processors=self.number_of_processors, 
			max_iterations=1,
			**self.kwargs )
		runner.run()

		self.current_iteration = 1
		self.current_runner = runner
		
		#parse a3m output
		a3mfile = runner.get_final_iteration_file()
		msa = MSA( blast_result_fp=open(a3mfile), buildalia3m=True )
		msa.parse_a3m_fasta()

		self.current_msa = msa

		#msa.set_combine_hsps() #not sure to work.
		
		#purge the result
		if self.use_overlapping_purging :
			if echo :
				print('purging overlapping regions...')
			#build pssm
			#for initial iteration, BLOSUM62 matrix is used.
			pssm = ScoreMat()
			pssm.set_blosum_mat()
			if self.number_of_processors > 1 :
				msa.purge_overlapping_hsps_multithreading( self.inserted_positions, 
					pssm, self.number_of_processors)
			else :
				msa.purge_overlapping_hsps( self.inserted_positions, pssm )

		if self.use_backblast_purging :
			for neighbor_msa in self.neighboring_msas :
				backblastpurger = BackblastPurger( msa, neighbor_msa )	

		self.current_msa_output = os.path.join( self.temp_dir, '%s.1.a3m'%self.output_base )
		if os.path.exists( self.current_msa_output ) :
			raise Exception( self.current_msa_output + " already exists!" )

		msa.build_a3m_fasta( open(self.current_msa_output, 'w') )
		self.a3m_files.append( self.current_msa_output )

		return msa


	def next_iteration( self, echo=True ) :
		'''
		Run next buildali and prepare output
		'''
		runner = BuildAliRunner( msa_input_fn=self.current_msa_output, 
			number_of_processors=self.number_of_processors, 
			max_iterations=1,
			**self.kwargs )
		runner.run()

		self.current_iteration += 1
		self.previous_runner = self.current_runner
		self.previous_msa_output = self.current_msa_output
		self.previous_msa = self.current_msa
		self.current_runner = runner

		i = self.current_iteration
		
		#parse a3m output
		a3mfile = runner.get_final_iteration_file()
		msa = MSA( blast_result_fp=open(a3mfile), buildalia3m=True )
		msa.parse_a3m_fasta()
		self.current_msa = msa

		#msa.set_combine_hsps() #not sure to work.
		
		#purge the result
		if self.use_overlapping_purging :
			if echo :
				print('purging overlapping regions...')
			#build pssm
			#for initial iteration, BLOSUM62 matrix is used.
			pssm = ScoreMat()
			pssm.build_pssm( self.previous_msa )
			if self.number_of_processors > 1 :
				msa.purge_overlapping_hsps_multithreading( self.inserted_positions, 
					pssm, self.number_of_processors)
			else :
				msa.purge_overlapping_hsps( self.inserted_positions, pssm )

		if self.use_backblast_purging :
			for neighbor_msa in self.neighboring_msas :
				backblastpurger = BackblastPurger( msa, neighbor_msa )	

		self.current_msa_output = os.path.join( self.temp_dir, '%s.%d.a3m'%(self.output_base,self.current_iteration) )
		if os.path.exists( self.current_msa_output ) :
			raise Exception( self.current_msa_output + " already exists!" )

		msa.build_a3m_fasta( open(self.current_msa_output, 'w') )
		self.a3m_files.append( self.current_msa_output )

		return msa
