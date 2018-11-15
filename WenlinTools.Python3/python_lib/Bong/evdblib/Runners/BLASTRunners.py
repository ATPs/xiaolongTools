'''
This module contains wrappers for the blastpgp program.

BLASTRunner class will run protein blast or psiblast.
'''
import os, sys, time, tempfile, shutil, copy
from subprocess import Popen, PIPE
from io import StringIO
from evdblib.Utils.Parsers import FASTA, BLAST
from evdblib.Utils import parse_sequence_filename, find_command_in_path, build_profile_filename

from . import Runner

formatdb = 'formatdb'
if find_command_in_path( formatdb ) :
	from evdblib.DBTools import Settings
	formatdb = Settings.get( 'formatdb'  )

verbose = 0
default_max_iterations_for_neighbors = 1

class BLASTRunner :
	'''
	BLASTRunner class will run blastpgp program.
	'''

	def __init__( self, input=None, output=None, database=None, evalue_cutoff=None, 
			input_string=None, input_alignment=None, profile_output=None, 
			checkpoint=None,
			summary_table_size=None, number_of_hits_to_display=None, 
			output_format=None, number_of_processors=None,
			effective_database_length=None,
			use_my_default=True, max_rerun = None, range=None,
			show_gi=True, blastcmd='blastpgp', timeout=0 ) :
		'''
		BLAST runner accepts the following options.

		input: BLAST input fasta file (-i option)
		input_string: Piped input of fasta file as a python string object.
		output: BLAST result file (-o option)
		database: -d option
		evalue_cutoff: cutoff value for display a hit. (-e option)
		summary_table_size: number of hits display in the summary table (-b option)
		number_of_hits_to_display: number of hits to dispaly in the alignment section (-v)
		output_format: Output file format. (-m)
		number_of_processors: default 1 (-a)
		show_gi: default True (-I)
		
		blastcmd: default 'blastpgp' in the path
		timeout: timeout value to kill the blast run in seconds.
		max_rerun: how many times to rerun the BLAST command 
			if the command does not produce output.
		range: SequenceRange object for running BLAST only
			on the given region.

		Among the above options input and output are the two options 
		that are requied to run blastpgp program.

		Note that the pipe is not allowed currently.
		'''

		self.input = input
		self.input_string = input_string

		self.input_alignment = input_alignment
		self.profile_output = profile_output
		self.checkpoint = checkpoint
		self.output = output
		self.database = database
		self.evalue_cutoff = evalue_cutoff

		self.summary_table_size =  summary_table_size
		self.number_of_hits_to_display = number_of_hits_to_display

		self.output_format  = output_format
		self.number_of_processors =  number_of_processors
		self.effective_database_length = effective_database_length
		self.show_gi = show_gi

		self.blastcmd = blastcmd 
		self.timeout = timeout

		self.sequencerange = range

		##############################################
		#prepare input fasta!
		self.inputfasta = FASTA.FASTA()

		if input_string :
			self.inputfasta.read_from_string( input_string )
		else :
			fp = open( self.input )
			self.inputfasta.read( fp )
			fp.close()

		self.inputfasta.sequencerange = self.sequencerange

		self.processed_inputfasta = self.inputfasta.extract_fastas()
		##############################################
		#processed_inputfasta should be the main interface 
		#to get the input file!!
		###############################################

		if max_rerun == None :
			max_rerun = 0
		self.max_rerun = int( max_rerun )

		##########################
		#setting My default values
		##########################
		self.use_my_default = use_my_default
		self.set_default()

		self.command_line = self.get_command_line()

		if os.path.exists( self.blastcmd ) :
			pass
		elif find_command_in_path( self.blastcmd ) :
			self.blastcmd = find_command_in_path( self.blastcmd )
		else :
			raise BLASTOptionError( "BLAST %s is not found." % self.blastcmd )

	def _initiate_runner( self ) :
		'''
		returns runner class.
		'''
		return Runner.TimedRunner( self.command_line, timeout=self.timeout, stdin=PIPE, stdout=PIPE, stderr=PIPE )



	def set_default( self ) :
		'''
		puts Non-BLAST default but default but for my setting.

		This functions sets the default evalue cutoff into 0.001,
		number of alignments to show into 30,000
		and summary table of hits into 100.
		'''
		if not self.use_my_default :
			return

		if self.evalue_cutoff == None :
			self.evalue_cutoff = 0.001
		if self.number_of_hits_to_display == None :
			self.number_of_hits_to_display = 30000
		if self.summary_table_size == None  :
			self.summary_table_size = 100
	

	def get_command_line( self ) :
		'''
		returns formatted command line arguments.

		Note that the missing input and output options will
		raise BLASTOptionError Exception.
		'''
		fp = StringIO()

		print(self.blastcmd, end=' ', file=fp)
		
		if self.output == None :
			raise BLASTOptionError( "No output file is given." )
		else :
			print('-o', self.output, end=' ', file=fp)

		if self.database != None :
			print('-d', self.database, end=' ', file=fp)

		if self.evalue_cutoff !=None :
			print('-e', self.evalue_cutoff, end=' ', file=fp)

		if self.summary_table_size != None :
			print('-v', self.summary_table_size, end=' ', file=fp)

		if self.number_of_hits_to_display != None :
			print('-b', self.number_of_hits_to_display, end=' ', file=fp)

		if self.output_format :
			print('-m', self.output_format, end=' ', file=fp)

		if self.number_of_processors != None :
			print('-a', self.number_of_processors, end=' ', file=fp)

		if self.input_alignment :
			print('-B', self.input_alignment, end=' ', file=fp)

		if self.profile_output :
			print('-Q', self.profile_output, end=' ', file=fp)
		if self.checkpoint :
			print('-C', self.checkpoint, end=' ', file=fp)

		if self.show_gi :
			print('-I T', end=' ', file=fp)
		else :
			print('-I F', end=' ', file=fp)

		return fp.getvalue()


	def run( self, background=False, echo=True ) :
		'''
		Start to run the task.

		If the background is True then the task will run in backgrounds.
		If not, the task will run foreground.

		Note that wait() can be used to syncronize the task.

		Note that the max_rerun arugment can only be used when
		background is false, since it is hard to check with background status.
		'''
		if echo :
			print(str(self.processed_inputfasta).strip())
			print(self.command_line)

		self.runner = self._initiate_runner()
		self.runner.run() 
		stdout, stderr = self.runner.communicate( input=str(self.processed_inputfasta) )
		
		#p = Popen( self.command_line.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE )
		#stdout, stderr = p.communicate( input=str(self.processed_inputfasta) )
		print(stdout)
		print(stderr, file=sys.stderr)

		if not background :
			run_count = 1

			while not (os.path.exists( self.output ) and os.path.getsize( self.output )) and run_count <= self.max_rerun :
				self.runner=self._initiate_runner() #needs to be reinitiated!
				self.runner.run() #input=str(self.processed_inputfasta) )
				#self.runner.communicate( input=str(self.processed_inputfasta) )
				#p = Popen( self.command_line.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE )
				#stdout, stderr = p.communicate( input=str(self.processed_inputfasta) )
				print(stdout)
				print(stderr, file=sys.stderr)

				run_count += 1
	

	def parse( self ) :
		'''
		Parse blast result and returns BLAST object.
		
		For details of the BLAST object,
		refer evdblib.Utils.Parsers.BLAST package.
		'''
		blast = BLAST.BLAST( query=self.processed_inputfasta, blast_result_fn=self.output )
		return blast
		

class PSIBLASTRunner( BLASTRunner ) :
	'''
	This class will run PSIBLAST using blastpgp program.
	'''

	def __init__( self, max_iterations=None, inclusion_cutoff=None, save_dir=None, profile_output=None, **kwargs ) :

		if max_iterations :
			try :
				self.number_of_iterations = int(max_iterations)
			except ValueError :
				raise BLASTOptionError( "Number of iteration should be given as integer value")

		output = kwargs.get( 'output' )
		if output == None and save_dir != None and kwargs.get('input') :
			dir, basename, ext = parse_sequence_filename( kwargs.get('input') )
			kwargs['output'] = os.path.join( save_dir , basename + '.psi' )

		self.output = output
				
		self.inclusion_cutoff = inclusion_cutoff
		self.save_dir = save_dir
		self.profile_output = profile_output
	
		if self.save_dir != None :
			if os.path.exists( save_dir ) :
				if not os.path.isdir(save_dir) :
					raise BLASTOptionError( "Save directory is a file, not a directory. " + self.save_dir )
			else :
				os.makedirs( save_dir )

		BLASTRunner.__init__( self, **kwargs )

		#iteration: aln_file dictionary.
		#this dictionary is important for easy retrieving of the 
		#alignment files for further processing..
		self.aln_files = {} 

		self.get_command_line = self.get_command_line()
		#self.runner = Runner.TimedRunner(self.command_line, timeout=self.timeout, stdin=PIPE, stderr=PIPE, stdout=PIPE )

	

	def get_command_line( self ) :
		fp = StringIO()

		print(BLASTRunner.get_command_line( self ), end=' ', file=fp)
		if self.number_of_iterations != None :
			print('-j', self.number_of_iterations, end=' ', file=fp)
		if self.profile_output != None :
			print('-Q', self.profile_output, end=' ', file=fp)
		if self.inclusion_cutoff != None :
			print('-h', self.inclusion_cutoff, file=fp)
	
		return fp.getvalue()

	def _save_aln_results( self, save_dir=None ) :
		'''
		first parse the blast result file,
		and then save alignment files into the save_dir.
		
		if save_dir is not defined,
		it does not perform any operations.
		'''

		if self.save_dir == None and save_dir == None:
			return
		elif save_dir == None :
			save_dir = self.save_dir
			
		blast = self.parse()
		self.msas = {}

		for i, msa in enumerate(blast) :
			iter = i+1
			dir, basename, ext = parse_sequence_filename( self.input )
			profile_name = build_profile_filename( save_dir, basename, iter, '.aln' )
			profile_fp = open( profile_name, 'w' )
			msa.build_psiblast_alignment_input( profile_fp )
			profile_fp.close()

			self.aln_files[iter] = profile_name
			self.msas[iter] = msa


	def get_aln_file( self, iteration ) :
		'''
		return alignment file for the given iteration number.
		If the iteration file is not available, returns None.
		'''
		
		if not self.aln_files :
			self._save_aln_results()
		
		return self.aln_files.get(iteration)


	def get_msa( self, iteration ) :
		'''
		returns MSA for the iteration
		'''
		if iteration == -1 :
			keys = list(self.msas.keys())
			keys.sort()
			return self.msas[ keys[-1] ]

		return self.msas.get( iteration )


		


class IterativeBLASTRunner :
	'''
	This class will run the blastpgp program iteratively.
	And generate BLAST result equivalent to that of PSIBLAST final iteration
	result.

	Note that this class is running BLASTRunner Iteratively.
	'''

	def __init__( self, input=None, output=None, range=None,
			max_iterations=None, 
			temp_dir=None,
			save_dir=None, **kwargs ) :
		'''
		Runner class of variant of PSIBLAST named as Iterative BLAST.
		
		To manipulate the procedure building MSA, 
		each iteration of BLAST is exposed.

		This class is proof of the concept for other classes like HangOut.
		This class is also a base class for most of other modified 
		blast runs.
		The initial_iteration and next_iteration functions needs to be rewritten
		to extend the functionality of this class.
		

		Most of keyword arguments are same as BLASTRunner,
		because mainly BLASTRunner will be used to build the profile.
		'''

		self.kwargs = kwargs
		self.current_iteration = 1 #initialize

		#getting_multiple processor info
		if 'number_of_processors' in self.kwargs :
			self.number_of_processors = self.kwargs['number_of_processors']
		else :
			self.number_of_processors = 1

		#setting temporary working directory
		self.temp_dir_to_be_removed = None
		if temp_dir == None :
			temp_dir = tempfile.mkdtemp() #need to be removed by deconstructor.
			self.temp_dir_to_be_removed = temp_dir
		self.temp_dir = temp_dir
		self.max_iterations = int(max_iterations)

		intermediate_result_dir = save_dir
		#setting up intermediate result directory.
		#if the directory is made, the result will be saved in it.
		if intermediate_result_dir != None :
			if not os.path.exists( intermediate_result_dir ) :
				os.makedirs( intermediate_result_dir )
				if not os.path.exists( intermediate_result_dir ) :
					raise BLASTOptionError( "Intermediate result saving directory cannot be made.") 
		else :
			BLASTOptionError( "Alignment saving directory (save_dir) should be given!" )
		self.intermediate_result_dir = intermediate_result_dir


		#######################
		#processing INPUT file
		#######################
		if input == None :
			raise BLASTOptionError( "Input file should be given." )
		elif input and not os.path.exists( input ) :
			raise BLASTOptionError( "Input file %s does not exists."%self.input )
		
		#analyze the input filename
		self.input=input
		self.input_dir, self.input_base, self.input_ext = parse_sequence_filename( self.input )

		#converted into FASTA.
		self.inputfasta = FASTA.parse( self.input )  #raw input
		self.inputfasta.sequencerange = range

		self.processed_inputfasta = self.inputfasta.extract_fastas()

		#process input file
		self.sequencerange = range
		if self.sequencerange != None :
			self.processed_input_string = self.inputfasta.extract_fasta_string()
			self.processed_input = os.path.join( self.temp_dir, self.input_base+self.input_ext )
			fasta = FASTA.FASTA()
			fasta.read_from_string( self.processed_input_string )
			fasta.save( self.processed_input )

			if verbose :
				print("IterativeBLASTRunner.__init__")
				print("input fasta")
				print(self.inputfasta)
				print("processed_inputfasta")
				print(self.processed_inputfasta)
				print("processed inputfasta string")
				print(self.processed_input_string)
				print("processed_input")
				print(open(self.processed_input).read())
		else :
			self.processed_input = self.input

		if verbose :
			print('Finally use fasta...')
			print(self.processed_input)
			print(open(self.processed_input).read())

		#the input file should be used is the self.processed_input

		#######################
		#processing Output file
		#######################

		self.output=output
		if output == None :
			dir, fnbase, ext = parse_sequence_filename( self.input )
			self.output_base, self.output_ext = fnbase, '.bla'
		else :
		#analyze the output filename
			dir, self.output_base, self.output_ext = parse_sequence_filename( self.output )

		self.aln_files = {} #iteration:aln_filename
		

	#This function will make sure the 
	#deletion of the temporary directory.
	def __del__( self ) :
		'''
		Deconstructor.
		This function simply deletes the temporary directory for clean up.
		'''
		if self.temp_dir_to_be_removed :
			shutil.rmtree( self.temp_dir_to_be_removed )

	def check_convergence( self, previous, current ) :
		'''
		returns true if the previous msa has more element than current msa.
		This function mgiht be rewritten to inherit this class
		'''
		if len(previous) >= len(current) :
			self.converged = True
			return True
		else :
			self.converged = False
			return False

	def initial_iteration( self ) :
		'''
		Runs BLAST in a modified way as a initial iteration.
		Modify this function for inherited class.
		'''
		temp_output = os.path.join( self.temp_dir, '%s.1%s'%(self.output_base,self.output_ext  ) )
		runner = BLASTRunner( input=self.processed_input, output=temp_output, **self.kwargs )
		runner.run()

		blast = BLAST.BLAST( self.processed_input, runner.output )

		msa_output = os.path.join( self.temp_dir,  '%s.1.aln'%(self.output_base) )
		fp = open( msa_output, 'w' )
		blast[0].build_psiblast_alignment_input( fp )
		fp.close()

		self.current_runner=runner
		self.current_msa_output=msa_output
		self.current_iteration = 1
		#self.current_parser = blast

		self.msa_files.append( msa_output )
		return blast[0] #msa

	def next_iteration( self ) :
		'''
		Runs BLAST as a sebsequent iteration after initial iteration.
		Modify this function for inherited class.
		'''
		self.current_iteration += 1
		self.previous_runner = self.current_runner
		self.previous_msa_output = self.current_msa_output
		#self.previous_parser = self.current_parser

		i = self.current_iteration
		temp_output = os.path.join( self.temp_dir, '%s.%s%s'%(self.output_base,str(i),self.output_ext  ) )
		runner = BLASTRunner( input=self.processed_input, output=temp_output, input_alignment=self.previous_msa_output, **self.kwargs )
		runner.run()

		blast = BLAST.BLAST( self.processed_input, runner.output )

		msa_output = os.path.join( self.temp_dir,  '%s.%s.aln'%(self.output_base,str(i)) )
		fp = open( msa_output, 'w' )
		blast[0].build_psiblast_alignment_input( fp ) #need to be implemented
		fp.close()

		#finally
		self.current_runner=runner
		self.current_msa_output = msa_output
		#self.current_parser = blast

		self.msa_files.append( msa_output )
		return blast[0] #msa

	def cleanup( self ) :
		'''
		Puts necessary data files into right places.
		'''
		msa_files = self.msa_files
		if self.converged :
			if self.intermediate_result_dir :
				for i, msa_file in enumerate(msa_files[:-1]) :
					try :
						shutil.move( msa_file, self.intermediate_result_dir )

						dir, fn = os.path.split( msa_file )
						self.aln_files[i+1] = os.path.join( self.intermediate_result_dir, fn )
					except IOError :
						print("WARNING: Intermediate result, %s, cannot be written." % msa_file)
						
			if msa_files and self.output :
				shutil.move( self.previous_runner.output, self.output )
		else :
			if self.intermediate_result_dir :
				for i, msa_file in enumerate( msa_files ) :
					try :
						shutil.move( msa_file, self.intermediate_result_dir )
						dir, fn = os.path.split( msa_file )
						self.aln_files[i+1] = os.path.join( self.intermediate_result_dir, fn )
					except IOError :
						print("WARNING: Intermediate result, %s, cannot be written." % msa_file)
						
			if msa_files and self.output :
				shutil.move( self.current_runner.output, self.output )
				

	def run( self, echo=True ) :
		'''
		Main function to initiate the iterative blast runs.

		Most behaviors are defined as member functions.
		1. check_convergence()
		2. initial_iteration()
		3. next_iteration()
		4. cleanup()
		'''
		self.msa_files = [] #results to save
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


	def get_aln_file( self, iteration ) :
		'''
		returns alignment filename for the iteration.
		'''
		return self.aln_files.get( iteration )

	def get_msa( self, iteration ) :
		'''
		returns MSA for the iteration
		'''
		if iteration == -1 :
			keys = list(self.msas.keys())
			keys.sort()
			return self.msas[ keys[-1] ]

		return self.msas.get( iteration )


class FragmentJoiningBLASTRunner( IterativeBLASTRunner ) :
	'''
	Not yet implemented.
	'''
	pass


class BackblastingPurgedBLASTRunner( IterativeBLASTRunner ) :
	def __init__( self, neighboring_msas=None, **kwargs ) :
		'''
		Removes the erroneously extended fragments by running profile
		got from another blast result against each hit.

		If neighboring_msas are given (list of MSAs for neighboring region),
		'''
		##########################################################
		#Did not implemented!!
		##########################################################
		#Use IterativePurgedOverlappingHSPsBLASTRunner 
		#with use_overlapping_purging option False.
		#Since that class have the backblast component in it, 
		#that option replaced this class.
		##########################################################
		pass

class IterativeRemovedOverlappingHSPsBLASTRunner( IterativeBLASTRunner ) :
	'''
	Removes the erroneously extended fragments by removing all HSPs 
	have overlaps within a hit.

	Also known as RemoveHit in the HangOut paper.
	'''
	
	def initial_iteration( self ) :
		'''
		Runs BLAST in a modified way as a initial iteration.
		'''
		temp_output = os.path.join( self.temp_dir, '%s.1%s'%(self.output_base,self.output_ext  ) )
		runner = BLASTRunner( input=self.processed_input, output=temp_output, **self.kwargs )
		runner.run()

		blast = BLAST.BLAST( self.processed_input, runner.output )

		msa = blast[-1]
		
		#remove hits
		msa.remove_overlapping_hsps()
		msa.set_combine_hsps()

		msa_output = os.path.join( self.temp_dir,  '%s.1.aln'%(self.output_base) )
		fp = open( msa_output, 'w' )
		msa.build_psiblast_alignment_input( fp )
		fp.close()

		self.current_runner=runner
		self.current_msa_output=msa_output
		self.current_iteration = 1
		#self.current_parser = blast

		self.msa_files.append( msa_output )
		return msa


	def next_iteration( self ) :
		'''
		Runs BLAST as a sebsequent iteration after initial iteration.
		'''
		self.current_iteration += 1
		self.previous_runner = self.current_runner
		self.previous_msa_output = self.current_msa_output
		#self.previous_parser = self.current_parser

		i = self.current_iteration
		temp_output = os.path.join( self.temp_dir, '%s.%s%s'%(self.output_base,str(i),self.output_ext  ) )
		runner = BLASTRunner( input=self.processed_input, output=temp_output, input_alignment=self.previous_msa_output, **self.kwargs )
		runner.run()

		blast = BLAST.BLAST( self.processed_input, runner.output )
		msa = blast[-1]

		#remove hits
		msa.remove_overlapping_hsps()
		msa.set_combine_hsps()

		msa_output = os.path.join( self.temp_dir,  '%s.%s.aln'%(self.output_base,str(i)) )
		fp = open( msa_output, 'w' )
		msa.build_psiblast_alignment_input( fp ) #need to be implemented
		fp.close()

		#finally
		self.current_runner=runner
		self.current_msa_output = msa_output
		#self.current_parser = blast

		self.msa_files.append( msa_output )
		return msa

	

class IterativePurgedOverlappingHSPsBLASTRunner( IterativeBLASTRunner ) :
	def __init__( self, neighboring_fastas=None, neighboring_msas=None,
			use_backblast_purging=True, use_overlapping_purging=True, **kwargs ) :
		'''
		Removes the erroneously extended fragments first by removing
		backblasting, and then overlapping region.
		Also Known as HangOut.
		For detils see HangOut paper published in Bioinformatics.

		If neighboring_fastas (FASTA files fro neighboring sequences) 
		or neightboring_msas (Multiple sequence alignments for neighboring sequences) are given,
		those informations are used.

		If not, neighboring information is extracted from the sequencerange.
		'''
		self.kwargs = kwargs

		IterativeBLASTRunner.__init__( self, **kwargs )
		self.use_backblast_purging = use_backblast_purging
		self.use_overlapping_purging = use_overlapping_purging
		self.inserted_positions = []

		if not self.sequencerange :
			raise BLASTOptionError( "Sequence range or insertion positions should be given!" )
		else :
			self.processed_inputfasta = self.inputfasta.extract_fastas()
			processed_sequencerange = self.processed_inputfasta.sequencerange
			self.inserted_positions = processed_sequencerange.get_breaking_positions()

			#adjust for 0 based index.
			self.inserted_positions = [i-1 for i in self.inserted_positions]

			if verbose :
				print("IterativePurgedOverlappingHSPsBLASTRunner.__init__():")
				print("\tprocessed_fasta        :", self.processed_inputfasta)
				print("\tprocessed_sequencerange:", processed_sequencerange)
				print('\tself.inserted_positions:', self.inserted_positions)

		if not self.inserted_positions :
			print("No insertion positions found! Overlapping purging cannot be done!", file=sys.stderr)
			self.use_overlapping_purging = False

		#build neibhoring msas
		if self.use_backblast_purging :
			if not neighboring_fastas  and not neighboring_msas  :
				if verbose : print("Inferring neighboring sequences for backblasting...")
				self.neighboring_msas = build_neighboring_msas( input_fasta=self.inputfasta, **kwargs )

			elif neighboring_fastas :
				if verbose : print("Using neighboring sequences for backblasting...", len(neighboring_fastas))
				self.neighboring_msas = build_neighboring_msas( neighboring_fastas=neighboring_fastas, **kwargs )

			elif neighboring_msas :
				self.neighboring_msas = neighboring_msas

			if not self.neighboring_msas :
				print("No neighboring positions found! Backblast purging cannot be done!", file=sys.stderr)
				self.use_backblast_purging = False

			if verbose :
				print("Backblast purgin preparation is done!")
				print("# of neighboring msas:", len(self.neighboring_msas))


	def initial_iteration( self, echo=True ) :
		'''
		Run initial blast and prepare output
		'''
		#run blast
		temp_output = os.path.join( self.temp_dir, '%s.1%s'%(self.output_base,self.output_ext  ) )
		runner = BLASTRunner( input=self.processed_input, output=temp_output, **self.kwargs )
		runner.run()

		#parse blast output 
		blast = BLAST.BLAST( self.processed_input, runner.output )
		msa = blast[-1]
		#flagging for combining HSPs for B input alignment.
		msa.set_combine_hsps() 

		#purge the result
		if self.use_overlapping_purging :
			if echo :
				print('purging overlapping regions...')

			#build pssm 
			#for initial iteration, BLOSUM 62 matrix is used!
			pssm = ScoreMat() # len(self.processed_inputfasta) )
			pssm.set_blosum_mat()
			#start to purge the matrix
			if self.number_of_processors > 1 :
				msa.purge_overlapping_hsps_multithreading( self.inserted_positions, pssm, self.number_of_processors )
			else :
				msa.purge_overlapping_hsps( self.inserted_positions, pssm )

		if self.use_backblast_purging :
			for neighbor_msa in self.neighboring_msas :
				#msa.psiblast_purge( neighbor_msa )
				backblastpurger = BackblastPurger( msa, neighbor_msa, **self.kwargs )

		msa_output = build_profile_filename( self.temp_dir,  self.output_base, 1, '.aln' )
		msa_output_fp = open( msa_output, 'w' )
		msa.build_psiblast_alignment_input( msa_output_fp )
		msa_output_fp.close()

		self.current_runner=runner
		self.current_msa_output=msa_output
		self.current_iteration = 1
		self.current_parser = blast

		self.msa_files.append( msa_output )
		return msa #msa


	def next_iteration( self, echo=True ) :
		self.current_iteration += 1
		self.previous_runner = self.current_runner
		self.previous_msa_output = self.current_msa_output
		self.previous_parser = self.current_parser

		i = self.current_iteration
		temp_output = os.path.join( self.temp_dir, '%s.%s%s'%(self.output_base,str(i),self.output_ext  ) )
		runner = BLASTRunner( input=self.processed_input, output=temp_output, input_alignment=self.previous_msa_output, **self.kwargs )
		runner.run()

		blast = BLAST.BLAST( self.processed_input, runner.output )
		msa = blast[-1]
		msa.set_combine_hsps()

		if self.use_overlapping_purging :
			if echo :
				print('purging overlapping regions...')
			pssm = ScoreMat()
			pssm.build_pssm( self.previous_parser[-1], **self.kwargs )
			msa.purge_overlapping_hsps( self.inserted_positions, pssm )

		if self.use_backblast_purging :
			for neighbor_msa in self.neighboring_msas :
				#msa.psiblast_purge( neighbor_msa )
				backblastpurger = BackblastPurger( msa, neighbor_msa )

		#make msa output
		msa_output = build_profile_filename( self.temp_dir,  self.output_base, i, '.aln' )
		msa_output_fp = open( msa_output, 'w' )
		msa.build_psiblast_alignment_input( msa_output_fp )
		msa_output_fp.close()

		#finally
		self.current_runner=runner
		self.current_msa_output = msa_output
		self.current_parser = blast

		self.msa_files.append( msa_output )

		return msa


class BackblastPurger :
	def __init__( self, msa, neighboring_msa, echo=True, **kwargs ) :
		'''
		Builds PSIBLAST db for neihboring MSA 
		and purge the alignment.

		This class is designed to be used
		to purge a single MSA using neighboring MSAs.

		For backblasting procedure, 
		BackblastingPurgedBLASTRunner should be used.
		'''

		self.del_names = []

		self.msa = msa
		self.neighboring_msa = neighboring_msa 
		
		#build Hit DB 
		hit_db_fp = tempfile.NamedTemporaryFile()
		hit_db = hit_db_fp.name
		self.build_hit_db( hit_db )
		self.del_names.append( hit_db )

		#build input MSA for blast
		input_alignment_fp = tempfile.NamedTemporaryFile()
		input_alignment = input_alignment_fp.name
		self.neighboring_msa.build_psiblast_alignment_input( input_alignment_fp )

		#prepare temporary output file for blast
		backblast_result_fp = tempfile.NamedTemporaryFile()
		backblast_result = backblast_result_fp.name

		if 'input' in kwargs :
			del kwargs['input']
		if 'input_string' in kwargs :
			del kwargs[ 'input_string' ]
		if 'database' in kwargs :
			del kwargs[ 'database' ]
		if 'output' in kwargs :
			del kwargs[ 'output' ]
		if 'effective_database_length' in kwargs :
			del kwargs[ 'effective_database_length' ]
		if 'input_alignment' in kwargs :
			del kwargs[ 'input_alignment' ]
		if 'profile_output' in kwargs :
			del kwargs[ 'profile_output' ]
		if 'range' in kwargs :
			del kwargs[ 'range' ]

		if echo :
			print("\nStarting Backblast Purging...")
		
		#run blast
		runner = BLASTRunner( input_string=str(self.neighboring_msa.query), output=backblast_result, database=hit_db, input_alignment=input_alignment, effective_database_length=5000, **kwargs )
		runner.run()

		blast = BLAST.BLAST( query=self.neighboring_msa.query, blast_result_fn=backblast_result )
		contaminant_msa = blast[-1]

		if kwargs.get( 'number_of_processors' ) :
			nthreads = kwargs.get( 'number_of_processors' )
		else :
			nthreads = 1

		if nthreads > 1 :
			self.msa.remove_overlapped_hit_regions_multithreading( contaminant_msa, nthreads )
		else :
			self.msa.remove_overlapped_hit_regions( contaminant_msa )

	def build_hit_db( self, hit_db ) :
		'''
		builds hit database for self.msa.
		This function generates blast database for 
		hits in self.msa.

		Then the hit db will be used for 
		purging based on the neighboring msa.
		'''
		fp = open( hit_db, 'w' )
		self.msa.write_hit_fasta( fp )
		fp.close()

		os.system( formatdb+' -i %s'%hit_db)
		

		
	def __del__( self ) :
		for fn in self.del_names :
			os.system( "rm -f %s*"%fn )

class ScoreMat :
        def __init__( self, query_len=0 ) :
                self.length = query_len

                self.columns = [] #PSSM
                self.blosum_mat = {} #Blosum

                self.residue_map = {}

                self.gop = 11 #gap opening panelty
                self.gep = 1  #gap extending panelty


        def build_pssm( self, msa, profile_fn=None, db=None, **kwargs ) :
                '''
                Generate PSSM for a given blast msa.
                1. Run psiblast for a dummy db.
                2. read in the pssm file.

		if profile_fn is given, the given file is used.
		otherwise, make a temporary file and use it.

		If dummy db is given, the file is used,
		otherwise make a temporary dummy file, which is a db containing 
		query sequence, is used.
                '''

		if verbose :
                	print("Building profile for \n", msa.query, "...")

		if profile_fn == None :
			profile_fp = tempfile.NamedTemporaryFile()
			profile_fn = profile_fp.name


		self.db_fn = None
		if db == None :
			db_fp = tempfile.NamedTemporaryFile()
			db_fn = db_fp.name

			msa.query.save( fp=db_fp )
			os.system( formatdb + ' -i '+ db_fn )
			db = db_fn
			self.db_fn = db_fn

		input_alignment_fp = tempfile.NamedTemporaryFile()
		input_alignment = input_alignment_fp.name
		msa.build_psiblast_alignment_input( input_alignment_fp )

		output_fp = tempfile.NamedTemporaryFile()
		output = output_fp.name

		if 'input' in kwargs :
			del kwargs[ 'input' ]
		if 'input_string' in kwargs :
			del kwargs[ 'input' ]
		if 'output' in kwargs :
			del kwargs[ 'output' ]
		if 'range' in kwargs :
			del kwargs[ 'range' ]
		if 'database' in kwargs :
			del kwargs[ 'database' ]
		if 'profile_output' in kwargs:
			del kwargs[ 'profile_output' ]
		if 'input_alignment' in kwargs :
			del kwargs[ 'input_alignment' ]

		#setting length of the query.
		self.length = len(msa.query)
	
                        
                runner = BLASTRunner( 
			input_string=str(msa.query), 
			output=output, 
			profile_output=profile_fn, 
			input_alignment=input_alignment, 
			database=db, 
			**kwargs  )

                runner.run()

                self.read_pssm_fn( profile_fn )

		#cleanup
		if self.db_fn :
			db_fp.close()
			os.system( 'rm -rf %s*'%self.db_fn )
			


        def combine( self, score_mat ) :
                '''
                Add given score_mat as a following part of PSSM
                '''
                pass


        def read_pssm_fn( self, fn ) :

                print('#Start to read profile', fn, file=sys.stderr)

                fp = open( fn )

                fp.readline() #read one line
                line = fp.readline()
                header = fp.readline().split()[:20]

                for i, a in enumerate(header) :
                        self.residue_map[ a.upper() ] = i
                        self.residue_map[ a.lower() ] = i

                for i, line in enumerate(fp) :
                        if line == '\n' :
                                break
                        if verbose :
                                print(line)

                        line = line.split()

                        index, residue = line[:2]
                        if verbose :
                                print('index:', index, 'residue:', residue)

                        if i != int(index)-1 :
                                print("WARNING: ", fn, "has problem!", i, index, file=sys.stderr)

                        self.columns.append( [ int(s) for s in line[2:22] ] )
                        #vector = []
                        #for s in line [2:22] :
                                #vector.append( int(s) )
                        #self.columns.append( vector )

                if self.length != len(self.columns) :
                        raise IndexError( 'The length of query and profile does not match! %s %s'% ( str(len(self.columns)), str(self.length) ) )


        def pssm( self, residue_number, query_residue_type, residue_type ) :
                #query residue type is not used!!
                if residue_type in self.residue_map :
                        s = self.columns[residue_number][ self.residue_map[residue_type] ]
                        if verbose :
                                print('%s pssm[%s][%s]=%s' % (residue_number,query_residue_type, residue_type,s))
                        return s
                else :  
                        print('Warning! residue_number:', residue_number, 'query_residue_type:', query_residue_type, 'residue_type:', residue_type, file=sys.stderr)
                        return 0.0



        def cal_pssm_score( self, hsp, start_pos=0, end_pos=0, matrix='pssm' ) :
                '''
                This function calculate pssm score for given hsp.
                Start and end are not query based residue numbering.
                '''

                score_mat = self.pssm
                if matrix == 'blosum' :
                        score_mat = self.blosum
                elif self.length == 0 :
                        score_mat = self.blosum

		if not hsp :
			return -1000.0 #return arbitarily negative value

                #By default tihs function calculate score for whole area.
                if start_pos == 0 and end_pos == 0 :
                        start_pos = hsp.query_start_index
                        end_pos = hsp.query_end_index

                start = hsp.convert_position_to_real_index( start_pos )
                end = hsp.convert_position_to_real_index(  end_pos  )

                if verbose :
                        print("Calculating scores using pssm...", start, '/', len(hsp.query), end, '/', len(hsp.hit))
                        print('query:', hsp.query[start:end])
                        print('hit  :', hsp.hit[start:end])

                gap_flag = 0

                #if the start position is gap opening, mark the gap_start
                if hsp.hit[start] == '-' or hsp.query[start] == '-' :
                        if ( start != 0 and hsp.hit[start-1] != '-' ) or ( start != 0 and hsp.query[start-1] != '-' ) :
                                gap_flag = 0
                        else :  
                                gap_flag = 1


                score = 0
                query_resnum = hsp.query_start_index - 1
                for q, h in zip(hsp.query[start:end],hsp.hit[start:end]) :
                        if q.isalpha() :
                                query_resnum += 1
                        if q == '-' or h == '-' :
                                if gap_flag :
                                        score -= self.gep
                                else :  
                                        score -= self.gop
                                        gap_flag = 1
                        else :  
                                score += score_mat( query_resnum, q, h )
                                gap_flag = 0

                        if verbose :
                                print('score:', score)


                return score



        def combine_optimal_by_pssm( self, hsp1, hsp2 ) :
                '''
                This function finds optimal intersection point between two intersecting hsps.
                '''
                new_hsp = copy.deepcopy( hsp1 )

                pass

        def blosum( self, resnum, res1, res2 ) :
                #resum is not used for blosum
                #place holoder for pssm score. :)
                score = 0.0
                try :   
                        score = self.blosum_mat[res1.upper(), res2.upper()]
                except :
                        print("WARNING:", res1.upper(), "or", res2.upper(), "not found in BLOSUM Mat.", file=sys.stderr)
                return score

        def set_blosum_mat( self, fn='' ) :
                mat_str = ''
                if fn : 
                        mat_str = open(fn).read()

                else : #default 
                        mat_str = '''#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
'''
                lines = mat_str.split('\n')

                aa_list = lines[6].split()

                for i, l in enumerate(lines[7:7+len(aa_list)]) :
                        l = l.split()
                        a = l[0]
                        scores = l[1:]
                        for j, score in enumerate(scores) :
                                self.blosum_mat[ a, aa_list[j] ] = float(score)

"""

def build_blast_db( input_file, db_name, formatdb='formatdb' ) :
	'''
	Run formatdb program to build a blast database.
	'''
	os.system( formatdb + ' -i ' + self.dummy_db + ' -f ' + db_name )


def build_pssm( msa, profile_file, **kwargs ) :
	'''
	Run PSIBLAST to get the position specific matrices.
	'''
	fp = tempfile.NamedTemporaryFile()
	temporary_db_name = fp.name

	build_blast_db( msa.query_file, temporary_db_name )
	output = temporary_db_name + ".psiblast_output"

	input_alignment = temporary_db_name + ".psiinputaln"
	msa.write_input_alignment( input_alignment )

	
	runner = PSIBLASTRunner( number_of_iterations=1, input=msa.query_file, input_alignment=input_alignment, output=output, profile_file=profile_file, database=temporary_db_name )
	runner.run()

	os.system( 'rm -f %s*'%temporary_db_name )

def read_pssm( profile_file ) :
	#print >>sys.stderr, '#Start to read profile', fn

	fp = open( profile_file )

	fp.readline() #read one line
	line = fp.readline()
	header = fp.readline().split()[:20]

	for i, a in enumerate(header) :
		self.residue_map[ a.upper() ] = i
		self.residue_map[ a.lower() ] = i

	for i, line in enumerate(fp) :
		if line == '\n' :
			break
		if verbose :
			print line

		line = line.split()

		index, residue = line[:2]
		if verbose :
			print 'index:', index, 'residue:', residue

		if i != int(index)-1 :
			print >>sys.stderr, "WARNING: ", fn, "has problem!", i, index

		self.columns.append( [ int(s) for s in line[2:22] ] )
		#vector = []
		#for s in line [2:22] :
			#vector.append( int(s) )
		#self.columns.append( vector )

	if self.length != len(self.columns) :
		print >>sys.stderr, 'WARNING! the length of query and profile does not match!'
"""


def build_neighboring_msas( input_fasta=None, neighboring_fastas=None, max_iterations_for_neighbors=None, echo=True, **kwargs ) :
	'''
	Run PSIBLAST to get the MSA for neighboring positions.
	'''

	if verbose :
		print("neighboring_fastas:", neighboring_fastas)
		print("input_fasta")
		print(input_fasta)
		print("input_fasta range:", input_fasta.sequencerange)
		
	if not input_fasta and not neighboring_fastas :
		return []

	elif input_fasta  :
		sequencerange = input_fasta.sequencerange
		if neighboring_fastas :
			neighboring_fastas.extend( input_fasta.extract_fastas( split_fragments=True, inverse=True ) )
		else :
			neighboring_fastas = input_fasta.extract_fastas( split_fragments=True, inverse=True )

	if verbose :
		print("neighboring_fastas:", neighboring_fastas)

	if not neighboring_fastas :
		return []

	neighboring_msas = []

	###############
	#remove all input and output settings possibly screw things up.
	#so remove them.!!
	if kwargs.get( 'input' ) :
		kwargs['input'] = None
	if kwargs.get( 'output' ) :
		kwargs['output'] = None
	if 'input_string' in kwargs :
		del kwargs['input_string']
	if kwargs.get( 'save_dir' ) :
		kwargs['save_dir'] = None
	if kwargs.get( 'range' ) :
		kwargs['range'] = None
	###############

	#################################
	#need to set up the arguments.
	if max_iterations_for_neighbors != None :
		kwargs['max_iterations' ] = max_iterations_for_neighbors
	else :
		kwargs['max_iterations' ] = default_max_iterations_for_neighbors
	#################################
		

	for i,fasta in enumerate(neighboring_fastas) :

		if echo :
			print('Building neighboring MSA', i ,'...') 

		namedfp = tempfile.NamedTemporaryFile()
		tempout = namedfp.name
		runner = PSIBLASTRunner( input_string=str(fasta), output=tempout, **kwargs)
		runner.run()
		
		blast = BLAST.BLAST( query=fasta, blast_result_fn=tempout )
		neighboring_msas.append( blast[-1] )

	return neighboring_msas



class BLASTOptionError( Exception ) :
	pass


if __name__ == '__main__' :
	import sys
	#r = BLASTRunner( input_string=open(sys.argv[1]).read(), output=sys.argv[2] )
	fasta = FASTA.FASTA()
	fasta.parse( sys.argv[1] )
	r = PSIBLASTRunner( input_string=str(fasta), output=sys.argv[2] )
	r.run()
	#print r.poll()
