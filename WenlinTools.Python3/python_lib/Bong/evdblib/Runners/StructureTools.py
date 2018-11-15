import os, sys, tempfile, shutil

from evdblib.Utils import parse_sequence_filename
from evdblib.DBTools import Settings
from subprocess import Popen, PIPE
from evdblib.Utils.Parsers.PDB.PDB import PDB

class DaliLiteDATRunnerError( Exception ) :
	pass


class DaliLiteDATRunner :
	def __init__( self, pdbfn, cmd=None, 
			save_dir=None, 
			output_fn=None, 
			identifier=None, 
			echo=True ) :
		'''
		Run DaliLite command and generate DAT file
		for large scale DaliLite structure comparisons.

		identifier is quite similar to PDBID (with chainID followed) 
		except that DaliLite differentiate lower vs upper characters and 
		the first character does not have to be an numeric character.

		Note that the settings in this DAT generator is somewhat different
		from the original settings in DaliLite program.
		In DaliLite the filename of DAT file should be (I guess)
		matching with the identifier and chainID.
		
		echo is a flag that controls output from the program
		for logging purposes and the actual command run to generate
		the DAT file.

		WARNING: This class depends on a slightly modified version
		of DaliLite that saves DAT files into the DAT directory in 
		the current working directory.
		'''

		self.cmd = cmd
		self.pdbfn = pdbfn
		self.save_dir = save_dir
		self.output_fn = output_fn
		self.identifier = identifier
		self.echo = echo

		self.temp_dir = tempfile.mkdtemp()

		if self.cmd == None :
			self.cmd = Settings.get( "dalilite" )

		if self.save_dir :
			if not os.path.exists( self.save_dir ) :
				os.makedirs( self.save_dir )
		else :
			self.save_dir = os.getcwd()

		#convert the data into absolute_dir.
		if self.pdbfn :
			shutil.copy( self.pdbfn, self.temp_dir )
			self.pdbfn = os.path.join( self.temp_dir, os.path.basename(self.pdbfn) )

		if not self.output_fn :
			dir, basename, ext = parse_sequence_filename( os.path.basename( self.pdbfn ) )
			self.output_fn = basename+'.dat'

		if os.path.exists( os.path.join(self.save_dir, self.output_fn)  ) :
			raise IOError( "File already exists.", self.output_fn )

		#read settings file to get the default idenfier
		#this is important for internal consistency
		if self.identifier == None :
			self.identifier = Settings.get( "default_dali_id" )

		#original DAT file that will be produced by the DaliLite program!
		self.dat_fn = os.path.join( self.temp_dir, "DAT", self.identifier + ".dat" )
		self.dssp_fn =  os.path.join( self.temp_dir, self.identifier[:4] + ".dssp" )

	def __del__( self ) :
		if os.path.exists( self.temp_dir ) :
			shutil.rmtree( self.temp_dir )

	def get_command_line( self ) :
		return [ self.cmd, "-readbrk", self.pdbfn, self.identifier[:4] ]
		
	def run( self ) :
		cmdline = self.get_command_line()
		if self.echo :
			print(cmdline)

		proc = Popen( cmdline, stdout=PIPE, stderr=PIPE, cwd=self.temp_dir ) 

		stdout, stderr = proc.communicate()
		if self.echo :
			sys.stdout.write( stdout )
			sys.stderr.write( stderr )

		if not os.path.exists(self.dat_fn) :
			raise DaliLiteDATRunnerError( self.pdbfn, self.identifier, "DAT file generation failed!" )

		if not os.path.exists( self.dssp_fn ) :
			raise DaliLiteDATRunnerError( self.pdbfn, "DSSP file generation failed!" )

		shutil.copy( self.dat_fn, os.path.join(self.save_dir, self.output_fn) )

class MaxSproutRunnerError( Exception ) :
	'''
	This Exception will be raise when
	the code running is failed or input/output parameter is not set correctly!
	'''
	pass

class MaxSproutRunnerSequenceChangeError( Exception ) :
	'''
	This exception will be raise when
	the input and output PDB has different sequences!
	'''
	pass


class MaxSproutRunner :
	def __init__( self, inputpdb=None, 
			read_cmd=None, buildbackbone_cmd=None, dglp_list=None,
			save_dir=None, output_fn=None, 
			echo=True ) :
		'''
		Run MaxSprout command buildbackbone to generate a new PDB file
		filling backbone atoms for CA only residues.
		This class does not full model building but only builds backbones.
		Full sidechain optimization has some problems.

		echo is a flag that controls output from the program
		for logging purposes and the actual command run to generate
		the full backbone PDB file.

		Note that inputpdb is assumed to be a single chain PDB file.
		If the input PDB sequence and output PDB sequence does not match,
		MaxSproutRunnerSequenceChangeError will be raised!

		Note that the finally backbone built model is actually not
		a complete PDB file format.
		It has missing records like occupancies and B-factors.
		'''
		self.read_cmd = read_cmd
		self.buildbackbone_cmd = buildbackbone_cmd
		self.dglp_list = dglp_list #Necessary input param for buildbackbone cmd.

		self.pdbfn = inputpdb
		self.save_dir = save_dir
		self.output_fn = output_fn
		self.echo = echo

		self.temp_dir = tempfile.mkdtemp()

		if self.read_cmd == None :
			self.read_cmd = Settings.get( "maxsprout_readbrk" )
		if self.buildbackbone_cmd == None :
			self.buildbackbone_cmd = Settings.get( "maxsprout_buildbackbone" )
		if self.dglp_list == None :
			self.dglp_list = Settings.get( "maxsprout_dglp_list" )

		if not (self.read_cmd and self.buildbackbone_cmd and self.dglp_list) :

			for k,v in Settings.settings.items() :
				print(k, ":", v)

			raise MaxSproutRunnerError( "Commands or dglp.list info was not retrieved!", self.read_cmd, self.buildbackbone_cmd, self.dglp_list)

		if self.save_dir :
			if not os.path.exists( self.save_dir ) :
				os.makedirs( self.save_dir )
		else :
			self.save_dir = os.getcwd()

		#convert the data into absolute_dir.
		if self.pdbfn :
			shutil.copy( self.pdbfn, self.temp_dir )
			self.pdbfn = os.path.join( self.temp_dir, os.path.basename(self.pdbfn) )

		if not self.output_fn :
			dir, basename, ext = parse_sequence_filename( os.path.basename( self.pdbfn ) )
			self.output_fn = basename+'.maxsprout'

		if os.path.exists( os.path.join(self.save_dir, self.output_fn)  ) :
			raise IOError( "Maxsprouted file already exists.", self.output_fn )

	def __del__( self ) :
		if os.path.exists( self.temp_dir ) :
			shutil.rmtree( self.temp_dir )

	def get_read_command_line( self ) :
		read = [self.read_cmd, "-pdb", os.path.basename(self.pdbfn), "-rd", self.temp_dir, "-wd", self.temp_dir]
		return read

	def get_buildbackbone_command_line( self ) :
		pdb_basename = parse_sequence_filename(self.pdbfn)[1] #basename
		buildbackbone = [ self.buildbackbone_cmd, "-pdb", pdb_basename, "-pl", self.dglp_list, "-d", "Y" ]

		return buildbackbone


	def copy_result( self, source, target ) :
		'''
		The output of buildbackbone command in maxsprout
		is not a complete PDB file.
		It has missing Atom numbers.

		This function fill in the atom number by just putting simple
		continuous number from 1.
		'''
		content = open( source ).readlines()
		#sys.stdout.writelines( content )
		new_content = []
		anum = 1
		for l in content :
			if l[:6] == 'ATOM  ' :
				nl = l[:6] + "%5d"%anum + l[11:]
				anum += 1
				new_content.append( nl )

		fp = open(target, 'w' )
		fp.writelines( new_content )
		fp.close()

 
	def check_built_model( self, inputpdb=None, outputfile=None ) :
		'''
		raise exceptions if the results are not correct!
		'''
		if not inputpdb :
			inputpdb = self.pdbfn

		if not outputfile :
			outputfile = self.output_fn

		inputpdb = PDB( inputpdb )
		if not inputpdb :
			raise MaxSproutRunnerError( "input file is not in correct PDB format", inputpdb )
		if len(inputpdb.get_model().chain_ids) != 1 :
			raise MaxSproutRunnerError( "More than one chain in the input PDB", self.pdbfn )
		inputchainid = inputpdb.get_model().chain_ids[0]

		outputpdb = PDB( outputfile )
		if not outputpdb :
			raise MaxSproutRunnerError( "Output file is not in correct PDB format", outputfile )
		if len(outputpdb.get_model().chain_ids) != 1 :
			raise MaxSproutRunnerError( "More than one chain in the output PDB", outputfile )
		outputchainid = outputpdb.get_model().chain_ids[0]

		inputchain = inputpdb[inputchainid]
		outputchain = outputpdb[outputchainid]

		inputseq = inputchain.extract_sequence()
		outputseq = outputchain.extract_sequence()

		if inputseq != outputseq :
			print(inputseq)
			print(outputseq)
			raise MaxSproutRunnerSequenceChangeError( self.pdbfn, outputfile )

		
	def run( self, check_output=True ) :
		'''
		Run two commands
		'''
		##############
		#reading pdb file
		##############
		cmdline = self.get_read_command_line()
		if self.echo :
			print(cmdline)

		proc = Popen( cmdline, stdout=PIPE, stderr=PIPE, cwd=self.temp_dir ) 

		stdout, stderr = proc.communicate()
		if self.echo :
			sys.stdout.write( stdout )
			sys.stderr.write( stderr )

		returncode = proc.poll()
		if returncode : #non-zero return code meaning something is not right!
			print(stderr, file=sys.stderr)
			raise MaxSproutRunnerError( "Failed to read PDB file", self.pdbfn, returncode )

		##############
		#building backbone
		##############
		cmdline = self.get_buildbackbone_command_line()
		if self.echo :
			print(cmdline)

		proc = Popen( cmdline, stdout=PIPE, stderr=PIPE, cwd=self.temp_dir ) 

		returncode = proc.poll()
		if returncode : #non-zero return code meaning something is not right!
			print(stderr, file=sys.stderr)
			raise MaxSproutRunnerError( "Failed to read PDB file", self.pdbfn, returncode )

		stdout, stderr = proc.communicate()
		if self.echo :
			sys.stdout.write( stdout )
			sys.stderr.write( stderr )

		result_fn = parse_sequence_filename(self.pdbfn)[1] + ".brk_mod" #as defined in the buildbackbone program

		#copy the result file in the temporary directory
		#into the save_dir with the filename defined as output_fn.
		self.copy_result( os.path.join(self.temp_dir, result_fn), os.path.join(self.save_dir, self.output_fn) )

		if check_output :
			self.check_built_model()


#############################
#Legacy code to run maxsprout
#############################
def build_backbone( pdb_content, id='junk' ) :
        temp_dir = tempfile.mkdtemp()
        cur_dir = os.getcwd()

        os.chdir( temp_dir )

        pdb_fn_base = id
        pdb_fn = id + '.pdb'
        open( pdb_fn, 'w' ).writelines( pdb_content ) #save pdb_content as junk.pdb
        
        maxsprout_read_cmd = "/home/kim/local/maxsprout_new/readbrk -pdb %s -rd ./ -wd ./ >& /dev/null"
        maxsprout_buildbackbone_cmd = '/home/kim/local/maxsprout_new/buildbackbone -pdb %s -pl /home/kim/local/maxsprout_new/dglp.list -d Y >& /dev/null'

        os.system( maxsprout_read_cmd % pdb_fn )
        os.system( maxsprout_buildbackbone_cmd % pdb_fn_base )

        content = open( pdb_fn_base + ".brk_mod" ).readlines()
        #sys.stdout.writelines( content )
        new_content = []
        anum = 1
        for l in content :
                if l[:6] == 'ATOM  ' :
                        nl = l[:6] + "%5d"%anum + l[11:]
                        anum += 1
                        new_content.append( nl )

        os.chdir( cur_dir )
        shutil.rmtree( temp_dir,1 )

        return new_content

