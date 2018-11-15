'''
This module contains the SGEJobScript class that manages
a submitted job to a SGE Job queue
and the SGEJobQueue class that controls how many jobs are submitted
to the SGE queue.
'''
import os, tempfile, shutil, sys, shlex, io
from subprocess import Popen, PIPE
from time import sleep

verbose = 0

class NoRemoteHostAvailableError( Exception ) :
	pass

from evdblib.Utils.SimpleThreadingTools import ThreadPool
class ComputeNodePool ( ThreadPool ) :
	'''
	Pool of compute nodes that runs commands remotely from
	a queue.
	'''
	def __init__(self, available_hosts, low_limit=8 ) :

		#default is checking for low limit is 8.
		#If available hosts are lower than this number
		#then there is no need to run remotely.
		if len(available_hosts) < low_limit :
			raise NoRemoteHostAvailableError( 'available_hosts:', len(available_hosts), '<', low_limit )

		self.tasks = Queue( len(available_hosts) )
		for i in range( len(available_hosts) ) :
			WorkerNode( self.tasks, available_hosts[i] )

	def add_command( self, command_line, **kwargs ) :
		self.tasks.put( (command_line, kwargs) )


class ComputeNodes :
	def __init__(self, qhost='qhost' ) :
		'''
		A container class for host list
		and other informations extracted through qhost command.

		Note that compute node will not automatically update the host information.
		If this class is used as a monitoring tool for the remote hosts,
		periodical update is required!
		'''
		self.hostinfos = {}
		self.qhost = qhost

		if self.qhost :
			self.update()
		
	def __iter__ (self):
		return self.hostinfos.__iter__()

	def __contains__ (self, y) :
		return y in self.hostinfos

	def __len__(self) :
		return len( self.hostinfos )

	def get_ncpu( self, hostname ) :
		if hostname in self.hostinfos :
			return self.hostinfos[hostname].get('ncpu')

	def get_load( self, hostname ) :
		if hostname in self.hostinfos :
			return self.hostinfos[hostname].get('load')

	def get_hostname_list( self ) :
		return list(self.hostinfos)

	def update( self, strict=True ) :
		'''
		Update host information
		'''
		p = Popen( self.qhost, stdout=PIPE, stderr=PIPE )
		info_txt = p.communicate()[0]
		infos = {}
		for l in info_txt.split('\n') :
			if l.startswith( 'HOSTNAME' ) :
				continue
			elif l.startswith( '-' ) :
				continue
			elif l.startswith( 'global' ) :
				continue
			elif not l :
				continue

			#expected compute node name starts with 
			#compute-...
			if not l.startswith( 'compute' ) and strict :
				raise SunGridEngineError( l, "Node name supposed to start with 'compute'")
			
			line = l.split() 
			try :
				hostname, arch, ncpu, load, memtot, memuse, swapto, swapus = line
			except :
				raise SunGridEngineError( "Node info line format cannot be understood!" )

			if load == '-' :
				#node is down!!
				continue

			infos[hostname] = { 
					"hostname": hostname, 
					'arch': arch,
					'ncpu': ncpu,
					'load': load,
					'memtot': memtot,
					'memuse': memuse,
					'swapto': swapto,
					'swapus': swapus
			}

		self.hostinfos = infos

	def _remote_cmd_runner( self, remote_runner='ssh', hostname='', args=[], verbose=True, stdout_descriptor=sys.stdout, stderr_descriptor=sys.stderr, **kwargs ) :
		'''
		wrapper for remote run of the command.
		'''
		if remote_runner and hostname and args :
			pass
		else :
			return 0

		if verbose :
			print('running...', hostname, ' '.join( args ))

		p = Popen( [remote_runner, hostname] + args, **kwargs )
		stdout, stderr = p.communicate()

		stdout_descriptor.write( stdout )
		stderr_descriptor.write( stderr )

		return 1
	

	def run_remote_command( self, command_line, parallele=8, verbose=True ) :
		'''
		Run a command in the command line on each node.

		If parallele option is false, 
		it behaves much like the cluster-fork SGE tool in fork.
		Note that the command_line is a list of shell arguments 
		supposed to be used as is in Popen in suprocess module.

		This function is a easy to use tool for 
		regularizing or standardizing environments in nodes.
		'''
		if parallele :
			pool = ThreadPool( parallele )
			for hostname in sorted( self.hostinfos ) :
				pool.add_task( self._remote_cmd_runner, hostname=hostname, args=command_line, stdout=PIPE, stderr=PIPE, verbose=verbose )
			pool.wait_completion()

		else :
			for hostname in sorted( self.hostinfos ) :
				self._remote_cmd_runner( hostname=hostname, args=command_line, stdout=PIPE, stderr=PIPE, verbose=verbose )


	def run_remote_command_thread_saving( self, command_line, pool, verbose=True ) :
		'''
		Run a command in the command line on each node.
		This function is written for avoid thread error "can't start new thread".

		pool is a ThreadPool instance from evdblib.Utils.
		'''
		for hostname in sorted( self.hostinfos ) :
			pool.add_task( self._remote_cmd_runner, hostname=hostname, args=command_line, stdout=PIPE, stderr=PIPE, verbose=verbose )
		pool.wait_completion()


	def run_remote_command2( self, command_line, parallele=-1, verbose=True ) :
		'''
		Run a command in the command line on each node.
		Almost same as run_remote_cmmand except that
		this command returns two lists of stdout and stderr 
		from each node.
		'''

		stdouts = {}
		stderrs = {}

		if parallele == -1 :
			parallele = len( self.hostinfos )

		if parallele :
			pool = ThreadPool( parallele ) 
			for hostname in sorted( self.hostinfos ) :
				stdout = io.StringIO()
				stderr = io.StringIO()

				stdouts[hostname] = stdout
				stderrs[hostname] = stderr

				pool.add_task( self._remote_cmd_runner, hostname=hostname, args=command_line, stdout_descriptor=stdout, stderr_descriptor=stderr, stdout=PIPE, stderr=PIPE, verbose=verbose )
		else :
			for hostname in sorted( self.hostinfos ) :
				stdout = io.StringIO()
				stderr = io.StringIO()

				stdouts[hostname] = stdout
				stderrs[hostname] = stderr

				self._remote_cmd_runner( hostname=hostname, args=command_line, stdout_descriptor=stdout, stderr_descriptor=stderr, stdout=PIPE, stderr=PIPE )

		return stdouts, stderrs
				

	def sort_host_by_availability( self ) :
		self.update()
		host_list = []
		for hostname in self.get_hostname_list() :
			host_list.append( (int(self.get_ncpu(hostname)) - float(self.get_load(hostname)), hostname) )

		host_list.sort()
		host_list.reverse()

		return host_list

	def select_available_nodes( self, number=None, remove_overheaded=True ) :
		hosts = self.sort_host_by_availability()

		if remove_overheaded :
			available_hosts = [ h for a,h in hosts if a>=0.9 ]
		else :
			available_hosts = [ h for a,h in hosts ]
	
		if number and number > 0 :
			return available_hosts[:number]
		else :
			return available_hosts

	
	def run_commands_remotely( self, command_lines, parallele=-1, chunk_size=1, DoNotUseOverheadedNodes=True ) :
		'''
		Run command_lines remotely.
		command_lines is a list of strings where each string ia a command line to run
		remotely.

		Note that each command lines will be run on nodes round-robin style,
		until all commands are finished.

		By default parallele is -1 which means use all available nodes,
		and chunk size is 1. For large number of commands, the chunk size
		is better to be proportionally bigger to avoid too much communication
		overhead.

		WARNING: Due to implementation limitation, chunk_size option does not
		work and ignored!
		'''
		#select
		hosts = self.select_available_nodes( number=parallele, remove_overheaded=DoNotUseOverheadedNodes )

		pool = ComputeNodePool( hosts )
		for line in command_lines :
			pool.add_command( line, stdout=PIPE, stderr=PIPE )
		pool.wait_completion()


	"""
	#need to be implemented!
	def run_commands_remotely2( self, command_lines, parallele=-1, chunk_size=1, DonNoteUseOverheadedNodes=False ) :

		'''
		Almost same as run_commands_remotely.
		Except that this function returns two lists 
		of captured Standard Outputs and Standard Errors
		in the form of strings.
		'''
		hosts = self.select_available_nodes( number=parallele, remove_overheaded=DoNotUseOverheadedNodes )
		cStringIO.StringIO()
		pool = ComputeNodePool( hosts )
		for line in command_lines :
			pool.add_command( line, stdout=PIPE, stderr=PIPE )
		pool.wait_completion()


	"""

from evdblib.Utils import build_sequence_filename, build_profile_filename
from evdblib.DBTools.UpdateTools.QualityControl import check_profile_integrity
from evdblib.DBTools import Settings

class JobScriptBuilder :
	'''
	Abstract Job Script Writer class
	'''
	def __init__( self, cmd='', domain_informations=[], job_dir=None, sgejobqueue=None ) :
		self.script_files = []

		self.cmd=cmd

		if job_dir == None :
			raise JobScriptWriterError( "No job script directory is given." )

		if not os.path.exists( job_dir ) :
			os.makedirs( job_dir )
		self.job_dir = job_dir


		if not domain_informations :
			raise JobScriptWriterError( "No domain information given." )
		
		for dominfo in domain_informations :
			if verbose :
				print("adding the following domain...")
				print(dominfo)

			filename = self.build_file(dominfo)
			if not filename :
				continue
			self.add_script_file( filename )

	def __iter__( self ) :
		return self.script_files.__iter__()

	def __getitem__(self, i ) :
		return self.script_files[i]

	def build_file( self, domid ) :
		pass

	def add_script_file( self, filename ) :
		self.script_files.append( filename )
		

class ProfileJobScriptBuilder( JobScriptBuilder ) :
	'''
	This class write job scripts for profile generation jobs.
	This is a tailored class for profile job script writing.

	This also shows an implementation example of the JobScriptBuilder class.
	'''
	def build_file( self, dominfo ) :
		filename = os.path.join( self.job_dir, 'prfb%s.job'%dominfo['uniqueid'] )

		if not self.cmd :
			self.cmd = Settings.get( 'profile_builder' )
			
		command = self.cmd 
		if dominfo.get( 'profile_sequence_range' ) :
			command += ' -r ' + str(dominfo['profile_sequence_range'])

		if Settings.get( 'profile_type' ) :
			command += ' -m %s ' % Settings.get('profile_type' )

		if Settings.get( 'blast_db' ) :
			command += ' -d %s ' % Settings.get( 'blast_db' )

		command += ' %(profile_sequence_file)s %(domain_path)s' % dominfo
		
		if os.path.exists( filename ):
			raise JobScriptWriteError( "Profile job script %s already exists!"%filename )
		
		fp = open( filename, 'w' )
		print(command, file=fp)
		fp.close()

		dominfo['progress'] = 2

		return filename


class ProfileSearchJobScriptBuilder( JobScriptBuilder ) :
	'''
	This class write job scripts for profile search jobs.
	This is a tailored class for profile search job script writing.
	'''
	def __init__( self, alignment_method='', search_db='', search_db_size=None, iteration='', **kwargs ) :
		from evdblib.DBTools import Settings
		self.alignment_method = alignment_method
		self.search_db = search_db
		self.search_db_size = search_db_size
		self.iteration=int(iteration)

		if self.search_db_size == None :
			raise ValueError ( "The search DB size is required!" )

		if alignment_method == 'COMPASS' :
			self.profile_extention = Settings.get( "compass_suffix" )
		elif alignment_method == 'HHsearch' :
			self.profile_extention = Settings.get( "hhm_suffix" )
		else :
			if alignment_method in Settings.get( "profile_comparison_methods" ) :
				raise ValueError ( "The method is not implemented!" )
			else :
				raise ValueError ( "Unknown Profile comparison method!" )

		self.alignment_extention = Settings.get( "alignment_suffix" )

		JobScriptBuilder.__init__( self, **kwargs )


	def build_file( self, dominfo ) :
		from evdblib.DBTools import Settings
		filename = os.path.join( self.job_dir, 'prfs%s.job'%dominfo['uniqueid'] )

		if not self.cmd :
			self.cmd = Settings.get( 'profile_searcher' )
			
		command = self.cmd 
		domain_path = dominfo['domain_path']
		profile_search_method = self.alignment_method
		profile_search_queryid = dominfo['uniqueid']
		iteration = self.iteration 

		query_iteration = min( self.iteration, check_profile_integrity(dominfo) )

		profile_search_query = build_profile_filename( domain_path, dominfo['uniqueid']+'.prof', query_iteration, self.profile_extention ) 
		profile_search_db = self.search_db
		profile_search_db_size = self.search_db_size

		profile_search_output = build_sequence_filename( domain_path, dominfo['uniqueid'], self.alignment_extention )

		command += ' -q %(profile_search_queryid)s -j %(iteration)s -m %(profile_search_method)s -u -s %(profile_search_db_size)s -d %(profile_search_db)s %(profile_search_query)s %(profile_search_output)s' % locals()
		
		if os.path.exists( filename ):
			raise JobScriptWriteError( "Profile search job script %s already exists!"%filename )
		
		fp = open( filename, 'w' )
		print(command, file=fp)
		fp.close()

		return filename


class StructureSearchJobScriptBuilder( JobScriptBuilder ) :
	'''
	This class write job scripts for structure search jobs.
	This is a tailored class for structure search job script writing.
	'''
	def __init__( self, alignment_method='', search_db='', search_db_size=None, **kwargs ) :
		from evdblib.DBTools import Settings
		self.alignment_method = alignment_method
		self.search_db = search_db
		self.search_db_size = search_db_size
		if self.search_db_size == None :
			raise ValueError ( "The search DB size is required!" )

		if alignment_method == 'DaliLite' :
			self.data_file_extention = Settings.get( "dali_data_suffix" )
		elif alignment_method == 'FAST' or alignment_method == 'TMalign' :
			self.data_file_extention = Settings.get( "processed_ca_structure_suffix" )
		else :
			if alignment_method in Settings.get( "structure_comparison_methods" ).split() :
				raise ValueError ( "The method is not implemented!", alignment_method )
			else :
				raise ValueError ( "Unknown structure comparison method!", alignment_method )

		self.alignment_extention = Settings.get( "alignment_suffix" )

		JobScriptBuilder.__init__( self, **kwargs )


	def build_file( self, dominfo ) :
		from evdblib.DBTools import Settings
		filename = os.path.join( self.job_dir, 'strs%s.job'%dominfo['uniqueid'] )

		if not self.cmd :
			self.cmd = Settings.get( 'structure_searcher' )
			
		command = self.cmd 
		domain_path = dominfo['domain_path']
		search_method = self.alignment_method
		search_queryid = dominfo['uniqueid']

		search_query = build_sequence_filename( domain_path, dominfo['uniqueid'], self.data_file_extention ) 
		search_db = self.search_db
		search_db_size = self.search_db_size

		search_output = build_sequence_filename( domain_path, dominfo['uniqueid'], self.alignment_extention )

		command += ' -q %(search_queryid)s -m %(search_method)s -u -s %(search_db_size)s -d %(search_db)s %(search_query)s %(search_output)s' % locals()
		
		if os.path.exists( filename ):
			raise JobScriptWriteError( "Structure search job script %s already exists!"%filename )
		
		fp = open( filename, 'w' )
		print(command, file=fp)
		fp.close()

		return filename




class SGEJobScript :
	def __init__( self, script_file=None, name=None, submit_cmd='qsub', delete_cmd='qdel' , use_script_dir = True) :

		'''
		This class is encapsulating job queue submission processes.
		If use_script_dir option is true,
		then the submission will be done 
		in the directory where the job script is!
		'''
		self.submit_cmd = submit_cmd 
		self.delete_cmd = delete_cmd

		self.job_id = None
		self.name = name

		self.script_file = script_file
		if not self.script_file :
			raise SGEJobScriptError( "No script file is given!" )
		elif not os.path.exists( self.script_file ) :
			raise SGEJobScriptError( "Script file does not exist!" )

		self.old_job_ids = []

		self.use_script_dir = use_script_dir

	def is_submitted( self ) :
		if self.job_id == None :
			return False

		#job_id = int(self.job_id)
		return True

	def submit( self ) :
		'''
		Submit a job script to the queue and get the SGE Queue Job ID.
		'''
		command_line = [ self.submit_cmd ]
		if self.name :
			command_line.append( '-N' )
			command_line.append( self.name )

		if Settings.get( 'sge_job_queue' ) :
			command_line.append( '-q' )
			command_line.append( Settings.get( 'sge_job_queue' ) )

		command_line.append( "-cwd" ) #run in the current directory
		command_line.append( self.script_file )

		#sorting out submit directory issue 
		#if the submission directory is not the same directory
		#of the script.
		if self.use_script_dir :
			script_dir = os.path.dirname( self.script_file )
		else :
			script_dir = None

		output = Popen( command_line, stdout=PIPE, cwd=script_dir ).communicate()[0]

		#saving previous job ids.
		if self.job_id :
			self.old_job_ids.append( self.job_id )

		try :
			self.job_id = output.split()[2]
		except IndexError :
			raise SunGridEngineError( "Submission %s failed!"%self.script_file )


	def delete( self ) :
		if self.job_id :
			popen = Popen( [self.delete_cmd, self.job_id] )
			self.old_job_ids.append( self.job_id )
			self.job_id = None
			return popen.wait()
		

class SGEJobQueue :
	'''
	'''
	def __init__( self, check_interval=120 ) :
		self.jobscripts = []
		self.check_interval=60
		self.qstat = SGEqstat()

	def __iter__( self ) :
		return self.jobscripts.__iter__()

	def __getitem__( self, i ) :
		return self.jobscripts[i]

	def job_ids2job_list( self, job_ids ) :
		jobids2job = {}
		for job in self :
			if job.job_id :
				jobids2job[job.job_id] = job
		
		jobs = []
		for job_id in job_ids :
			if job_id in jobids2job :
				jobs.append( jobids2job[job_id] )
		return jobs

	def add( self, scriptfile, submit=True ) :
		'''
		Adds a script file into the queue.
		By default the jobs is submtted.
		'''
		sgejobscript = SGEJobScript( script_file=scriptfile )
		self.jobscripts.append( sgejobscript )
		if submit == True :
			sgejobscript.submit()

	def get_submitted_job_ids( self ) :
		jobs = []
		for job in self :
			if job.is_submitted() :
				jobs.append( job.job_id )
		return jobs

	def get_currently_running_job_ids( self ) :
		submitted_ids = set( self.get_submitted_job_ids() )

		self.qstat.update_queue()
		job_ids_in_queue = set()
		job_ids_in_queue.update( self.qstat[ ['running', 'waiting'] ] )
		#job_ids_in_queue.update( self.qstat['waiting'] )

		return submitted_ids.intersection( job_ids_in_queue )

	def get_errored_job_ids( self ) :
		submitted_ids = set( self.get_submitted_job_ids() )

		self.qstat.update_queue()
		job_ids_in_queue = set()
		job_ids_in_queue.update( self.qstat['error'] )

		return submitted_ids.intersection( job_ids_in_queue )


	def wait( self ) :
		'''
		Waits until all jobs in the system finishes.
		Note that it returns list of jobscript objects that had problems in running.
		'''
		#get jobs not submitted.
		#and submit
		for job in self :
			if not job.is_submitted() :
				job.submit()
		
		#get jobs ids in the queue
		running_job_ids = self.get_currently_running_job_ids()
		while running_job_ids :
			sleep( self.check_interval )
			running_job_ids = self.get_currently_running_job_ids()

		erroneous_job_ids = self.get_errored_job_ids()
		return self.job_ids2job_list( erroneous_job_ids )

		
	def kill( self ) :
		'''
		Delete all running jobs.
		'''
		for jobscript in self :
			jobscrtip.delete()

class JobScriptWriteError( Exception ) :
	pass

class SCGJobScriptError( Exception ) :
	pass



class SGEqstat :
	def __init__( self, qstat_cmd='qstat' ) :
		self.qstat_cmd = qstat_cmd
		self.running_job_ids = []
		self.error_job_ids = []
		self.waiting_job_ids = []

	def __getitem__( self, status ) :
		job_ids = []
		if 'running' in status :
			job_ids.extend( self.running_job_ids )

		if 'waiting' in status :
			job_ids.extend( self.waiting_job_ids )

		if 'error' in status :
			job_ids.extend( self.error_job_ids )

		return job_ids

	def update_queue( self ) :
		'''
		returns list of job IDs in dictionary of list
		'''
		raw_output = Popen( [self.qstat_cmd], stdout=PIPE ).communicate()[0]
		job_ids = {'errors':[], 'running':[], 'waiting':[]}
		for line in raw_output.split('\n')[2:] :
			cols = line.split()
			if not cols :
				continue
			#debug
			#print cols

			job_id = cols[0]
			job_name = cols[2]
			owner = cols[3]
			job_state = cols[4]

			if 'E' in job_state or 's' in job_state or 'S' in job_state or 'T' in job_state :
				job_ids[ 'errors' ].append( job_id )

			elif 'w' in job_state or 'h' in job_state :
				job_ids['waiting'].append( job_id )

			elif 'r' in job_state or 't' in job_state :
				job_ids[ 'running' ].append(job_id)

		self.error_job_ids = job_ids['errors']
		self.running_job_ids = job_ids['running']
		self.waiting_job_ids = job_ids['waiting']
	

	def check_jobs_in_queue( self, job_id_list=None, stat_list=['errors', 'waiting', 'running'] ) :
		'''
		returns job ids that is not in the job queue in designated in the stat_list.
		'''
		if not job_id_list :
			raise SunGridEngineError( "Jod ID list should contain valid Jod IDs!" )
		job_ids = self.update_queue()

		job_ids_in_queue = set()
		for status in stat_list :
			job_ids_in_queue.update( job_ids[status] )

		job_id_set = set(job_id_list)
		return list(job_id_set - job_ids_in_queue)


class HostnameNotFoundError( Exception ) :
	pass

class RemotepathNotFoundError( Exception ) :
	pass

class RemoteCopyError( Exception ) :
	pass

def print_rmtree_error( func, path, exec_info ) :
	print('onerror', file=sys.stderr)
	print(func, file=sys.stderr)
	print(path, file=sys.stderr)
	print(exec_info, file=sys.stderr)

class RemoteFileCopier :
	def __init__( self, remotepath2hostname={}, save_dir=None, compute_nodes=None, copy_cmd='scp -q', remotepath2hostname_mapping_file=None ) :
		'''
		This class encapsulate remote file copying procedure.
		This class is supposed to use filename2hostname dictionary 
		if it is available. If not, just rely on the user's input.

		Note that the local copy of the files will be deleted when
		the class is destroyed, unless save_dir is specified. 
		'''

		self.copy_cmd = copy_cmd
		self.temp_dir = None
		self.save_dir = save_dir
		self.compute_nodes = compute_nodes #computename id iterator

		if not self.save_dir :
			self.temp_dir = tempfile.mkdtemp()
			self.save_dir = self.temp_dir

		self.remotepath2hostname = remotepath2hostname
		self.basename2hostname = {}
		self.basename2remotepath = {}

		if self.remotepath2hostname :
			for path, hostname in self.remotepath2hostname.items() :
				if self.compute_nodes and not hostname in self.compute_nodes :
					raise HostNameNotFoundError( hostname )

				self.basename2hostname[ os.path.basename(path) ] = hostname
				self.basename2remotepath[ os.path.basename(path) ] = path

		self.remotepath2hostname_mapping_file = remotepath2hostname_mapping_file
		self.add_remote_path_to_host_mapping_file( self.remotepath2hostname_mapping_file )


	def add_remote_path_to_host_mapping_file( self, remotepath2hostname_mapping_file ) :
		#adding additional mapping in the file.
		#the file is simply structured by two column lines
		#of hostname and remotepath.
		if remotepath2hostname_mapping_file :
			fp = open( remotepath2hostname_mapping_file )
			for l in fp :
				h, r  = l.split()
				self.remotepath2hostname[ r ] = h

		self.update_path_related_stuff()


	def update_path_related_stuff( self ) :
		if self.remotepath2hostname :
			for path, hostname in self.remotepath2hostname.items() :
				if self.compute_nodes and not hostname in self.compute_nodes :
					raise HostNameNotFoundError( hostname )

				self.basename2hostname[ os.path.basename(path) ] = hostname
				self.basename2remotepath[ os.path.basename(path) ] = path

	def __del__( self ) :
		'''
		destroys local copy directory, unless
		save_dir was given in the initializing process.
		'''

		if hasattr(self, 'temp_dir' ) and self.temp_dir :
			shutil.rmtree( self.temp_dir, ignore_errors=True  )
		'''
	rmtree(path, ignore_errors=False, onerror=None)
        Recursively delete a directory tree.
        
        If ignore_errors is set, errors are ignored; otherwise, if onerror
        is set, it is called to handle the error with arguments (func,
        path, exc_info) where func is os.listdir, os.remove, or os.rmdir;
        path is the argument to that function that caused it to fail; and
        exc_info is a tuple returned by sys.exc_info().  If ignore_errors
        is false and onerror is None, an exception is raised.
		'''

	def is_file_available( self, filename ) :
		return filename in self.basename2hostname or filename in self.remotepath2hostname
	

	def _run_copy_cmd( self, cmd, verbose_flag ) :
		if verbose_flag :
                        p = Popen( shlex.split( cmd ), stdout=PIPE, stderr=PIPE )
                else : 
                        p = Popen( shlex.split( cmd ) )

                stdout, stderr = p.communicate()
                if stdout :
                        sys.stdout.write( stdout )
                if stderr :
                        sys.stderr.write( stderr )


	def copy_file( self, remotePathOrFile, remotehost=None, verbose=False, retry = 3, retry_interval=3 ) :
		'''
		This function copies the given remote path or filename
		into the save_dir or (if not temporary directory),
		and returns the path of local file.
		'''
		
		if not remotehost :
			#if no remote host information is given
			#this function will find it..
			if remotePathOrFile in self.remotepath2hostname :
				remotepath = remotePathOrFile
				remotehost = self.remotepath2hostname[ remotePathOrFile ]
			elif remotePathOrFile in self.basename2hostname :
				remotehost = self.basename2hostname[ remotePathOrFile ]
				remotepath = self.basename2remotepath[ remotePathOrFile ]
			else :
				raise RemotepathNotFoundError( remotePathOrFile )

		local_path = os.path.join( self.save_dir, os.path.basename(remotepath) )
		infos = {'copy_cmd':self.copy_cmd, 'hostname': remotehost, 'remotepath':remotepath, 'localpath':local_path}
		cmd = "%(copy_cmd)s %(hostname)s:%(remotepath)s %(localpath)s" % infos

		self._run_copy_cmd( cmd, verbose )

		while (not os.path.exists(local_path)) and retry > 0 :
			retry -= 1
			sleep(retry_interval) #invertal
			self._run_copy_cmd( cmd, verbose )

		if os.path.exists( local_path ) :
			return local_path
		else :
			raise RemoteCopyError( remotehost, remotePathOrFile, local_path )


	def delete_file( self, remotePathOrFile, verbose=False ) :
		'''
		This function removes the given remote path or filename
		from the save_dir or (temporary directory),
		returns True of the file is in the directory,
		otherwise returns False.
		'''

		basename = os.path.basename(remotePathOrFile)
		delete_fn = os.path.join( self.save_dir, basename )
		if os.path.exists( delete_fn ) :
			os.remove( delete_fn )
			return True
		else :
			return False

	#def copy_files needs to be implemented!
	
			
from queue import Queue
from evdblib.Utils.SimpleThreadingTools import Worker
class WorkerNode( Worker ) :
	'''
	Executing command over the network!
	'''
	def __init__(self, tasks, hostname ) :
		self.hostname = hostname
		Worker.__init__(self, tasks )

	def run( self ) :
		while True :
			command_line, kwargs = self.tasks.get()
			try :
				self._remote_cmd_runner( command_line=command_line, **kwargs)
			except Exception as e: print(e)
			self.tasks.task_done()

	def _remote_cmd_runner( self, remote_runner='ssh', command_line='', **kwargs ) :
		'''
		wrapper for remote run of the command.
		'''
		print('running...', self.hostname, command_line.strip())

		p = Popen( [remote_runner, self.hostname, command_line], **kwargs )
		stdout, stderr = p.communicate()

		if stdout :
			sys.stdout.write( stdout )
		if stderr :
			sys.stderr.write( stderr )

		return 1


class RemoteFileManager :
	def __init__( self, compute_nodes=None ) :
		self.compute_nodes = compute_nodes
		if not self.compute_nodes :
			self.compute_nodes = ComputeNodes()

	def get_remote_mapping( self, remote_dir ) :
		'''
		Convenient fuction to prepare remote path to host mapping.
		'''
		stdouts, stderrs = self.compute_nodes.run_remote_command2( ["ls", '-1', remote_dir] )
		
		r2h={}

		for hostname, stdout in stdouts.items() :
			#stdout should be an instance of StringIO.StringIO
			#or at least an file like object!
			files = stdout.getvalue().split()
			for fn in files :
				r2h[fn] = hostname

		return r2h


class SunGridEngineError( Exception )  :
	pass
class JobScriptWriterError( Exception )  :
	pass

if __name__ == '__main__' :
	print("An examplary usage of the module")

	print("building a job")
