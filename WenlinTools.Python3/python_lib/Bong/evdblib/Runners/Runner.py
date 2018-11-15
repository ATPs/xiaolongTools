import os, signal, shlex, time

from time import sleep
from subprocess import Popen, PIPE
from threading import Thread, Timer

module_verbose = 0

class ThreadedPopen(Thread) :
	'''
	ThreadedPopen is a subprocess opener wrapped by Thread class.
	This class is designed to be used in Runner class.
	'''
	def __init__( self, args_popen, timer=None, **kwargs_popen ) :
		Thread.__init__(self)

		self.args_popen = args_popen
		self.kwargs_popen = kwargs_popen
		self.pid = None
		self.subproc = None
		self.start_time = None #time.time()
		self.end_time = None

		self.timer=timer

	def kill(self) :
		if not self.subproc :
			return
		
		#check is needed since the python v2.5 does not have 
		#kill function
		if hasattr( self.subproc, "kill" ) :
			self.subproc.kill()
		#instead the kill function in os module should be used!
		#note that this os.kill is only available in unix/linux.
		else :
			os.kill( self.subproc.pid, signal.SIGKILL )


	def run(self) :
		if self.args_popen :

			if module_verbose : 
				print("starting...", ' '.join( self.args_popen ))

			self.start_time = time.time()
			self.subproc = Popen( self.args_popen, **self.kwargs_popen )
			self.pid = self.subproc.pid
			#self.stdout, self.stderr = self.subproc.communicate( input=input )
			self.end_time = time.time()

			#important to turn off the Timer class
			#to end the process
			if self.timer :
				self.timer.cancel() 

			if module_verbose : 
				print("finished...", ' '.join( self.args_popen ))
		else :
			if module_verbose : 
				print("No command line is given")


	def get_walltime( self ) :
		'''
		returns wall clock time to run the command.

		If the process is still running,
		it returns None.
		'''
		if self.start_time != None and self.end_time != None :
			return self.end_time - self.start_time
		

	def communicate( self, input=None ) :
		'''
		put input stream in string format
		and returns standard output data and standard error data
		in a tuple of strings (stdoutdata, stderrdata).
		'''
		if self.subproc :
			return self.subproc.communicate( input=input )
		else :
			if module_verbose :
				print("ThreadedPopen no subproc:", input) 
				print("was not delivered!")

			#somhow the start of subproc can be slower than communicate
			#when the communicate function is shooted right after the run function.
			#so this function need to wait until the subproc is available.
			stdin, stdout, stderr = self.kwargs_popen.get('stdin'), self.kwargs_popen.get('stdout'), self.kwargs_popen.get( 'stderr' )
			if stdin==PIPE or stdout==PIPE or stderr==PIPE  :
				while not self.subproc :
					sleep( 0.5 )
				return self.subproc.communicate( input=input )
			else :
				return None, None
					

	def poll( self ) :
		if self.subproc :
			return self.subproc.poll()

	def wait( self ) :
		if self.subproc :
			return self.subproc.wait()



class TimedRunner :
	'''
	A class to facilatate the timed running of the given command.
	
	Initializing is simple and need to just issue a command as in a string.
	If the runtime of the command need to be limited, just put the keyword argument
	timeout in seconds.

	Note that by default this class run the subprocess to fork independently.
	So, if you want to wait until the run finishes,
	wait() function should be called.

	example:
	runner = Runner( "ls -l" ) #this will make the class
	runner.run()

	runner2 = Runner( "sleep 100", timeout=10 ) #suppose to stops after 10sec.
	runner2.run()
	runner2.wait() #this function makes the main process waits until the runner2 run finishes or timout is reached..

	'''
	def __init__( self, command_line, timeout=0, **kwargs_popen ) :
		self.command_line = command_line
		self.kwargs_popen = kwargs_popen

		self.args_popen = shlex.split( command_line )
		if not self.args_popen :
			raise RunnerNullCommand()

		#setting timer thread
		if timeout > 0 :
			self.timer = Timer( timeout, self.kill )
		elif timeout < 0 :
			self.timer = None
			raise ValueError( "Timeout value should be positive integer or 0 (No timeout)." )
		else :
			self.timer = None

		self.popen = ThreadedPopen( self.args_popen, timer=self.timer, **kwargs_popen )

	def kill( self ) :
		if self.popen and self.popen.poll() == None :
			if module_verbose : print('Killing', self.popen.pid, '...')
			self.popen.kill()
		elif self.popen.poll() != None :
			if module_verbose : print('Already dead', self.popen.pid, '...')
			pass
		elif not self.popen :
			raise NotRunningError( '' )

	def wait( self ) :
		self.popen.join()
		

	def run( self ) :
		self.popen.start()
		if self.timer :
			self.timer.start()
		
	def communicate( self, input=None ) :
		'''
		delivers input string into the running process
		and returns file pointers of standard output and error.

		Note that stdin, stdout, stderr arguments should have been	
		given as PIPE to be used.
		'''
		if self.popen :
			try :
				return self.popen.communicate( input=input ) 
			except OSError :
				pass

		return None, None

	def get_walltime( self ) :
		'''
		returns wall clock time of the process running time in seconds.
		Note that it returns None when the process is still running.
		'''
		if self.popen :
			if module_verbose : self.command_line, self.popen.pid, 'walltime:', self.popen.get_walltime()
			return self.popen.get_walltime()

	def poll( self ) :
		if self.popen :
			return self.popen.poll()

if __name__ == '__main__' :
	#thread1 = MsgThread( 'thread1' )
	#thread2 = MsgThread( 'thread2' )
	'''
	thread1 = ThreadedPopen( ['sleep','100'], stdout=PIPE, stderr=PIPE ) 
	timer1 = Timer( 3, thread1.kill )
	thread2 = ThreadedPopen( ['sleep','90'] , stdout=PIPE, stderr=PIPE ) 
	timer2 = Timer( 3, thread2.kill )

	

	thread1.start()
	timer1.start()
	thread2.start()
	timer2.start()
	
	thread1.join()
	thread2.join()

	print thread1.communicate()
	print thread2.communicate()

	print 'thread1.poll()', thread1.poll()
	print 'thread2.poll()', thread2.poll()

	if thread1.isAlive()  :
		thread1.kill()
	if thread2.isAlive() :
		thread2.kill()
	#sleep(3)

	#thread1 = Runner( 'sleep 100', 4,  stdout=PIPE, stderr=PIPE ) 
	thread1 = TimedRunner( 'sleep 100', 4 ) 
	#thread2 = Runner( 'sleep 101', 4,  stdout=PIPE, stderr=PIPE ) 
	thread2 = TimedRunner( 'sleep 101', 4,  stdout=PIPE, stderr=PIPE ) 
	thread1.run()
	thread2.run()

	print "polling_thread1:", thread1.poll()
	print "polling_thread2:", thread2.poll()

	print "output from runner1"
	print thread1.communicate( )
	#print thread1.communicate( )[1]

	print "output from runner2"
	print thread2.communicate( )
	#print thread2.communicate( )[1]

	
	print 'walltime:',thread1.get_walltime()
	print 'walltime:',thread2.get_walltime()

	thread1.wait()
	thread2.wait()

	print 'walltime:',thread1.get_walltime()
	print 'walltime:',thread2.get_walltime()
	print "polling_thread1:", thread1.poll()
	print "polling_thread2:", thread2.poll()
	'''

	input_string = '''>12asB ASPARAGINE SYNTHETASE 11.83 2.2 X-RAY 1-30
AYIAKQRQISFVKSHFSRQLEERLGLIEVQ
	'''
	blastrun = TimedRunner( 'blastpgp -o 12asB.psi -e 0.001 -b 100 -v 30000 -I T -j 3', stdin=PIPE, stdout=PIPE, stderr=PIPE )
	blastrun.run()

	blastrun.communicate( input=input_string )
