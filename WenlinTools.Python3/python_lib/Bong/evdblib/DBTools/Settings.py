#This module has all settings about the EvolvableDatabase.
import os, sys
import configparser

settings = {} #container for all configurations
db_home_dir =  os.path.abspath(
		os.path.join( 
			os.path.dirname( __file__ ),
			'..', '..' 
	      		)
		)

config_filename = 'db.config'

# use local config first
config_file = os.path.join( os.getcwd(), 'db.config' ) 
if not os.path.exists( config_file ) :
	config_file = os.path.join( db_home_dir , 'db.config' )

#print config_file

def get( property, silent=False ) :
	'''
	Reads settings dictionary where all the configuration is stored
	and returns the setting value.
	
	If property is not in the dictionary KeyError exception will occur,
	unless the silent flag is True

	Silent flag makes the function return Null String when the property is not
	defined.
	'''
	if silent :
		if property in settings :
			return settings[property]
		else :
			return ''
	else :
		return settings[property]

def set( property, value, silent=True ) :
	'''
	Set settings directionary value.

	If silent value is True, 
	overwrite already existing property value is allowed!
	'''
	if not silent :
		if property in settings :
			raise PropertyOverwritingError()

	settings[property] = value


import time
from evdblib.Utils import string2pathname
def prepare_intermediate_result( update_name='' ) :
	'''
	makes intermediate result directories and generate appropriate settings.

	The given update_name will be used to identify previous updates.
	Note that if the update_name is not given or same as previous update, this function
	will raise a UpdateNameError Exception.

	If the update name is not given, mm.dd.yyyy format of update name will be used.
	'''
	if not update_name :
		update_name = time.strftime('%m.%d.%Y')

	intermediate_result_root_dir = get( 'intermediate_result_root_dir' )
	if os.path.exists( intermediate_result_root_dir )  :
		pass
	else :
		os.makedirs( intermediate_result_root_dir ) 

	##############################################
	#checking if the root directory exist or not!!
	##############################################
	if not os.path.exists( intermediate_result_root_dir ) :
		raise UpdateNameError( "No intermediate result root directory can be prepared! %s"%
			intermediate_result_root_dir )

	#################################################
	#checkinf if the update name contains 
	#not allowed characters
	#################################################
	if update_name != string2pathname( update_name ) :
		raise UpdateNameError( "Name of an update better not contain special characters since the name is used in building directory names. Use this instead %s."% string2pathname( update_name ) )

	intermediate_result_dir = os.path.join( intermediate_result_root_dir, update_name )

	#################################################
	#add the add hoc setting value "intermediate_result_dir".
	#################################################
	settings[ 'intermediate_result_dir' ] = intermediate_result_dir
	settings[ 'intermediate_domain_info_db' ] = os.path.join( intermediate_result_dir, 'domain_info.db' ) 
	settings[ 'intermediate_classification_info_db' ] = os.path.join( intermediate_result_dir, 'classification_info.db' ) 
	
	if os.path.exists(intermediate_result_dir) :
		raise UpdateNameError( "Update name already exists! Try another name!" )
	else :
		os.makedirs( intermediate_result_dir )
		if not os.path.exists( intermediate_result_dir ) :
			raise UpdateNameError( "Intermediate result directory cannot be prepared.", intermediate_result_dir )


class UpdateNameError( Exception ) :
	'''
	This error will be raised when the update name has some problems most likely
	the given update name already exists.
	'''
	pass


class SettingsManager:
	'''
	Manage all settings from the config file.
	'''

	def __init__( self ) :
		self.config_file = config_file

		#print self.config_file
		if not os.path.exists( self.config_file ) :
			print("No config file is found!", self.config_file, file=sys.stderr)
			print("Run configuration initilize program to set up directories.", file=sys.stderr)


	def question( self, question, current_value='' ) :
		'''
		Asks a question and read the user response and
		return the user response as String.

		#This funciton is written for configure_initial_settings function.
		'''
		print(question)
		print("Current:", current_value)
		new_value = input( "New:" )

		return new_value
				

	def configure_initial_settings( self ) :
		'''
		Reads initial settings interactively from user.
		'''
		print("Starting initial setting of Evolvable Database...\n")
		print("This function will guide you to set up all necessary conditions.")
		
		pass
		print(settings)

		#save the setting data
		print("Setting is done!")
		

	def read( self, config_file='' ) :
		'''
		Reads settings from the config file.
		'''
		if not config_file :
			config_file = self.config_file

		parser = configparser.ConfigParser()
		parser.read( config_file )
		parser.set( 'settings', 'database_dir', db_home_dir )

		settings.update( parser.items( 'settings' ) )

	def save( self, config_file='' ) :
		pass


###############################
#reading configuration file
###############################
manager = SettingsManager()
manager.read()
