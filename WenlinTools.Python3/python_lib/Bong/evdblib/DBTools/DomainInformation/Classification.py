'''
This module contains a class that manage protein classification
'''
from evdblib.DBTools import Settings
 
class Classification :
	'''
	Classification class has functionalities of 
	1) reading saved classification, 
	2) change single classification, and
	3) batch change the classification.
	'''

	def __init__( self, ) :
		#get setting class
		print(Settings.settings)

	def read( self, filename ) :
		pass
