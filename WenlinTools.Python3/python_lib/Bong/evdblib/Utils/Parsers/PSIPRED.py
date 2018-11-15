import os

verbose = 1

class PSIPREDParseError( Exception ) :
	pass

class PSIPRED :
	def __init__( self, horiz=None, query=None ) :
		'''
		Parser for PSIPRED result .horiz file.

		Optional query fasta object can be also given for 
		possible parsing error check!
		'''

		self.horiz = horiz
		self.query = query

		self.conf_marker = 'Conf:'
		self.pred_marker = 'Pred:'
		self.aa_marker   = '  AA:'

		if self.horiz :
			self.parse_horiz()


	def get_confidence( self ) :
		return self.confidence

	def get_prediction( self ) :
		return self.prediction

	def parse_horiz( self, horiz=None, query=None ) :
		'''
		Parse PSIPRED result .horiz file.
		'''
		if horiz == None :
			horiz = self.horiz

		if query == None :
			query = self.query

		if not horiz :
			return 

		if not os.path.exists( horiz )  :
			raise PSIPREDParseError( "PSIPRED result (horiz) file not given!" )

		if verbose :
			print("result file:", horiz)
			print("query sequence", query.sequence)

		conf = ''
		pred = ''
		aa = ''
		
		fp = open( horiz )
		for line in fp :
			if line.startswith( self.conf_marker ) :
				line = line[len(self.conf_marker):]
				conf += line.strip()
			elif line.startswith( self.pred_marker ) :
				line = line[len(self.pred_marker):]
				pred += line.strip()
			elif line.startswith( self.aa_marker ) :
				line = line[len(self.aa_marker):]
				aa += line.strip()
			else :
				continue

		self.confidence = conf
		self.prediction = pred
		self.sequence = aa
		
		if query :
			if aa == query.sequence :
				pass
			else :
				raise PSIPREDParseError( 'Query sequence do not match with the sequence parsed. \n%s\n%s'%(self.query.sequence, self.sequence ) )
				

		if verbose :
			print(conf)
			print(pred)
			print(aa)
