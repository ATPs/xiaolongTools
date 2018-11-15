"""
This module contains a class for A3M file type.
A3M is a type of multiple sequence alignment format,
with slightly condensation by putting lower case letters,
used in Soding's HHsearch package.
"""
from evdblib.Utils.Parsers.FASTA import parse_multiple
from io import StringIO

class BuildAliA3M :
	'''
	class for A3M parsing.
	'''
	def __init__( self, fn=None ) :
		'''
		fn: a3m filename
		buildali: when the a3m file is from buildali program in HHsearch.
		'''
		self.query = None
		self.hits = []
		self.additional_records = []
		self.buildaliflag = True
		self.fn = fn

		if self.fn and self.buildaliflag :
			self.parse_buildali_a3m( self.fn )

	def parse_buildali_a3m( self, fn ) :
		self.fastas = parse_multiple( file=fn, checksequence=False  )

	def __str__( self ) :
		s = StringIO()
		for f in self.fastas :
			s.write( f.__str__() )
		return s.getvalue()

