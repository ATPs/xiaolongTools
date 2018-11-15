'''
This module contains tools related to
mapping PDB sequence to various sequences
dervied from PDB structure.

This module is mainly focused on the clean
process of PDB sequence using PDB sub-package.
'''
import sys, os, io
verbose = 0

from evdblib.Utils.Parsers.Range import PDBRange, SequenceRange
#Histag checking class
from evdblib.Utils import TagChecker
histagchecker = TagChecker()

from evdblib.Utils.Algorithms.SequenceAlignment import align_sequences_with_affine_gap
from evdblib.Runners.StructureTools import MaxSproutRunner

class PDBMap :
	def __init__( self, pdb=None, pdbrange=None, reconstitute_backbone=False, reference_type='bc', atom_bio_diff_cutoff=15, remove_histag=True ) :
		'''
		Generate and keep mappings between various sequence
		derived from PDB files.

		If recontitute_backbone flag is True,
		reconstitution of backbone atoms will be done,
		using maxsprout program.

		Note that the backbone reconstitution part is not yet 
		implemented!

		reference_type can be [b][q,a,c,f].
		"b" represents filtered by biological sequence defined by
		DBREF record.
		"q" represents seqres record.
		"a", "c", "f" represent residues having any ATOM records, 
		CA ATOM records, or full backbone ATOM records.
		The default value is set to "bc" when the difference between "a" and "ba"
		is different less than 15 residues (atom_bio_diff_cutoff).
		atom_bio_diff_cutoff is used only 
		if referency_type includes "b" or biological sequences.
		If the difference is bigger than the cutoff,
		"b" option will be ignored.

		If remove_histag is True,
		the seqres sequence is checked if it contains a histidine tag.
		And the tag region will be removed from the reference sequence.
		'''
		self.seqres = None
		self.histag = None
		self.bioseq = None
		self.atmseq = None
		self.caonly = None
		self.fullbb = None
		self.refseq = None

		self.missing_character = '='
		self.sequencerange = None

		#using maxsprout program to reconstitiute
		#backbone coordinates
		self.mxspbb = None 
		self.mxspbb_records = None 
		#this variable is for marking 
		#maxsprouting or not.
		self.reconstitute_backbone = reconstitute_backbone

		self.pdb = pdb
		if pdb != None and pdb.fn :
			self.pdbfile = os.path.abspath(pdb.fn) #this is something can be used when read from the path.
		else :
			self.pdbfile = None

		self.pdbrange = pdbrange

		self.reference_type = reference_type
		#cutoff for not using bioseq range.
		self.atom_bio_diff_cutoff = atom_bio_diff_cutoff 
		self.remove_histag = remove_histag

		#if the PDB and PDBRange object are given,
		#run mapping file.
		if self.pdb and self.pdbrange :
			self.generate_mapping( self.pdb, self.pdbrange, self.reconstitute_backbone )

			self.set_histag()
			self.set_refseq()

	def check_atom_bioatom_difference( self ) :
		'''
		returns true if the difference
		is smaller than self.atom_bio_diff_cutoff.
		Otherwise returns False.
		'''

		#getting the atom&bio sequence.
		atombioseq = self.get_intersection( self.atmseq, self.bioseq ).replace( self.missing_character, '' )
		atmseq = self.atmseq.replace( self.missing_character, '' ) 

		if len(atmseq)-len(atombioseq) > self.atom_bio_diff_cutoff :
			return False
		else :
			return True


	def set_refseq( self ) :
		'''
		Setting the variable for reference sequence (refseq).
		Note that the reference sequece setting can be controlled
		by three parameters supposed set by the __init__ function.
		1. reference_type
		2. atom_bio_diff_cutoff
		3. remove_histag

		Note that this function raises ValueErorr when the
		refseq value becomes None, due to bugs or nogo reference type
		settings.
		'''
		reference_type = self.reference_type
		remove_histag = self.remove_histag

		if not self.check_atom_bioatom_difference() :
			#turn off "b" or biological sequence
			reference_type = reference_type.replace( "b", '' ) 
		
		refseq = None

		if 'a' in reference_type :
			refseq = self.atmseq
		elif 'c' in reference_type :
			refseq = self.caonly
		elif 'f' in reference_type :
			refseq = self.fullbb
		elif 'm' in reference_type :
			refseq = self.mxspbb
		elif 'q' in reference_type :
			refseq = self.refseq

		if refseq == None :
			raise ValueError( reference_type, "Reference type setting should be done correctly!" )

		if 'b' in reference_type :
			refseq = self.get_intersection( refseq, self.bioseq )

		if remove_histag :
			inversed_histag = self.get_inverse( self.histag )
			refseq = self.get_intersection( refseq, inversed_histag )
		
		if refseq == None :
			raise ValueError( refseq, "Reference sequence value should not be None!" )

		self.refseq = refseq



	def get_intersection( self, seq1, seq2  ) :
		'''
		returns the residues common 
		(alphabet characters at the same postions) 
		in both seq1 and seq2.
		For positions that are not aligned,
		the missing character (default: =)
		will be placed.

		Note that seq1 and seq2 should be pre-aligned.
		'''
		seq_fp = io.StringIO()
		for a, b in zip( seq1, seq2 ) :
			if a.isalpha() and b.isalpha() :
				seq_fp.write( a )
			else :
				seq_fp.write( self.missing_character )

		return seq_fp.getvalue()

	def get_union( self, seq1, seq2 ) :
		'''
		returns the union of the two sequences.

		Similar to get_intersection.
		'''
		seq_fp = io.StringIO()
		for a, b in zip( seq1, seq2 ) :
			if a.isalpha() :
				seq_fp.write( a )
			elif b.isalpha() :
				seq_fp.write( b )
			else :
				seq_fp.write( self.missing_character )

		return seq_fp.getvalue()

	def get_inverse( self, seq1, seq2=None ) :
		'''
		returns the sequence of residues not present in seq1
		from seq2.
		If the residues are present in seq1, then
		change them into self.missing_character
		By default, seq2 will be the seqres sequence.

		Note that the seq1 and seq2 are pre-aligned.
		'''

		if seq2== None :
			seq2 = self.seqres

		fp = io.StringIO()
		for a, b in zip( seq1, seq2 ) :
			if a.isalpha() :
				fp.write( self.missing_character )
			else :
				fp.write( b )

		return fp.getvalue()

		


	def get_intersection_indices( self, seq1, seq2, reference=1 ) :
		'''
		returns indices for the common positions in seq1 and seq2.
		Reference parameter 1 or 2 designates seq1 or seq2, respectively.
		Note that seq1 and seq2 should be pre-aligned.
		'''
		indices = []
		i = -1
		for a, b in zip( seq1, seq2 ) :
			if reference == 1 :
				if a.isalpha() :
					i += 1
			elif reference == 2 :
				if b.isalpha() :
					i += 1
				
			if a.isalpha() and b.isalpha() :
				indices.append( i )

		return indices 

	def get_equivalent_positions( self, seq1, seq2 ) :
		'''
		returns indices for the common positions in seq1 and seq2
		relative to themelves.
		'''
		eq1, eq2 = [], []
		i, j = -1, -1
		for a, b in zip( seq1, seq2 ) :
			if a.isalpha() :
				i += 1
			if b.isalpha() :
				j += 1

			if a.isalpha() and b.isalpha() :
				eq1.append(i)
				eq2.append(j)

		return [eq1, eq2]


	def get_difference_indices( self, seq1, seq2 ) :
		'''
		returns the indices of the resdiues available 
		in seq1 but not in seq2.
		Note that seq1 and seq2 should be pre-aligned.
		'''
		indices = []
		i = -1
		for a, b in zip( seq1, seq2 ) :
			if a.isalpha() :
				i += 1
				
			if a.isalpha() and not b.isalpha() :
				indices.append( i )

		return indices 


	def get_difference( self, seq1, seq2 ) :
		'''
		returns the resdiues available in seq1 but not in seq2.
		Note that seq1 and seq2 should be pre-aligned.
		'''
		indices = self.get_difference_indices( seq1, seq2 )
		seq_fp = io.StringIO()
		for i in indices :
			seq_fp.write( seq1[i] )

		return seq_fp.getvalue()

	
	def set_custom_sequence( self, seq ) :
		'''
		Set the new custom sequence supposed to be partial or slightly different
		from the SEQRES sequence.
		The sequence will be aligned to the seqres sequence 
		to prepare to filter other alingned sequences same as the custom sequence.
		Then the common positions between the custom sequence and reference sequence
		will the extracted and saved.
		
		'''
		alignment = align_sequences_with_affine_gap( self.seqres, seq )
		score = alignment[2]
		if score < len(alignment[0]) :
			raise ValueError( score, seq, "Sequences are too different!" )

		eq = self.get_equivalen_positions( alignment[0], alignment[1] )
		filtered_eq = [[], []] #mapping between reference and custom
		filtered_ca_eq = [[],[]] #mapping between reference/ca and custom
		for i, j in zip( eq[0], eq[1] ) :
			refseq_residue = self.refseq[i]
			if refseq_residue.isalpha() :
				filtered_eq[0].append(i)
				filtered_eq[1].append(j)

			if self.caonly[i].isalpha() :
				filtered_ca_eq[0].append(i)
				filtered_ca_eq[1].append(j)

		self.cstseq = seq
		self.refseq_index2customseq_index = filtered_eq
		self.refseqca_index2customseq_index = filterd_ca_eq

	def filter_aligned_sequence( self, aln, offset=0, filter_type="ref" ) :
		'''
		returns filtered aligned_sequence.

		This function requires the setting of aligned sequence positions
		using set_custom_sequence function.

		filter_type can be one of the two:
		"ref" for filter against the reference sequence
		or "caref" for filter against the CA only and the reference sequence.

		synopsis:
		pdbmap = PDBMapping( ... )
		pdbmap.set_custom_sequence( custom_seq )
		pdbmap.filter_aligned_sequence( aln )

		#for filter the custom sequence
		filtered_custom_seq = pdbmap.filter_aligned_sequence( custom_seq )
		'''

		feq = None
		if filter_type == 'ref' :
			feq = self.refseq_index2customseq_index
		elif filter_type == 'caref' :
			feq = self.refseqca_index2customseq_index

		#make filtered index set
		filtered_indices = set(feq[1])
		new_aln_fp = io.StringIO()
		i = -1
		for a in aln :
			if a.isalpha() :
				i += 1
				if i + offset in filtered_indices :
					new_aln_fp.write( a )
				else :
					new_aln_fp.write( '-' )
			else :
				new_aln_fp.write( a )

		return new_aln_fp.getvalue()


	def filter_caonly_records( self ) :
		'''
		returns filtered CA only PDB ATOM records.

		This function requires "set_custom_sequence" function already
		run.
		'''

		caonly_records = self.caonly_records.strip().split('\n')
		if len(caonly_records) != len(self.caonly) :
			raise IndexError( "CAonly records should be same CAonly sequence." )

		fp = io.StringIO()
		for i in self.refseqca_index2customseq_index[0] :
			print(caonly_records[i], file=fp)

		return fp.getvalue()
		

	def set_histag( self ) :
		'''
		check histidine tag on sequence based on seqres sequence.
		And set the histag variable.
		'''
		tagmatches = histagchecker.check( self.seqres )
		histag = self.missing_character * len(self.seqres) 

		for tagmatch in tagmatches :
			start, end = tagmatch.span()
			newhistag = self.missing_character * start + tagmatch.group() + self.missing_character * (len(self.seqres)-end)
			histag = self.get_union( histag, newhistag )
			
		self.histag = histag


	def generate_mapping( self, pdb, pdbrange, reconstitute_backbone=False ) :
		'''
		generates five types of sequences.
		'''
		seqres_residues = pdb.extract_seqres_residues( pdbrange )

		bioseq_residues = []
		atmseq_residues = []
		#caonly_residues = []
		#fullbb_residues = []

		seqres_fp = io.StringIO()
		bioseq_fp = io.StringIO()
		atmseq_fp = io.StringIO()
		caonly_fp = io.StringIO()
		fullbb_fp = io.StringIO()

		mc = self.missing_character
	
		#first time run through the loops. :)
		for residue in seqres_residues :
			standard_aa_letter = residue.get_standard_amino_acid_name_letter()

			seqres_fp.write( standard_aa_letter )

			if residue.is_biological_residue() :
				bioseq_residues.append( residue )
				bioseq_fp.write( standard_aa_letter )
			else :
                                bioseq_fp.write( mc )

			if not residue.is_missing_residue() :
				atmseq_fp.write( standard_aa_letter )
				atmseq_residues.append( residue )
			else :
				atmseq_fp.write( mc )

			if residue.has_ca() :
				#caonly_residues.append( residue )
				caonly_fp.write( standard_aa_letter )
			else :
				caonly_fp.write( mc )

			if residue.has_full_backbone_atoms() :
				#fullbb_residues.append( residue )
				fullbb_fp.write( standard_aa_letter )
			else :
				fullbb_fp.write( mc )


		if verbose :

			print('DOMRNG', self.pdbrange, pdb.pdbrange2sequencerange( self.pdbrange, biological=False, atomrecord=False ))
			print('SEQRES', seqres_fp.getvalue())
			print('BIOSEQ', bioseq_fp.getvalue())
			print('ATMSEQ', atmseq_fp.getvalue())
			print('CAONLY', caonly_fp.getvalue())
			print('FULLBB', fullbb_fp.getvalue())


		self.sequencerange = pdb.pdbrange2sequencerange( self.pdbrange, biological=False, atomrecord=False )
		self.seqres = seqres_fp.getvalue()
		self.bioseq = bioseq_fp.getvalue()
		self.atmseq = atmseq_fp.getvalue()
		self.caonly = caonly_fp.getvalue()
		self.fullbb = fullbb_fp.getvalue()


		#####################
		#clean up the strings.
		caonly_fp.reset()
		caonly_fp.truncate()

		fullbb_fp.reset()
		fullbb_fp.truncate()
		#end of cleaning.
		#####################

		#second time run through the loop.
		#it is to build ATOM records. 
		ca_serial_number = 1
		bb_serial_number = 1
		for i, residue in enumerate( atmseq_residues ) :
			resnum = i+1

			if residue.has_ca() :
				caonly_fp.write( residue.get_ca_record( resnum=resnum, alanine_fall_back=True, number=ca_serial_number, chainid='A', icode=' ', altloc=' ' ))
				ca_serial_number += 1

			if residue.has_full_backbone_atoms() :
				fullbb_fp.write( residue.get_atom_records( resnum=resnum, alanine_fall_back=True, number=bb_serial_number, chainid='A', icode=' ', altloc=' ' ))
				bb_serial_number += len(residue)

		if verbose :
			print(caonly_fp.getvalue())
			print(fullbb_fp.getvalue())

		self.caonly_records = caonly_fp.getvalue()
		self.fullbb_records = fullbb_fp.getvalue()


	def save_mapping_file( self, basename, save_dir=None, mapping_suffix='.map', caonly_suffix='.ca', fullbb_suffix='.bb', mxspbb_suffix='.msbb', update=False ) :
		'''
		save mapping files and CAonly file.
		'''

		if save_dir :
			if os.path.exists( save_dir ) :
				pass
			else :
				os.makedirs( save_dir )

			basename = os.path.join( save_dir, basename )

		mapping_fn = basename+mapping_suffix
		caonly_fn = basename+caonly_suffix
		fullbb_fn = basename+fullbb_suffix
		mxspbb_fn = basename+mxspbb_suffix

		if not update and os.path.exists( mapping_fn ) :
			raise IOError( "Mapping file %s already exists." % mapping_fn )

		mapping_fp = open( mapping_fn, 'w' )
		if self.pdbfile :
			print('PDB', self.pdbfile, file=mapping_fp)
		print('DOMRNG', self.pdbrange, self.sequencerange, file=mapping_fp)
		print('SEQRES', self.seqres, file=mapping_fp)
		print('HISTAG', self.histag, file=mapping_fp)
		print('BIOSEQ', self.bioseq, file=mapping_fp)
		print('ATMSEQ', self.atmseq, file=mapping_fp)
		print('CAONLY', self.caonly, file=mapping_fp)
		print('FULLBB', self.fullbb, file=mapping_fp)
		if self.mxspbb :
			print('MXSPBB', self.mxspbb, file=mapping_fp)
		print('REFSEQ', self.refseq, file=mapping_fp)
		mapping_fp.close()

		if not update and os.path.exists( caonly_fn ) :
			raise IOError( "CAonly file %s already exists." % caonly_fn )
		caonly_fp = open( caonly_fn, 'w' ) 
		caonly_fp.write( self.caonly_records )
		caonly_fp.close()

		if not update and os.path.exists( fullbb_fn ) :
			raise IOError( "Full Backbone file %s already exists." % fullbb_fn )
		fullbb_fp = open( fullbb_fn,'w' ) 
		fullbb_fp.write( self.fullbb_records )
		fullbb_fp.close()

		#stop here if max sprouting savin
		#is not happening.
		if not mxspbb_fn :
			return

		if self.mxspbb_records :
			if not update and os.path.exists( mxspbb_fn ) :
				raise IOError( "Full Backbone file %s already exists." % mxspbb_fn )
			mxspbb_fp = open( mxspbb_fn,'w' ) 
			mxspbb_fp.write( self.mxspbb_records )
			mxspbb_fp.close()


	def parse_mapping_file( self, mapping_fn=None, caonly_fn=None, fullbb_fn=None, mxspbb_fn=None ) :
		'''
		read mapping file and other coordinate files.
		'''
		fp = open( mapping_fn )
		line = fp.readline()
		while line :
			header = line.split()[0]
			if header == 'PDB' :
				pdbfile = ' '.join( line.split()[1:] )
				self.pdbfile = pdbfile
				
			elif header == 'DOMRNG' :
				pdbrange, sequencerange = line.split()[1:]
				self.pdbrange = PDBRange()
				self.pdbrange.parse(pdbrange)

				self.sequencerange = SequenceRange()
				self.sequencerange.parse( sequencerange )
				
			elif header == 'SEQRES' :
				self.seqres = line.split()[1]
			elif header == 'HISTAG' :
				self.histag = line.split()[1]
			elif header == 'ATMSEQ' :
				self.atmseq = line.split()[1]
			elif header == 'BIOSEQ' :
				self.bioseq = line.split()[1]
			elif header == 'CAONLY' :
				self.caonly = line.split()[1]
			elif header == 'FULLBB' :
				self.fullbb = line.split()[1]
			elif header == 'MXSPBB' :
				self.mxspbb = line.split()[1]
			elif header == 'REFSEQ' :
				self.refseq = line.split()[1]

			line = fp.readline()
		fp.close()

		#read CAonly ATOM records
		if caonly_fn :
			fp = open( caonly_fn )
			self.caonly_records = fp.read()
			fp.close()

		#read Full Backbone ATOM records
		if fullbb_fn :
			fp = open( fullbb_fn )
			self.fullbb_records = fp.read()
			fp.close()

		if mxspbb_fn :
			fp = open( mxspbb_fn )
			self.mxspbb_records = fp.read()
			fp.close()


	def reconstitute_missing_backbone( self ) :
		'''
		Run MaxSprout program developed by Liisa Holm.
		and check the generated backbone atoms.
		'''
		maxsproutrunner = MaxSproutRunner( inputpdb=self.pdb )
		maxsproutrunner.run()
		maxsprouted_pdb = PDB( maxsrpoutrunner.output_fn )
