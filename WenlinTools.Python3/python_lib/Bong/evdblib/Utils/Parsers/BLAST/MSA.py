from evdblib.Utils.Parsers.BLAST import *
from evdblib.Utils.Parsers.FASTA import FASTA, parse_a3m
from .Hit import Hit
from evdblib.Utils.SimpleThreadingTools import ThreadPool

class MSAInitializationError( Exception ) :
	pass

class MSA :
	'''
	A class encapsulating a BLAST result or PSIBLAST iteration.
	'''
	def __init__( self, query=None, blast_result_fp=None, combine_hsps=False, buildalia3m=False ) :
		'''
		buildalia3m: a flag for parsing A3M file from buildali program.
		'''

		#if both objects are not available,
		#MSA class cannot contain any meaningful data.
		if buildalia3m :
			if not blast_result_fp :
				raise MSAInitializationError( "BuildAli A3M file pointer is needed for parsing")

		elif not query or not blast_result_fp  :
			raise MSAInitializationError( "query and blast_result_fp are needed for parsing.")

		#blast hits.
		#main container in this class
		self.hits = []

		#query information
		self.query = query

		#convinient access of the hit.
		#?would this be necessary?
		self.gi2hit = {}

		#?saving parsing result ?
		self.parsing_status = None
		
		#saving the given file pointer
		self.fp = blast_result_fp
		self.initial_fp_position = self.fp.tell()

		self.combine_hsps = combine_hsps  #default value
		#?would it be useful?
	
		#calling main parser of the class
		#self.parse()

		########################################
		#?this should be defined in the BLAST class?
		self.effective_length_of_the_database = 0 
		########################################

		#deprecated!
		#self.annotation = ""
		#self.query_seq = ""
		#self.blast_result_fn = blast_result_fn
		#self.annotation, self.query_seq = self.get_query( seq_fn )
		#self.dummy_db = ''
	
		#buildalia3m 
		#specific variable
		self.non_protein_fastas = []

	def set_combine_hsps( self ) :
		'''
		sets the flag of combine hsps true.
		'''
		self.combine_hsps = True

	def unset_combine_hsps( self ) :
		self.combine_hsps = False


	def parse( self ) :
		while True:
			try :
				hit = Hit() 
				hit.parse( self.fp ) 

			except NoMoreIteration :
				if hit.has_data() :
					self.hits.append( hit )
				raise

			except NoMoreHit as e :
				if hit.has_data() :
					if verbose :
						print("Last hit is saved!")
					self.hits.append( hit )
				else :
					if verbose :
						print("Last hit is not saved!")

				if verbose :
					print(e)
				break

			except FormatError as e :
				print('format error occurred!', file=sys.stderr)
				print(e, file=sys.stderr)
				#without saving the hit.
			else :
				self.hits.append( hit )


	def parse_a3m_fasta( self ) :
		'''
		Parse BuildAli A3M file.
		'''
		fastas = parse_a3m( descriptor=self.fp )
		protein_fastas = []

		#detect secondary structure records
		#and other non-protein fastas
		for f in fastas :
			if f.is_protein_fasta() :
				protein_fastas.append(f)
			else :
				self.non_protein_fastas.append( f )

		gi2hit = {}
		#query fasta setting
		self.query = protein_fastas.pop(0)
		for f in protein_fastas :
			if f.gi in gi2hit :
				hit = gi2hit[f.gi]
			else :
				hit = Hit()

			hit.parse_a3m_fasta( f, self.query )
			if hit.has_data()  :
				pass
			else :
				if verbsoe :
					print('Last hit in a3m is not saved!')
				continue


			if not f.gi in gi2hit :
				self.hits.append( hit )


	def __len__( self ) :
		'''
		len( msa ) -> return number of hits
		'''
		return self.hits.__len__()

	def __getitem__( self, i ) :
		'''
		msa[i] -> return i'th hit object
		'''
		return self.hits[i]

	def __iter__( self ) :
		'''
		iteration over the hit list.
		'''
		return self.hits.__iter__()

	def __str__(self) :
		fp = cStringIO.StringIO()
		for i, hit in enumerate( self ) :
			fp.write('hit #%d\n'%(i+1))
			fp.write(str(hit))
			fp.write('\n')
		return fp.getvalue()


	def build_a3m_fasta( self, write_fp, remove_gap=False, combine_hsps=None ) :
		'''
		Builds A3M MSA file.
		This can be used for buildali.pl program input file.
		'''
		query_fasta = self.query
		query_seq = query_fasta.sequence
		query_fasta.save( fp=write_fp )
		for hit in self.hits :
			hit.write_a3m_fasta( write_fp, query_fasta.sequence, remove_gap=remove_gap, combine_hsps=combine_hsps )


	def build_psiblast_alignment_input( self, write_fp, remove_gap=True, combine_hsps=None ) :
		"""
		Builds PSIBLAST Input MSA file.
		"""
		query_fasta = self.query
		query_seq = query_fasta.sequence 
		annotation = '_'.join( query_fasta.header[1:].strip().split() )

		if combine_hsps != None :
			combine_hsps = self.combine_hsps

		if remove_gap == True :

			write_fp.write( '%-20s %s\n'%( annotation[:20], query_seq ) )
			for hit in self.hits :
				aligned_seq = hit.get_aligned_hit_sequence( len(query_seq), combine_hsps=self.combine_hsps )
				write_fp.write( aligned_seq )

		else :
			max_gap_info = {}
			for hit in self.hits :
				gap_info = hit.get_gapped_positions()
				for k,v in gap_info.items() :
					if k in max_gap_info :
						max_gap_info[k] = max( v, max_gap_info[k] )
					else :
						max_gap_info[k] = v

			write_fp.write( '%-20s %s\n'%( annotation[:20], query_seq ) )
			for hit in self.hits :
				#deliver information to hit
				#and then the hit will deilver to ....
				aligned_seq = hit.get_aligned_hit_sequence( len(query_seq), combine_hsps=self.combine_hsps, max_gap_info=max_gap_info )
				write_fp.write( aligned_seq )
				
		write_fp.flush()


	def purge_overlapping_hsps( self, positions, pssm ) :
                '''
                This function makes modification ont he psiblast result by removing
                possibly contaminated regions in the hsps.
                Check which portion of the HSP are good if the given hsps are overlapping.
                Experiments #4
                '''

                #print '###Start to purge overlapping hits...'
                for hit in self.hits[:] :
                        hit.purge_overlapping_hsps( positions, pssm )

                        ###################################
                        #remove the hit if selection failed
                        ###################################
                        if not hit :
                                self.hits.remove( hit )

	def purge_overlapping_hsps_multithreading( self, positions, pssm, nthreads ) :
		'''
		This function works as purge_overlapping_hsps 
		but using multithreading.
	
		The number of concurrent threads can be defined by nthreads.
		'''

		pool = ThreadPool( nthreads )
		for hit in self :
			pool.add_task( hit.purge_overlapping_hsps, positions, pssm )

		pool.wait_completion()

		for hit in self.hits[:] :
			self.hits.remove( hit )


	def remove_overlapped_hit_regions( self, blastmsa ) :
		"""
		removes hit region where overlap with the hits in blastmsa.
                """

                msa_gi_set1 = self.get_hit_gi_set()
                msa_gi_set2 = blastmsa.get_hit_gi_set()

                if verbose :
                        print("msa_gi_set1", msa_gi_set1)
                        print("")
                        print("msa_gi_set2", msa_gi_set2)
                        print("")

                common_gis = msa_gi_set2.intersection( msa_gi_set1 )
                if verbose :
                        print('common_gis:', common_gis)

                for gi in common_gis :
                        if verbose:
                                print(gi, "will be purged")
                        hit1 = self.get_hit( gi )
                        hit2 = blastmsa.get_hit( gi )

                        hit1.remove_overlapped_hit_region( hit2 )

			if not len(hit1) :
				self.hits.remove(hit1)


	def remove_overlapped_hit_regions_multithreading( self, blastmsa, nthreads ) :
		"""
		removes hit region where overlap with the hits in blastmss
		using multiple threads.
                """

                msa_gi_set1 = self.get_hit_gi_set()
                msa_gi_set2 = blastmsa.get_hit_gi_set()

                if verbose :
                        print("msa_gi_set1", msa_gi_set1)
                        print("")
                        print("msa_gi_set2", msa_gi_set2)
                        print("")

                common_gis = msa_gi_set2.intersection( msa_gi_set1 )


		pool = ThreadPool( nthreads )
                for gi in common_gis :
                        hit1 = self.get_hit( gi )
                        hit2 = blastmsa.get_hit( gi )

                        hit1.remove_overlapped_hit_region( hit2 )
			pool.add_task( hit1.remove_overlapped_hit_region, hit2 )

		pool.wait_completion()
		
		for gi in common_gis :
			hit1 = self.get_hit( gi )
			if not len( hit1 ) :
				self.hits.remove(hit1)


        def remove_overlapping_hsps( self ) :
                """
                This function makes modification on the psiblast result by removing
                possibly contaminated Hits.
                """
                for hit in self.hits[:] :
                        hit.remove_overlapping_hsps()
                        if not len(hit) :
                                self.hits.remove( hit )


	def get_hit_gi_set( self ) :
		'''
		returns a set of GIs from hits.
		'''
 	        id_set = set()
                for hit in self.hits :
                        id_set.add( hit.gi )

                return id_set

		
	def build_gi2hit( self ) :
                '''
                This function make unique gi index called gi2hit.
                As a byproduct, this function also make compaction in the hit list
                to make all hits are uniquely found by a give gi number.
                '''
                for hit in self.hits[:] :
                        if hit.gi in self.gi2hit :
                                self.gi2hit[ hit.gi ].add_hit_as_hsp( hit )
                                self.hits.remove( hit )
                        else :  
                                self.gi2hit[ hit.gi ] = hit

	def write_hit_fasta( self, fp ) :
		for hit in self :
			hit.write_hit_fasta( fp )

        def get_hit( self, gi ) :
		'''
		returns a hit deteremined by the GI value.
		'''
                if not self.gi2hit :
                        self.build_gi2hit()

                if gi in self.gi2hit :
                        return self.gi2hit[gi]
                else :  
                        return None

		
	'''
	#deprecated from HangOut original version
	def parse( self, use_fp = 0, psiblast_flag = 1 ) :
		blast_hit = BlastHit()
		fp = None
		if use_fp :
			fp = self.blast_result_fp
		else :
			fp = open( self.blast_result_fn )
		if psiblast_flag == 1 :
			best_iteration_start, selected_iteration = select_last_iteration( fp )
			print >>sys.stderr, "Iteration", selected_iteration, " is selected!"
			print >>sys.stderr, "Starts at", best_iteration_start, "!"
		
			fp.seek( best_iteration_start ) #pointing to the start position of final iteration
					
		if verbose :
			print '##### after psi-blast result parsing routine ######'
			print 'current file pointer, fp.tell():', fp.tell()
		parse_flag = blast_hit.parse( fp )
		if verbose :
			print 'first blast hit parsing is done!'
			print 'parse_flag:',parse_flag
			print 'blast_hit.gi:', blast_hit.gi
		while parse_flag == 1 :
			self.hits.append( blast_hit )
			blast_hit = BlastHit()
			parse_flag = blast_hit.parse( fp )
			if parse_flag == -2 :
				if verbose :
					print 'End of Blast Result is processed!'
				self.hits.append( blast_hit )
				break
			if verbose :
				print 'parse_flag @ BlastMSA parse :', parse_flag
			
			
		self.parsing_status = parse_flag
		if parse_flag == -1 :
			print >>sys.stderr, "WARNING: Erroneous Parsing of a Blast Hit.", parse_flag
			print >>sys.stderr, "This possibly indicate premature finish of a Blast/PSI Blast result."
			print >>sys.stderr, fp.readline()
			return False
		elif parse_flag == 0 or parse_flag == -2 :
			pass #everything is done
		elif parse_flag == 1 :
			#something is not complete
			print >>sys.stderr, "WARNING: Hit records are all correctly parsed, but the result may not contain all data."
		else :
			print >>sys.stderr, "Error! This message should not be seen!"
		#Need to parse effective database size
		#code will be here!!
			
		fp.close()
		return True
	'''
