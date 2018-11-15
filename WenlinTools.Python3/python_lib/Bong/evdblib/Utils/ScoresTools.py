import math, sys, os, io

NaN = 1e300000/1e300000
#NaN = float( 'nan' )

def get_variation_using_full_equation( score, mean ) :
        nrecord = len(score)
        rn = 1.0/nrecord ;
        s = 0.0
        for i in score :
                s += rn*(i-mean)**2 ;
        return s

def get_mean_and_stdev( scores ) :
        count = 0 #number of node with non-NA value scores!

        s_score = 0.0
        ss_score = 0.0
        for i in scores :
                if i == 'na' :
                        continue
                count += 1
                s_score += i
                ss_score += i*i

        rcount = 1.0/count
        mean = rcount * s_score
        var = rcount*ss_score - mean*mean
        if var == NaN :
                print("Warnign: NaN detected. Full calculation is used!", file=sys.stderr)
                var = get_variation_using_full_equation( scores, mean, count )

        return mean, math.sqrt( var )

def get_scaled_score( s ) :
        #convert vector s:
        #s[0] = s12
        #s[1] = s11
        #s[2] = s22
        #s[3] = s12r
        s12 = float(s[0])
        s11 = float(s[1])
        s22 = float(s[2])
        s12r = float(s[3])

        sself = (s11+s22)*0.5
        if sself - s12r :
                return  (s12-s12r)/(sself-s12r)
        else :
                return 0.0


def get_z_transformation( score, id1, id2_lists, mu, sigma ) :
	'''
	returns list of Z-scores transformed from score list.
	mu and signma is dictionaries of averages and standard devidations.

	id1 is the query id.
	id2_lists is a list of hit ids.
	'''
        zscore = []
        for id2, s in zip(id2_lists, score):
                mean = (mu[id1] + mu[id2])*0.5
                stdev = math.sqrt( ( sigma[id1]*sigma[id1] + sigma[id2]*sigma[id2] + mu[id1]*mu[id1] + mu[id2]*mu[id2] )*0.5 - mean*mean )
                zscore.append( 1 / (1+math.exp( -(s-mean)/stdev ) ) ) #new addition in transformation.
                                                                #scales Z-score into (0,1) range
        return zscore


def get_z_transformation_pair_id( score, pair_ids, mu, sigma, allow_omission=True ) :
	'''
	returns list of Z-scores transformed from score list.
	mu and signma is dictionaries of averages and standard devidations.

	id1 is the query id.
	id2_lists is a list of hit ids.
	If allow_omission is True,
	id2 is not found among the set, then the z score will be simple z score!
	'''
        zscore = []
        for pairid, s in zip(pair_ids, score):
		id1, id2 = pairid.split()
		if allow_omission and (id2 not in mu or id2 not in sigma) :
			mean = mu[id1] #standard Z score..
			stdev = sigma[id1]
		else :
                	mean = (mu[id1] + mu[id2])*0.5
                	stdev = math.sqrt( ( sigma[id1]*sigma[id1] + sigma[id2]*sigma[id2] + mu[id1]*mu[id1] + mu[id2]*mu[id2] )*0.5 - mean*mean )

                zscore.append( 1 / (1+math.exp( -(s-mean)/stdev ) ) ) #new addition in transformation.
                                                                #scales Z-score into (0,1) range
        return zscore



import pickle
def dump_mean_and_sd_directory( mean_and_sd_dir_root, score_list=[], id_list=[], mean_suffix='.mean', sd_suffix='.std', dump_fn='', verbose=True ) :
	'''
	Dumps all sub-directories mean and sd information
	input two pickled objects

	The mean_and_sd_dir contains list of directories same as 
	score names.
	Each score subdirectories contains list of *.mean and *.sd files.
	Each file cotnains corresponding mean or sd values 
	(single float value in ascii text).

	if score_list and id_list are given,
	they can be used to check the validity of the contents.
	Note that to check the validity 
	both score_name_list and id_list should be given.
	Otherwise, the checking cannot be done!

	If dump_fn is not given,
	this function will dump means and stdevs dictionaries into 
	<mean_and_sd_dir_root>/mean_and_sd.dump file.
	'''

	means = {}
	stdevs = {}

	dir_contents = os.listdir( mean_and_sd_dir_root )

	for pathname in dir_contents :
		#if the directory contains non directory paths,
		#they are ignored!
		if not os.path.isdir( os.path.join(mean_and_sd_dir_root,pathname) ) :
			continue

		if verbose : print("processing...", pathname)
			
		score_means = {}
		score_stdevs = {}
		means[ pathname ] = score_means
		stdevs[ pathname ] = score_stdevs

		score_dir = os.path.join(mean_and_sd_dir_root, pathname )
		for fn in os.listdir(score_dir) :
			#the mean value of scores for a query will be stored in a file
			#<queryid>.mean. And the standard deviation is similarly
			#in <queryid>.sd.
			if fn.endswith( mean_suffix ) :
				query_id = fn[:-len(mean_suffix)]
				fp = open( os.path.join(score_dir,fn) )
				mean_value = float(fp.read())
				score_means[ query_id ] = mean_value
				fp.close()
				#if verbose : print score_dir, fn, 'processed!'
				
			elif fn.endswith( sd_suffix ) :
				query_id = fn[:-len(sd_suffix)]
				fp = open( os.path.join(score_dir,fn) )
				sd_value = float(fp.read())
				score_stdevs[ query_id ] = sd_value
				fp.close()

				#if verbose : print score_dir, fn, 'processed!'

	try :
		for score_name in score_list :
			score_means = means[score_name]
			for query_id in id_list :
				if query_id in score_means :
					pass
				else :
					raise Exception()
	except :
		raise ScoresUtilError( "Mean check failed!" )

	
	try :
		for score_name in score_list :
			score_means = stdevs[score_name]
			for query_id in id_list :
				if query_id in score_means :
					pass
				else :
					raise Exception()
	except :
		raise ScoresUtilError( "Stdev check failed!" )

		
	if not dump_fn :
		dump_fn = os.path.join( mean_and_sd_dir_root, "mean_and_sd.dump" )
	
	if os.path.exists( dump_fn )  :
		raise ScoresUtilError( dump_fn, "is already exists!" )
		
	#dump all scores in one big file
	dump_fp = open( dump_fn, 'w' )
	pickle.dump( means, dump_fp )
	pickle.dump( stdevs, dump_fp )
	dump_fp.close()

	#dump individual scores
	for score_name in means :
		dump_fn = os.path.join( mean_and_sd_dir_root, score_name + '.mean' )
		if os.path.exists( dump_fn ) :
			raise ScoresUtilError( dump_fn + " is already exists!" )
		dump_fp = open( dump_fn, 'w' )
		pickle.dump( means[score_name], dump_fp )
		dump_fp.close()

	for score_name in stdevs :
		dump_fn = os.path.join( mean_and_sd_dir_root, score_name + '.std' )
		if os.path.exists( dump_fn ) :
			raise ScoresUtilError( dump_fn + " is already exists!" )
		dump_fp = open( dump_fn, 'w' )
		pickle.dump( means[score_name], dump_fp )
		dump_fp.close()

	return dump_fn


def load_mean_and_sd_directory( mean_and_sd_dir_root='', dump_fn='' ) :
	'''
	Load Pickled mean and sd values from pickled file.
	If dump_fn is given, the file is used to restore values,
	Otherwise, <mean_and_sd_dir_root>/mean_and_sd.dump file will 
	be used to restore mean and sd values.
	'''
	if not dump_fn :
		dump_fn = os.path.join( mean_and_sd_dir_root, 'mean_and_sd.dump' )

	fp = open(dump_fn )
	means = pickle.load( fp )
	stds = pickle.load( fp )
	return means, stds


class ScoresUtilError( Exception ) :
	pass

######################################
#Constants for ScoreExtractor
######################################
NO_SORT=0
SORT_BY_SCORE=1
SORT_BY_ID=2
######################################

class ScoreExtractor :
	'''
	This class extract scores from PairwiseAlignmentRecords object.

	Currently the sorting option is not implemented!
	'''

	def __init__( self, alignments=None, sortby=NO_SORT ) :
		'''
		'''
		#main object to extract scores from
		self.alignments = alignments
	
		self.sortby = sortby
		if self.sortby :
			raise NotYetImplemented( "Sortby option is not working yet!" )

		#maybe the option for ordering by the the given alignments
		#is needed in future!

	def extract_scores( self, method_name, multiple=True ) :
		'''
		Returns a list of lists.
		Returned list is formatted like the followings;
		[ [raw_score, normalized_score, id1, id2], ... ]

		If mulitple option is true,
		the method name is matched promiscuously,
		and the maximum normalized score record results will be
		returned.
		'''
		scores = []
		for id1 in self.alignments :
			queryrec = self.alignments[id1]
			for id2 in queryrec :
				hitrec = queryrec[id2]
				method_records = hitrec.get_method_records(method_name, exact=not multiple)

				temp_scores = []
				for method_rec in method_records :
					if method_rec.is_null() :
						continue
					temp_scores.append( ( float(method_rec.get_normalized_score()), float(method_rec.get_raw_score()) ) )

				
				if temp_scores :
					(norm_score, raw_score) = max(temp_scores)
					scores.append( [raw_score, norm_score, id1, id2, method_rec.additional_scores] )
		return scores


	def format( self, scores, fp=None, format="%(id1)s\t%(id2)\t%(raw_score)s\t%(norm_score)s" ) :
		'''
		Format the scores, list of list as returned from extract_scores function.

		The format string follows the convention of python
		string formatting with dictionary.
		And the keywords are "id1", "id2", "raw_score" and "norm_score".

		Note that the "new line" character should be placed at that end of the format string
		if you want a line for scores for the specified pair.
		'''
		if not fp :
			fp = io.StringIO()

		for score in scores :
			raw_score, norm_score, id1, id2, additional_scores = score

			if additional_scores :
				#to prevent errors happening by
				#None-type additional_scores
				format_dictionary = {}
				format_dictionary.update(additional_scores)
				format_dictionary.update(locals())
			else :
				format_dictionary = locals()

			print(format%format_dictionary, file=fp)

		return fp
