# -*- coding: utf-8 -*-

def parsePASAStringTieGff(filename):
    '''
    filename is a input of PASA gff output in 2018, which looks like
    scaffold100000_len3256_cov128	assembler	transcript	770	1109	.	-	.	gene_id "PASA_cluster_1"; transcript_id "align_id:836794|asmbl_1";
scaffold100000_len3256_cov128	assembler	exon	770	1109	.	-	.	gene_id "PASA_cluster_1"; transcript_id "align_id:836794|asmbl_1";

scaffold100000_len3256_cov128	assembler	transcript	2017	3200	.	-	.	gene_id "PASA_cluster_2"; transcript_id "align_id:836795|asmbl_2";
scaffold100000_len3256_cov128	assembler	exon	2017	3200	.	-	.	gene_id "PASA_cluster_2"; transcript_id "align_id:836795|asmbl_2";
    
    The same for StringTie gff file, which looks like
    # StringTie version 1.3.4d
scaffold100000_len3256_cov128	StringTie	transcript	2082	3185	1000	.	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; 
scaffold100000_len3256_cov128	StringTie	exon	2082	3185	1000	.	.	gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; exon_number "1"; 
scaffold100004_len29874_cov123	StringTie	transcript	2448	29874	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; 
scaffold100004_len29874_cov123	StringTie	exon	2448	2647	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "1"; 
scaffold100004_len29874_cov123	StringTie	exon	4233	4373	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "2"; 
scaffold100004_len29874_cov123	StringTie	exon	26868	26981	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "3"; 
scaffold100004_len29874_cov123	StringTie	exon	27785	27937	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "4"; 
scaffold100004_len29874_cov123	StringTie	exon	28317	28442	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "5"; 
scaffold100004_len29874_cov123	StringTie	exon	28841	29044	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "6"; 
scaffold100004_len29874_cov123	StringTie	exon	29294	29722	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "7"; 
scaffold100004_len29874_cov123	StringTie	exon	29807	29874	1000	+	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.1"; exon_number "8"; 
    return a dictionary, with scaffold_ids as key, and sub_dictionary as value
        the sub_dictionary has transcript_id as key, and a list of splicing sites as value
    for transcript id, use the part after transcript_id
    '''
    lines = open(filename).readlines()
    
    dc_scf2transcript = {}
    for n, line in enumerate(lines):
        line = line.strip()
        if len(line) < 2:
            continue
        if line[0] == '#':
            continue
        
        elements = line.split('\t')
        
        if len(elements) != 9:
            print('maybe something wrong with line',n,line)
            continue
        
        scf = elements[0]
        if scf not in dc_scf2transcript:
            dc_scf2transcript[scf] = {}
        transcript_id = elements[8].split(';')[1].split()[1].replace('"','')
        if transcript_id not in dc_scf2transcript[scf]:
            dc_scf2transcript[scf][transcript_id] = set()
        site_start = int(elements[3])
        site_end = int(elements[4])
        dc_scf2transcript[scf][transcript_id].add(site_start)
        dc_scf2transcript[scf][transcript_id].add(site_end)
    
    return dc_scf2transcript

def parseSpalnGff(filename):
    '''
    filename is a gff Spaln output, which looks like
    ##gff-version	3
##sequence-region	scaffold25630_len102164_cov122 1 102212
scaffold25630_len102164_cov122	ALN	gene	66473	79023	1589	-	.	ID=gene00001;Name=scaffold25630_len102164_cov122_72
scaffold25630_len102164_cov122	ALN	mRNA	66473	79023	1589	-	.	ID=mRNA00001;Parent=gene00001;Name=scaffold25630_len102164_cov122_72
scaffold25630_len102164_cov122	ALN	cds	77890	79023	1389	-	0	ID=cds00001;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 1 366 +
scaffold25630_len102164_cov122	ALN	cds	77680	77686	21	-	0	ID=cds00002;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 367 368 +
scaffold25630_len102164_cov122	ALN	cds	76288	76289	22	-	2	ID=cds00003;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 369 369 +
scaffold25630_len102164_cov122	ALN	cds	75328	75331	25	-	0	ID=cds00004;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 370 370 +
scaffold25630_len102164_cov122	ALN	cds	74628	74629	16	-	2	ID=cds00005;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 371 371 +
scaffold25630_len102164_cov122	ALN	cds	74570	74573	18	-	0	ID=cds00006;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 372 372 +
scaffold25630_len102164_cov122	ALN	cds	74486	74499	37	-	2	ID=cds00007;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 373 376 +
scaffold25630_len102164_cov122	ALN	cds	67408	67411	27	-	0	ID=cds00008;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 377 377 +
scaffold25630_len102164_cov122	ALN	cds	67297	67299	33	-	2	ID=cds00009;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 378 378 +
scaffold25630_len102164_cov122	ALN	cds	66473	66477	28	-	2	ID=cds00010;Parent=mRNA00001;Name=scaffold25630_len102164_cov122_72;Target=mvi0.2 379 379 +
    return a dictionary, with scaffold_ids as key, and sub_dictionary as value
        the sub_dictionary has transcript_id as key, and a list of splicing sites as value
    for transcript_id, use the part 'Parent=mRNA00001;' and 'Target=mvi0.2' combined, 'mRNA00001_mvi0.2' as transcript_id
    '''
    
    
    
    lines = open(filename).readlines()
    
    dc_scf2transcript = {}
    for n, line in enumerate(lines):
        line = line.strip()
        if len(line) < 2:
            continue
        if line[0] == '#':
            continue
        
        elements = line.split('\t')
        
        if len(elements) != 9:
            print('maybe something wrong with line',n,line)
            continue
        
        fragtype = elements[2]# type of fragments, exon, mRNA, gene or so
        fragtype = fragtype.upper()
        if fragtype != 'CDS' and fragtype != 'EXON':
            continue
        
        scf = elements[0]
        if scf not in dc_scf2transcript:
            dc_scf2transcript[scf] = {}
        descriptions = elements[8].split(';')
        mRNA_id = descriptions[1].split('=')[1]
        target_id = descriptions[3].split('=')[1].split()[0]
        transcript_id = mRNA_id+'_'+target_id
        if transcript_id not in dc_scf2transcript[scf]:
            dc_scf2transcript[scf][transcript_id] = set()
        site_start = int(elements[3])
        site_end = int(elements[4])
        dc_scf2transcript[scf][transcript_id].add(site_start)
        dc_scf2transcript[scf][transcript_id].add(site_end)
    
    return dc_scf2transcript

def get_common_transcripts(file_gff_genome, file_gff_spaln):
    '''
    
    '''
    dc_gff_spaln = parseSpalnGff(file_gff_spaln)
    dc_gff_genome = parsePASAStringTieGff(file_gff_genome)
    
    #filter two based on shared scaffolds
    scf_common = [e for e in dc_gff_genome if e in dc_gff_spaln]
    scf_common = set(scf_common)
    dc_gff_spaln_filter1 = {k:v.copy() for k,v in dc_gff_spaln.items() if k in scf_common}
    dc_gff_genome_filter1 = {k:v.copy() for k,v in dc_gff_genome.items() if k in scf_common}
    
    #delete transcript with only one exon
    dc_gff_spaln_filter2 = {}
    dc_gff_genome_filter2 = {}
    for scf in scf_common:
        dc_gff_spaln_filter2[scf] ={}
        dc_gff_genome_filter2[scf] = {}
        for transcript_id in dc_gff_genome_filter1[scf].keys():
            if len(dc_gff_genome_filter1[scf][transcript_id]) > 2:
                dc_gff_genome_filter2[scf][transcript_id] = dc_gff_genome_filter1[scf][transcript_id]
        for transcript_id in dc_gff_spaln_filter1[scf].keys():
            if len(dc_gff_spaln_filter1[scf][transcript_id]) > 2:
                dc_gff_spaln_filter2[scf][transcript_id] = dc_gff_spaln_filter1[scf][transcript_id]
    
    #for each transcript, change the values from set to sorted tuple. Do not include the first and last one
    for scf in scf_common:
        for transcript_id in dc_gff_genome_filter2[scf]:
            sites = dc_gff_genome_filter2[scf][transcript_id]
            sites = list(sites)
            sites.sort()
            sites = sites[1:-1]
            sites = tuple(sites)
            dc_gff_genome_filter2[scf][transcript_id] = sites
        for transcript_id in dc_gff_spaln_filter2[scf]:
            sites = dc_gff_spaln_filter2[scf][transcript_id]
            sites = list(sites)
            sites.sort()
            sites = sites[1:-1]
            sites = tuple(sites)
            dc_gff_spaln_filter2[scf][transcript_id] = sites
    
    #create dictionary, with scf as key, a set of tuples of splicing sites of different genes as value
    dc_scf2splice_genome = {}
    dc_scf2splice_spaln = {}
    for scf in scf_common:
        dc_scf2splice_genome[scf] = set()
        for sites in dc_gff_genome_filter2[scf].values():
            dc_scf2splice_genome[scf].add(sites)
        dc_scf2splice_spaln[scf] = set()
        for sites in dc_gff_spaln_filter2[scf].values():
            dc_scf2splice_spaln[scf].add(sites)
    
    # for each scaffold, get sites in both genome and spaln
    dc_sites_both = {}
    for scf in scf_common:
        dc_sites_both[scf] = [e for e in dc_scf2splice_genome[scf] if e in dc_scf2splice_spaln[scf]]
    
    # count the number of transcripts that identical in two method
    sites_keep = []
    for sites in dc_sites_both.values():
        sites_keep += sites
    print(len(sites_keep))
    
    # get gene_ids based on dc_sites_both
    ls_transcript_keep = []
    for scf in scf_common:
        _dc_sites_found = {e:0 for e in dc_sites_both[scf]} #only output a single transcript for one combination of sites
        for transcript_id,sites in dc_gff_genome_filter2[scf].items():
            if sites in _dc_sites_found:
                if _dc_sites_found[sites] == 0:
                    ls_transcript_keep.append(transcript_id)
                    _dc_sites_found[sites] = 1
                    
            
    # output the gff3 file with the ls_transcript_keep
    st_transcript_keep = set(ls_transcript_keep)
    lines = open(file_gff_genome).readlines()
    fout = open('good_transcript_genome_homology.gff','w')
    for n, line in enumerate(lines):
        line = line.strip()
        if len(line) < 2:
            continue
        if line[0] == '#':
            continue
        
        elements = line.split('\t')
        
        if len(elements) != 9:
            print('maybe something wrong with line',n,line)
            continue
        
        transcript_id = elements[8].split(';')[1].split()[1].replace('"','')
        if transcript_id in st_transcript_keep:
            fout.write(line+'\n')
    fout.close()

file_gff_genome = 'genome_based.gff'
file_gff_spaln = 'spaln.homolog.gff'
get_common_transcripts(file_gff_genome, file_gff_spaln)