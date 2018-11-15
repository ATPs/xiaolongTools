#this is a lib to deal with the synchronization problem for the barcode database
#here we use a sql machenizion to prevent overwrite
import os
import sqlite3 as sql
from fullname_lib import get_IDlist, get_SRNPnumber

def init_table():
    #this is the command has run to create the table
    #just put it here for archive
    conn = sql.connect(fall)
    c = conn.cursor()
    c.execute("""CREATE TABLE barcode_fasta (sampleID text, CSid text, isPCR INTEGER, genus text, species text, full_name text PRIMARY KEY, sequence text)""")
    conn.commit()
    conn.close()
    os.system('chmod a+w %s' % fall)

#the file saving all barcodes
fall = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/data/barcodes/all_barcodes.db'
if not os.path.exists(fall):
    print('init table...')
    init_table()


sequenced_set = set(get_IDlist())

CSid_set = set(get_SRNPnumber())


def parse_ID(ID):
    ID = str(ID)
    removed_labels = ['11-BOA-', 'NVG-']
    for label in removed_labels:
        ID = ID.replace(label, '')
    ID = ID.strip()
    return ID
sequenced_set = set([parse_ID(ID) for ID in sequenced_set])


def parse_name(name):
    removed_syntax = ['"', '\'', '?']
    replace_syntax = {
            ' ': '_',
            '\t': '_',
            }

    for char in removed_syntax:
        name = name.replace(char, '')

    for char1 in replace_syntax:
        char2 = replace_syntax[char1]
        name = name.replace(char1, char2)

    return name


def find_genus_and_species(Oname):
    genus, sp = 'NA', 'NA'

    name = Oname.split('|')[0]
    items = name.strip().replace(' ', '_').replace('\t', '_').split('_')
    genus = items[0]

    try:
        if parse_ID(genus) in sequenced_set:
            genus, sp = items[1:3]
        else:
            sp = items[1]
    except:
        if 'COI-5P' in Oname:
            try:
                items = Oname.split('|')[-2].split('_')[:2]
                if len(items) == 1:
                    genus = items[0]
                else:
                    genus, sp = items[:2]
            except:
                pass

            if 'JTRIO' in Oname and ('-14' in Oname):
                genus = 'Calephelis'
    return genus, sp


def find_IDs(Oname):
    sampleID, CSid = 'NA', 'NA'
    #currently assume only the second column has the ID
    try:
        ID = Oname.strip().split('|')[1]
    except:
        return sampleID, CSid

    ID = parse_ID(ID)
    if ID in sequenced_set:
        sampleID = ID

    if ID in CSid_set:
        CSid = ID

    if CSid == 'NA':
        for each in Oname.strip().split('|'):
            if '-SRNP-' in each:
                CSid = each
                break

    return sampleID, CSid



def lock_database():
    conn = sql.connect(fall, timeout=2333)
    return conn


def close_database(conn):
    purge_databases_into_text(conn)
    conn.commit()
    conn.close()


def select_barcodes(condition, conn):
    c = conn.cursor()
    #print "SELECT * FROM barcode_fasta WHERE %s" % condition
    c.execute("SELECT * FROM barcode_fasta WHERE %s" % condition)
    alist = c.fetchall()
    return alist


def add_barcode(sampleID, isPCR, CSid, name, seq, conn, isCheck=True, isCommit=False):
    #when adding the barcode, I would like to
    name = parse_name(name)
    genus, sp = find_genus_and_species(name)
    if len(seq) == 698:
        seq = seq[20:678]

    if len(seq) != 658:
        print('warning: this sequence length %s is not right for %s %s' % (len(seq), name, seq))

    c = conn.cursor()
    isInDb = False
    if isCheck:
        c.execute('SELECT * FROM barcode_fasta WHERE full_name = \'%s\'' % name)
        lines = c.fetchall()
        #print lines
        if len(lines) != 0:
            isError = True
            #TODO: check if sequence is the same
            if len(lines) == 1:
                seq2 = str(lines[0][-1])
                if seq == seq2:
                    isError = False
                    isInDb = True

            if isError:
                print('Error: %s exist in database, please decide what to do (default is skip)!' % name)
                return conn
    if not isInDb:
        c.execute("INSERT INTO barcode_fasta VALUES ('%s', '%s', %s, '%s', '%s', '%s', '%s') " % (sampleID, CSid, isPCR, genus, sp, name, seq))

    if isCommit:
        conn.commit()
    return conn


def delete_barcode_by_name(fullname, conn, isCheck=True, isCommit=False):
    c = conn.cursor()
    #firstly check if it is only find one
    if isCheck:
        c.execute('SELECT * FROM barcode_fasta WHERE full_name = \'%s\'' % fullname)
        lines = c.fetchall()
        if len(lines) == 0:
            print('Error: Can not find %s in database!' % fullname)
            return conn
        else:
            if len(lines) != 1:
                print('Error: more than 1 barcode found for %s' % fullname)
                return conn
    print('deleting %s from the database...' % fullname)
    c.execute('DELETE FROM barcode_fasta WHERE full_name = \'%s\'' % fullname)
    if isCommit:
        conn.commit()
    return conn


def merge_seq(baseSeq, seq):
    new = []
    hasMerge = False
    for i, char in enumerate(baseSeq):
        if char == 'N':
            char2 = seq[i]
            if char2 != 'N':
                new.append(char2)
                hasMerge = True
            else:
                new.append(char)
        else:
            new.append(char)
    return ''.join(new), hasMerge

def add_primer(seq):
    a = 'ACTAATCATAAAGATATTGG'
    b = 'TGATTTTTTGGTCATCCAGA'

    return '%s%s%s' % (a, seq, b)


def count_gap(x):
    return x.count('N') + x.count('-')


def collapse_sp(adict):
    mapping = []
    rlist = []
    for key in adict:
        if key.strip() == '':
            continue
        alist = adict[key]
        alist = [each for each in alist
                if len(each[1]) == 658]
        if len(alist) == 0:
            continue

        if len(alist) == 1:
            name, seq = alist[0]
            rlist.append((key, seq))
            mapping.append('%s\t%s' % (key, name))

        else:
            sorted_list = sorted(alist, key=lambda x: count_gap(x[1]))
            baseSeq = sorted_list[0][1].replace('-', 'N')
            namelist = [sorted_list[0][0]]
            for name, seq in sorted_list[1:]:
                seq = seq.replace('-', 'N')
                baseSeq, hasMerge = merge_seq(baseSeq, seq)
                if hasMerge:
                    namelist.append(name)

            mapping.append('%s\t%s' % (key, ','.join(namelist)))

            rlist.append((key, baseSeq))
    fout = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/data/barcodes/species_barcode.nameMapping'
    with open(fout, 'w') as dp:
        dp.write('\n'.join(mapping))
    os.system('chmod a+w %s' % fout)
    return rlist


def purge_databases_into_text(conn):
    print('reformating dataset...')
    #TODO: if length is not 658, don't use it in bait but use it in verification set
    fverify = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/data/barcodes/all_barcodes_4verify.fa'
    fsp = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    #TODO: a mapping table
    c = conn.cursor()
    all_barcodes = []
    spDict = {}
    for row in c.execute("SELECT * from barcode_fasta"):
        genus, sp, name, seq = row[-4:]
        name = parse_name(name)
        label = '%s_%s' % (genus, sp)
        if sp == 'NA':
            label = genus
        label = parse_name(label)

        all_barcodes.append((name, seq))
        try:
            spDict[label].append((name, seq))
        except KeyError:
            spDict[label] = [(name, seq)]

    output_barcodes(all_barcodes, fverify, skipBWA=True)

    #ignore length less than 658
    sp_barcodes = collapse_sp(spDict)

    output_barcodes(sp_barcodes, fsp, addPrimer=True)



def output_barcodes(barcodes, fn, addPrimer=False, skipBWA=False):
    with open(fn, 'w') as dp:
        for name, seq in barcodes:
            name = parse_name(name)
            seq = seq.replace('-', 'N')
            if addPrimer:
                seq = add_primer(seq)
            fasta = '>%s\n%s\n' % (name, seq)
            dp.write(fasta)
    cmd = 'cd /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/data/barcodes/; makeblastdb -in=%s -dbtype=nucl;' % fn
    if not skipBWA:
        cmd += 'module add bwa; bwa index %s;' % fn
    os.system(cmd)

    os.system('chmod a+w %s*' % fn)
