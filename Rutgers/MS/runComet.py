import os
import numpy as np
# protein_database = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/ProteinDB/comet/20200529contaminants.UniprotRefSeqGencodeHLAChess.StringTie.fa.idx'
# #protein_database2 = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/ProteinDB/comet/20200529contaminants.UniprotRefSeqGencodeHLAChess.StringTie.fa.DECOY.idx'
# parameter_file = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/ProteinDB/comet.params.NoDecoySearch.Trypsin.HCD'
# file_mzml = '/scratch/xc278/mzml/MSV000078777/SulfenM_RKO_LCA_A01.mzML.gz'
# folder_out = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/MSMS/MassIVE-KB/comet'
# prefix = 'test'

def runComet(file_mzml, protein_database, parameter_file, folder_out = None, prefix=None):
    '''
    run comet for the file_mzml with protein_database and parameter_file
    file_mzml is location of one file.
    '''
    # determine file_mzml type, if is .gz file
    if file_mzml.endswith('.gz'):
        file_type = 'gz'
    else:
        file_type = 'mzML'
    
    # output filename
    if prefix is None:
        prefix = os.path.basename(file_mzml)
        if prefix.endswith('.mzML.gz') or prefix.endswith('.mzml.gz'):
            prefix = prefix[:-8]
        else:
            prefix = prefix[:-5]
    
    # prepare output
    if folder_out is None:
        folder_out = os.getcwd()
    elif not os.path.exists(folder_out):
        os.makedirs(folder_out)
    
    # create a softlink of file_mzml to output folder
    if file_type == 'gz':
        file_mzml2 = os.path.join(folder_out, prefix + '.mzML.gz')
    else:
        file_mzml2 = os.path.join(folder_out, prefix + '.mzML')
    if not os.path.exists(file_mzml2):
        os.system('ln -s {} {}'.format(file_mzml, file_mzml2))
    
    if os.path.exists(f'{folder_out}/{prefix}.pin.gz'):
        print(file_mzml, 'already finished.')
        return True
    
    # move protein_database to /mnt/scratch if possible
    useTemp = 2 # 0, do not use temp folder. 1, use temp folder. 2, use temp folder on percentage_useTemp cases
    percentage_useTemp = 0.2
    if useTemp == 2:
        useTemp = np.random.random()
        if useTemp < percentage_useTemp:
            useTemp = 1
    if useTemp == 1:
        cashe_folder = '/scratch/'
        if os.path.exists(cashe_folder):
            stat = os.statvfs(cashe_folder)
            # check if there is over 60GB space in /mnt/scratch
            if stat.f_bfree*stat.f_bsize /1e9 > 60:
                if not os.path.exists(cashe_folder + 'xc278'):
                    os.makedirs(cashe_folder + 'xc278')
                protein_database_basename = os.path.basename(protein_database)
                protein_database_tmp = os.path.join(cashe_folder + 'xc278', protein_database_basename)
                print('will use {} folder to speed up'.format(cashe_folder + 'xc278'))
                if not os.path.exists(protein_database_tmp):
                    os.system(f'cp {protein_database} {protein_database_tmp}')
                    protein_database = protein_database_tmp
                elif os.path.getsize(protein_database) == os.path.getsize(protein_database_tmp):
                    protein_database = protein_database_tmp
        os.system(f'cat {protein_database} > /dev/null')


    # run comet
    cmd = f'cd "{folder_out}" && comet -D{protein_database} -P{parameter_file} "{file_mzml2}" && gzip "{prefix}.pin"'
    print(cmd)
    os.system(cmd)
    
    # clean up
    if os.path.islink(file_mzml2):
        os.remove(file_mzml2)

description = '''
    run comet with the settings
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--file_mzml', help = 'mzML or mzML.gz file', required=True)
    parser.add_argument('-d','--protein_database',help = 'protein database', required=True)
    parser.add_argument('-r','--parameter_file',help='comet parameter file', required=True)
    parser.add_argument('-f', '--folder_out', help='whether to put output file. default None, output to current folder', default=None)
    parser.add_argument('-p', '--prefix', help='prefix for ouput file. default, the same as the mzML file', default=None)

    f = parser.parse_args()
    runComet(file_mzml = f.file_mzml, protein_database = f.protein_database, parameter_file = f.parameter_file, folder_out = f.folder_out, prefix=f.prefix)

