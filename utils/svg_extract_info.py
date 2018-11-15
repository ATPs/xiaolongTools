
filename = r"C:\Users\K\OneDrive\Lab\UTSW\2018JunoniaProject_Xiaolong\20180918PhyloTree\20181010junoniaSeq_CDS_split.RAxMLbestTree.sum.newick.renamed.svg"

from bs4 import BeautifulSoup
soup = BeautifulSoup(open(filename).read(),'lxml')
txts = soup.find_all(name='text')

def isfloat(element):
    '''
    test if element can be converted to float
    '''
    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False

txt_label = [e for e in txts if not isfloat(e.text)]

ls_labels = [[e.text,float(e.attrs['transform'].split(',')[-1][:-1])] for e in txt_label]
ls_labels.sort(key=lambda x:x[1])
ls_samples = ['S'+e[0].split('_')[0] for e in ls_labels]
print('\n'.join(ls_samples))