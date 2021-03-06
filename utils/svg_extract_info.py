
filename = r"C:\Users\ATPs\OneDrive\Lab\UTSW\2018JunoniaProject_Xiaolong\20180918PhyloTree\20181010junoniaSeq_whole_max30Gap_split.RAxMLbestTree.sum.newick.renamed.svg"

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

print('number of labels', len(txt_label),'\n\n\n')

ls_labels = []
for e in txt_label:
    if 'transform' in e.attrs:
        ls_labels.append([e.text, float(e.attrs['transform'].split(',')[-1][:-1])])
    elif 'transform' in e.parent.attrs:
        ls_labels.append([e.text, float(e.parent.attrs['transform'].split(',')[-1][:-1])])
ls_labels.sort(key=lambda x:x[1])
ls_samples = ['S'+e[0].split('_')[0] for e in ls_labels]
print('\n'.join(ls_samples))

open(r"C:\Users\ATPs\Downloads\temp.txt",'w',newline='\n').write('\n'.join(e[1:] for e in ls_samples))

