file1 = r"C:\Users\ATPs\Google Drive\Lepidoptera_projects\All-plates.xlsx"
file2 = r"C:\Users\ATPs\Google Drive\Lepidoptera_projects\BugsInRNAlater-2015-04-11.xlsx"


import pandas as pd
df1 = pd.read_excel(file1, sheet_name='Sheet1', dtype=str,encoding = 'utf-8')
df2 = pd.read_excel(file2, sheet_name='Sheet1', dtype=str,encoding = 'utf-8')

df1.index = df1['DNA number'].apply(lambda x:x.replace('NVG-',''))
df2.index = df2['DNA Number'].apply(lambda x:x.replace('NVG-',''))

df1 = df1[['Taxon name', 'Type status', 'Sex', 'Country','State/Province', 'Date']]
df1.columns = ['Taxon name', 'Type status', 'Sex', 'Country','State', 'Date']
df2 = df2[['Taxon name', 'Type status', 'Sex', 'Country','State', 'Date']]
df = pd.concat([df1,df2])

fout = open(r"C:\Users\ATPs\Google Drive\Lepidoptera_projects\all_samples_name_tab.txt",'w',newline="\n")
for sample_id, row in df.iterrows():
    txt = (sample_id + '\t'+'_'.join(row).replace('_nan','')+'\n').replace(' ','_')
    txt = txt.replace('_00:00:00','')
    fout.write(txt)
fout.close()
