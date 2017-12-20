
# coding: utf-8

# This is a notebok that drafts a small script to pull in all the data for the WGS pathogen detection and microbiome data. It needs to pull in the following:
#     * sequence summary file from albacore
#     * two blast output files (specific database and NCBI)
#     * get the length of all porchoped reads
#     * get processing (T/F) for porchoped and nanolyze
#     * get this all into one data frame
#     
#     
#     

# In[1]:

import os
import pandas as pd
import subprocess as sub
from Bio import SeqIO
import argparse


# In[ ]:

parser = argparse.ArgumentParser(description='''
This is a notebok that drafts a small script to pull in all the data for the WGS pathogen detection and microbiome data. It needs to pull in the following:

    * sequence summary file from albacore
    * two blast output files (specific database and NCBI)
    * get the length of all porchoped reads
    * get processing (T/F) for porchoped and nanolyze
    * get this all into one data frame
''')

parser.add_argument("BASEDIR", help="base folder, supposed to have all the sub folders processed by WGS script")


# In[ ]:

args = parser.parse_args()


# In[2]:

#lets define the base folder
#BASEDIR = '/home/yiheng/data/Wagga_run1'
BASEDIR = args.BASEDIR
#this will become the only flag of argparse


# In[3]:

#write a quick check that looks for all the right folders
folder_list = 'basecalled_data  scripts  tracking  workspace'.split(' ')
for x in range(0,folder_list.count('')):
    folder_list.remove('')
#fix this test
if not set(os.listdir(os.path.abspath(BASEDIR))) == set (folder_list):
    'Print something=asdfasdf'


# In[4]:

seq_df_headers = ['read_id','passes_filtering', 'sequence_length_template', 'mean_qscore_template',                  'barcode_arrangement', 'barcode_score', 'kit', 'variant']


# In[5]:

#now get the sequencing_summary file
base_called_folder = os.path.join(BASEDIR, 'basecalled_data')
for thing in os.listdir(base_called_folder):
    if not os.path.isdir(os.path.join(base_called_folder, thing)):
        next
    else:
        seq_sum_file = os.path.join(base_called_folder, thing, 'sequencing_summary.txt')
        if not os.path.exists(seq_sum_file):
            print('No sequencing summary file. Please go check')
            continue
        seq_df = pd.read_csv(seq_sum_file, sep='\t')
        #capture the thing as the prefix of the fastq/fasta files in the barcode folders
        prefix = thing
        #might be a better way to only read in the wanted columns. Not subsetting afterwards.
        #please go check
        seq_df = seq_df.loc[:, seq_df_headers].copy()

        


# In[6]:

#now get all the rgblast_output databases done 
rg_blast_df_file_list = []
nt_blast_df_file_list = []
workspace = os.path.join(BASEDIR, 'workspace')
folder_counter = 0
for folder in os.listdir(workspace):
    folder = os.path.join(workspace,folder)
    if not os.path.isdir(folder):
        next
    folder_counter += 1
    for file in os.listdir(folder):
        if not file.endswith('.rgblast_output') or not file.endswith('.ntblast_output'):
            next
        if file.endswith('.rgblast_output'):
            rg_blast_df_file_list.append(os.path.join(folder, file))
        if file.endswith('.ntblast_output'):
            nt_blast_df_file_list.append(os.path.join(folder, file))
    if not len(nt_blast_df_file_list) == len(rg_blast_df_file_list) == folder_counter:
        print('Not all barcode folders have all blast output files')
    #print(folder)


# In[7]:

#now get all the names of the reads that survived the lyzing
workspace = os.path.join(BASEDIR, 'workspace')
folder_counter = 0
nl_readid_list = []
for folder in os.listdir(workspace):
    folder_long = os.path.join(workspace,folder)
    if not os.path.isdir(folder_long):
        next
    folder_counter += 1
    for file in os.listdir(folder_long):
        lyzed_fastq = (os.path.join(folder_long, '%s.%s.fastq' % (prefix, folder )))
        if not os.path.exists(lyzed_fastq):
            print("No lyzed fastq file present.")
        elif file == '%s.%s.fastq' % (prefix, folder ):
            print(lyzed_fastq)
            lyzed_fastq_name = '%s.name.tmp' % lyzed_fastq
            #cmd = "grep '^@' %s | cut -c1-37 > %s"%\
            # (lyzed_fastq, lyzed_fastq_name)
            #cmd = "grep '^@' %s | cut -d:-f1 > %s"%\
             #(lyzed_fastq, lyzed_fastq_name)
            cmd = r"grep '^@' %s | sed -e 's/@\([a-zA-Z0-9-]\+\).*/\1/' > %s"%            (lyzed_fastq, lyzed_fastq_name)
            #print(cmd)
            sub.check_output(cmd, shell=True, stderr=sub.STDOUT)
            #sub.run(cmd.split(' '))
            nl_readid_list.append(lyzed_fastq_name)
        else:
            next


# In[8]:

#now read in the tmp file of nanolyzed ids and add them to the dataframe as column
def get_read_ids(filename_list):
    #write a check that both _list and header is a list
    read_id_list = []
    for fn in filename_list:
        tmp_list = pd.read_csv(fn, sep='\t', header=None).loc[:,0].tolist()
        read_id_list= read_id_list + tmp_list
    return read_id_list


# In[9]:

nl_survied_list = get_read_ids(nl_readid_list)


# In[154]:

#do same for porechop + length
#this is pretty slow. May consider parallzing.
pc_length_dict = {}
porechop_survived_list = []
folder_counter = 0
for folder in os.listdir(workspace):
    folder_long = os.path.join(workspace,folder)
    if not os.path.isdir(folder_long):
        next
    folder_counter += 1
    porechoped_file = os.path.join(folder_long,'%s.chopped.%s.fastq'%(prefix, folder))
    if not os.path.exists(porechoped_file):
            print("Porechopped fastq missing for %s." % folder)
    else:
        for seq in SeqIO.parse(porechoped_file, 'fastq'):
            pc_length_dict[seq.id] = len(seq.seq)
            porechop_survived_list.append(seq.id)
    
#take for loop from above to loop over chopped files
#filename = '/home/yiheng/data/Wagga_run1/workspace/barcode01/Wagga_run1_albacore202.chopped.barcode01.fastq'
#for seq in SeqIO.parse(filename, 'fastq'):
#    pc_length_dict[seq.id] = len(seq.seq)


# In[155]:

porechop_survived_single_list = [x.split('_')[0] for x in porechop_survived_list]


# In[156]:

seq_df['pc_survived'] = seq_df['read_id'].isin(porechop_survived_single_list)


# In[157]:

seq_df.head()


# In[158]:

#add the nanolyze survived column
seq_df['nl_survived'] = seq_df['read_id'].isin(nl_survied_list)
#make porchoped survived column


# In[159]:

import numpy as np


# In[179]:

blast_header = ['qseqid',
 'sseqid',
 'evalue',
 'bitscore',
 'length',
 'pident',
 'nident',
 'sgi',
 'sacc',
 'staxids',
 'scomnames',
'sskingdoms']
# original headers: qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms
def make_all_blast_df(_list, header, chopped_len_dict):
    df = pd.DataFrame()
    #write a check that both _list and header is a list
    for x in _list:
        tmp_df = pd.read_csv(x, sep='\t',header=None, names=header)
        first_column = tmp_df.columns[0]
        tmp_df['read_id'] = tmp_df[first_column].apply(lambda x: str(x).split('_')[0])
        tmp_df['read_length_pc'] = tmp_df[first_column].apply(lambda x: chopped_len_dict[x])
        df = pd.concat([df, tmp_df.iloc[:,[0,1,2,4,5,6,8,9,10,12,13]]])
        
    #now reduce the columns to what we want
    
    print(first_column)
    #df['read_id'] = df.iloc[:,0].str.split('_')[0]
    return df


# In[180]:

nt_df = make_all_blast_df(nt_blast_df_file_list, [x +'_nt' for x in  blast_header], pc_length_dict)
rg_df = make_all_blast_df(rg_blast_df_file_list, [x +'_rg' for x in  blast_header], pc_length_dict)
#reduce column number of blast dataframe to what you want before you merge


# In[181]:

#need to take care of porchop split reads checkin if second last character of string is _
#make a new column of the blast_df that has the initial read_id
#need to check how merge behaves when getting a doublcate of value in one df.
tmp_df = pd.merge(seq_df, rg_df,how='outer',left_on= 'read_id', right_on='read_id')
final_df = pd.merge(tmp_df, nt_df,how='outer',left_on= 'read_id', right_on='read_id')


# In[187]:

final_df.iloc[:200,0:25]


# In[ ]:




# In[ ]:

def length_column(x):
    if x in pc_length_dict.keys():
        return pc_length_dict[x]
    else:
        return np.nan


# In[ ]:

#porchop length column once pc_length_dict is done
seq_df['pc_length'] = seq_df.read_id.apply(lambda x: length_column(x))


# In[ ]:

seq_df['pc_length'].unique()


# In[ ]:

#remove tmp files again
pc_length_dict


# In[ ]:




# In[ ]:




# In[ ]:




# In[127]:

#need to take care of porchop split reads checkin if second last character of string is _
#make a new column of the blast_df that has the initial read_id
#need to check how merge behaves when getting a doublcate of value in one df.
test_df = pd.merge(seq_df, rg_df,how='outer',left_on= 'read_id', right_on='read_id')


# In[139]:

test_2 = test_df[test_df['qseqid_rg'].str.contains('_') == True]


# In[144]:

test_df.iloc[[992299,992300],:]


# In[143]:

test_2.loc[:,['read_id', 'qseqid_rg']]


# In[ ]:

test_df[test_df.read_id=='f7e483ca-c6b8-46be-975f-3bc433b3b0f6']


# In[ ]:

test_df.loc[:, ['read_id', 'qseqid_rg']]


# In[ ]:

rg_df.columns


# In[ ]:

#now get the sequencing_summary file
base_called_folder = os.path.join(BASEDIR, 'basecalled_data')
for thing in os.listdir(base_called_folder):
    if not os.path.isdir(os.path.join(base_called_folder, thing)):
        next
    else:
        seq_sum_file = os.path.join(base_called_folder, thing, 'sequencing_summary.txt')
        if not os.path.exists(seq_sum_file):
            print('No sequencing summary file. Please go check')
        seq_df = pd.read_csv(seq_sum_file, sep='\t')
        #might be a better way to only read in the wanted columns. Not subsetting afterwards.
        #please go check
        seq_df = seq_df.loc[:, seq_df_headers].copy()


# In[ ]:




# In[ ]:



