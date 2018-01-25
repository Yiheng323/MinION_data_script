
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

# In[59]:

import os
import pandas as pd
import subprocess as sub
from Bio import SeqIO
import argparse


# In[60]:

parser = argparse.ArgumentParser(description='''
This is a notebok that drafts a small script to pull in all the data for the WGS pathogen detection and microbiome data. It needs to pull in the following:

    * sequence summary file from albacore
    * two blast output files (specific database and NCBI)
    * get the length of all porchoped reads
    * get processing (T/F) for porchoped and nanolyze
    * get this all into one data frame
''')

parser.add_argument("BASEDIR", help="base folder, supposed to have all the sub folders processed by WGS script. The same as Indir in YH_script2")


# In[61]:

args = parser.parse_args()


# In[62]:

#lets define the base folder
#BASEDIR = '/home/yiheng/test/20170617_FAH05731'
BASEDIR = args.BASEDIR
#this will become the only flag of argparse


# In[83]:

#write a quick check that looks for all the right folders
folder_list = 'basecalled_data  scripts  tracking  workspace'.split(' ')
for x in range(0,folder_list.count('')):
    folder_list.remove('')
#fix this test
if not set(os.listdir(os.path.abspath(BASEDIR))) >= set (folder_list):
    print("Something wrong with basefolder. check it please.")


# In[64]:

# define the columns that you want to pick up from sequencing summary file.
# Here is the columns I chose for plotting out data, enough information for me so I did not pick others.
seq_df_headers = ['read_id','passes_filtering', 'sequence_length_template', 'mean_qscore_template',                  'barcode_arrangement', 'barcode_score', 'kit', 'variant']


# In[65]:

# add a function to check the tar.gz file and its corresponding unzipped folder in the basecalled_data folder. 
# If it is not unzipped, unzip it.


# In[66]:

#now get the headers sequencing_summary file
base_called_folder = os.path.join(BASEDIR, 'basecalled_data')

# here added a function to check the tar.gz file and its corresponding unzipped folder in the basecalled_data folder.  
# If it is not unzipped, unzip it.
for thing in os.listdir(base_called_folder):
    judge_list = [os.path.isdir(os.path.join(base_called_folder, thing))]

if any(x == True for x in judge_list):
    seq_sum_file = os.path.join(base_called_folder, thing, 'sequencing_summary.txt')
    if not os.path.exists(seq_sum_file):
        print('No sequencing summary file from basecalled folder. Please go check')
    
    print("now getting the headers from %s" % seq_sum_file)
    seq_df = pd.read_csv(seq_sum_file, sep='\t')
    #capture the thing as the prefix of the fastq/fasta files in the barcode folders
    prefix = thing
    #might be a better way to only read in the wanted columns. Not subsetting afterwards.
    #please go check
    seq_df = seq_df.loc[:, seq_df_headers].copy()
    
elif all(x == False for x in judge_list):
        
    zipped_basecalled_file = os.path.join(base_called_folder, thing)
    if zipped_basecalled_file.endswith('tar.gz'):
        print("now unzipping file %s." % zipped_basecalled_file)
        tar = tarfile.open(zipped_basecalled_file)
        tar.extractall(base_called_folder.split(".")[0])
        tar.close()
    else:
        print("there is something strange in the basecalled folder, please check.")
else:
        print("there is something strange in the basecalled folder, please check.")

        


# In[67]:

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
    else:
        next
    #print(folder)


# In[68]:

#now get all the names of the reads that survived the lyzing
workspace = os.path.join(BASEDIR, 'workspace')
folder_counter = 0
nanolyze_readid_list = []
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
            print("now processing nanolyzed file %s" % lyzed_fastq)
            lyzed_fastq_name = '%s.name.tmp' % lyzed_fastq
            #cmd = "grep '^@' %s | cut -c1-37 > %s"%\
            # (lyzed_fastq, lyzed_fastq_name)
            #cmd = "grep '^@' %s | cut -d:-f1 > %s"%\
             #(lyzed_fastq, lyzed_fastq_name)
            cmd = r"grep '^@' %s | sed -e 's/@\([a-zA-Z0-9-]\+\).*/\1/' > %s"%            (lyzed_fastq, lyzed_fastq_name)
            #print(cmd)
            sub.check_output(cmd, shell=True, stderr=sub.STDOUT)
            #sub.run(cmd.split(' '))
            nanolyze_readid_list.append(lyzed_fastq_name)
        else:
            next


# In[69]:

#now read in the tmp file of nanolyzed ids and add them to the dataframe as column
def get_read_ids(filename_list):
    #write a check that both _list and header is a list
    read_id_list = []
    for filename in filename_list:
        tmp_list = pd.read_csv(filename, sep='\t', header=None).loc[:,0].tolist()
        read_id_list= read_id_list + tmp_list
    return read_id_list


# In[70]:

nl_survied_list = get_read_ids(nanolyze_readid_list)


# In[71]:

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
    print("now processing porechopped file %s." % porechoped_file)
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


# In[72]:

# some of the porechopped reads will be splited into two reads they will have an _ in their read id
porechop_survived_single_list = [x.split('_')[0] for x in porechop_survived_list]


# In[73]:

# add porechop survived column
seq_df['pc_survived'] = seq_df['read_id'].isin(porechop_survived_single_list)


# In[74]:

seq_df.to_csv


# In[75]:

#add the nanolyze survived column
seq_df['nl_survived'] = seq_df['read_id'].isin(nl_survied_list)
#make porchoped survived column


# In[76]:

import numpy as np


# In[77]:

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
    
    #print(first_column)
    #df['read_id'] = df.iloc[:,0].str.split('_')[0]
    return df


# In[78]:

nt_df = make_all_blast_df(nt_blast_df_file_list, [x +'_nt' for x in  blast_header], pc_length_dict)
print("now adding the nt output columns.")
rg_df = make_all_blast_df(rg_blast_df_file_list, [x +'_rg' for x in  blast_header], pc_length_dict)
print("now adding the rg output columns.")
#reduce column number of blast dataframe to what you want before you merge


# In[79]:

#need to take care of porchop split reads checkin if second last character of string is _
#make a new column of the blast_df that has the initial read_id
#need to check how merge behaves when getting a doublcate of value in one df.
tmp_df = pd.merge(seq_df, rg_df,how='outer',left_on= 'read_id', right_on='read_id')
final_df = pd.merge(tmp_df, nt_df,how='outer',left_on= 'read_id', right_on='read_id')
print("now creating the final dataframe!")


# In[80]:

analysis_foler = os.path.join(BASEDIR, 'analysis')
if not os.path.exists(analysis_foler):
    os.mkdir(analysis_foler)
final_df_fn = os.path.join(analysis_foler, 'summary_df.tab')
final_df.to_csv(final_df_fn, sep='\t', index=None)
print("All done. Congratulations!")


# In[81]:

#final_df.iloc[:200,0:25]


# In[ ]:




# In[24]:

#def length_column(x):
#    if x in pc_length_dict.keys():
#        return pc_length_dict[x]
#    else:
#        return np.nan


# In[25]:

#porchop length column once pc_length_dict is done
#seq_df['pc_length'] = seq_df.read_id.apply(lambda x: length_column(x))


# In[ ]:



