
# coding: utf-8

# In[16]:

import argparse
import os
import subprocess
import datetime


# In[17]:

parser = argparse.ArgumentParser(description="This is a script used for WGS data analysis from the field infected wheat samples")



# In[18]:

parser.add_argument("threads", type=int, help="number of threads")


# In[19]:

parser.add_argument("ram", help="expected ram usage per each threads")


# In[20]:

parser.add_argument("sep_barcodes", help="barcodes that will be separately put in individual folder. separate by ','")


# In[21]:

parser.add_argument("indir", help="folder storing all input data")


# In[22]:

parser.add_argument("blastdb_1", help="database for the first blast")


# In[23]:

parser.add_argument("blastdb_2", help="database for the second blast")


# In[ ]:




# In[24]:

args = parser.parse_args()


# In[25]:

t = args.threads
r = args.ram
s = args.sep_barcodes
d = args.indir
db1 = args.blastdb_1
db2 = args.blastdb_2


# This will be a dummy notebook for generating the analysis scripts for Yihengs WGS nanopore data analysis

# In[51]:

os.getcwd()


# In[52]:

#stuff you get from argparse
print ("INDIR = {}".format(d))
print ("sep_barcodes = {}".format(s))
print ("threads = {}".format(t))
print ("ram = {}".format(r))
print ("BLAST_DB_1 = {}".format(db1))
#print ("BLAST_DB_2 = {}".format(db2)
#BLAST_DB_2 = 'nt'



# In[53]:

#get the date
now = datetime.datetime.now()
date = now.strftime("%Y%m%d")


# In[54]:

#define the path for control DNA sequence for nanolyse
CDS_PATH = '/home/yiheng/bio/ncbi/db/DNA_CS.fasta.gz'


# In[55]:

date


# In[56]:

if not os.path.exists(INDIR):
    print('INDIR %s does not exsist. Exit!' % INDIR)
    #exit()


# In[57]:

os.listdir(INDIR)


# In[58]:

os.path.abspath(INDIR)


# In[59]:

SCRIPT_FOLDER = os.path.join(os.path.abspath(INDIR), 'scripts')
WORKSPACE = os.path.join(os.path.abspath(INDIR), 'workspace')
TRACKING = os.path.join(os.path.abspath(INDIR), 'tracking')
if not os.path.exists(SCRIPT_FOLDER):
    os.mkdir(SCRIPT_FOLDER)
if not os.path.exists(WORKSPACE):
    os.mkdir(WORKSPACE)
if not os.path.exists(TRACKING):
    os.mkdir(TRACKING)


# In[60]:

BASECALLFOLDER = os.path.join(os.path.abspath(INDIR),'basecalled_data' )


# In[61]:

#not check again if the basecalled data is there
if not os.path.exists(BASECALLFOLDER):
    print('No basecalled data in %s. Exit!' % INDIR)
    #exit()


# In[62]:

os.listdir(BASECALLFOLDER)


# In[63]:

if  len([os.path.join(BASECALLFOLDER, x) for x in os.listdir(BASECALLFOLDER) if x.endswith('tar.gz')]) == 1:
    tar_file = [os.path.join(BASECALLFOLDER, x) for x in os.listdir(BASECALLFOLDER) if x.endswith('tar.gz')][0]
    runid = tar_file.split('/')[-1].split('.')[0]
else:
        print('None or mulitiple tar files')


# In[64]:

#now untar the tar gz file
os.chdir(BASECALLFOLDER)
unzip_command = 'tar -xvf %s' %(tar_file)
print(unzip_command)
unzip_command_stderr = subprocess.check_output(unzip_command, shell=True, stderr=subprocess.STDOUT)


# In[65]:

#now we should have generated the folder with the basecallded data
BASECALLED_DATA_FOLDER = os.path.join(BASECALLFOLDER, tar_file.split('.')[0])
if not os.path.exists(BASECALLED_DATA_FOLDER):
    print("Something with the unzipping of the basecalled data went wrong.")


# In[66]:

#check the contents in the zipped folder
content = ['workspace', 'configuration.cfg', 'sequencing_summary.txt', 'pipeline.log']
if not set(os.listdir(BASECALLED_DATA_FOLDER)) == set(content):
    print("Something with the unzipping of the basecalled data went wrong.")

BASECALLED_DATA_WORKSPACE = os.path.join(BASECALLED_DATA_FOLDER, 'workspace')


# In[67]:

PASS_FOLDER = os.path.join(BASECALLED_DATA_WORKSPACE, 'pass')
os.listdir(PASS_FOLDER)


# In[68]:

FAIL_FOLDER = os.path.join(BASECALLED_DATA_WORKSPACE, 'fail')
os.listdir(FAIL_FOLDER)


# In[69]:

sep_barcodes_list = sep_barcodes.split(',')


# In[70]:

# sep_barcodes_list

# make a barcode 00 folder
# loop throug all folders that are not in the sep_barcode list and combine those
# add 00 to the sep_barcode list

# sep_barcodes_list
# In[71]:

#make the barcode00 folder to combine all the misidentified/unclassified barcode
if not os.path.exists(os.path.join(FAIL_FOLDER, 'barcode00')):
    os.mkdir(os.path.join(FAIL_FOLDER, 'barcode00'))

if not os.path.exists(os.path.join(PASS_FOLDER, 'barcode00')):
    os.mkdir(os.path.join(PASS_FOLDER, 'barcode00'))


# In[72]:

['barcode%s' % x for x in sep_barcodes_list]


# In[73]:

#cat all barcode misclassified/unclassified reads in the pass folder to one fastq file and redirect to barcode00 folder

if not os.path.isfile('%s/barcode00.fastq' % os.path.join(PASS_FOLDER, 'barcode00')):
    cat_fastq_command = 'cat'
    for folder in [y.split('/')[-1] for y in os.listdir(PASS_FOLDER)]:
        if folder not in ['barcode%s' % x for x in sep_barcodes_list]:
            tmp_folder = os.path.join(PASS_FOLDER, folder)
            #print(folder)
            for fastq in [os.path.join(tmp_folder, z) for z in os.listdir(tmp_folder)  if z.endswith('fastq')]:
                cat_fastq_command += ' %s' % fastq
                #print(fastq)
    cat_fastq_command +=' > %s/barcode00.fastq' % os.path.join(PASS_FOLDER, 'barcode00')
    cat_fastq_command_stderr = subprocess.check_output(cat_fastq_command, shell=True, stderr=subprocess.STDOUT)


# In[74]:

#cat all barcode misclassified/unclassified reads in the fail folder to one fastq file and redirect to barcode00 folder
if not os.path.isfile('%s/barcode00.fastq' % os.path.join(FAIL_FOLDER, 'barcode00')):
    cat_fastq_command = 'cat'
    for folder in [y.split('/')[-1] for y in os.listdir(FAIL_FOLDER)]:
        if folder not in ['barcode%s' % x for x in sep_barcodes_list]:
            tmp_folder = os.path.join(FAIL_FOLDER, folder)
            #print(folder)
            for fastq in [os.path.join(tmp_folder, z) for z in os.listdir(tmp_folder)  if z.endswith('fastq')]:
                cat_fastq_command += ' %s' % fastq
                #print(fastq)
    cat_fastq_command +=' > %s/barcode00.fastq' % os.path.join(FAIL_FOLDER, 'barcode00')
    cat_fastq_command_stderr = subprocess.check_output(cat_fastq_command, shell=True, stderr=subprocess.STDOUT)


# In[75]:

#elements of script
def print_begining_of_script(barcode, fn):
    """
    Takes the barcode and the file name to start off the script
    """
    begining_of_script = "#!/bin/bash\n\
#$ -M yiheng.hu@anu.edu.au\n\
#$ -m a\n\
#$ -cwd\n\
#$ -V\n\
#$ -j y\n\
#$ -pe threads %s\n\
#$ -l h_vmem=%s,virtual_free=%s\n\
#$ -N barcode_%s\n\
set -vx\
" %(threads, ram, ram, barcode)
    print(begining_of_script, file=fn)


# In[76]:

#now get the nanolyse going
def print_nanolyse(barcode, tar_file, fn):
    """Fill in something"""
    fastq_prefix = runid + '.barcode%s' % (barcode)
    nanolyse_step = "\n#now starts the combining of pass/fail data and nanolyze step\n"
    nanolyse_step += "cd %s\n" % BASECALLED_DATA_WORKSPACE
    BARCODE_WS = os.path.join(WORKSPACE, 'barcode%s' %barcode)
    nanolyse_step += "mkdir -p %s\n" % BARCODE_WS
    nanolyse_step += "cat pass/barcode%s/*.fastq fail/barcode%s/*.fastq > %s/%s.unlysed.fastq\n"    % (barcode, barcode, BARCODE_WS, fastq_prefix )
    nanolyse_step += "gzip %s/%s.unlysed.fastq\n" % (BARCODE_WS, fastq_prefix)
    nanolyse_step += "gunzip -c %s/%s.unlysed.fastq" % (BARCODE_WS, fastq_prefix)
    nanolyse_step += " | NanoLyse --reference %s " % (CDS_PATH)
    nanolyse_step += "| gzip > %s/%s.fastq.gz\n" % (BARCODE_WS, fastq_prefix)
    nanolyse_step += "gunzip %s/%s.fastq.gz\n" % (BARCODE_WS, fastq_prefix)


    print(nanolyse_step, file=fn)


# In[77]:

# def print_porechop(barcode, folder, runid, fn):
#    """
#    Fill in something
#    """
#    fastq = '%s.barcode%s.fastq'% (runid, barcode)
#    chopped_fastq = '%s.chopped.barcode%s.fastq'% (runid, barcode)
#    porechop_step ="\n#now the porechop step\n"
#    porechop_step += 'porechop -i %s/%s -o %s/%s --format fastq --middle_threshold 95'\
#    %(folder,fastq, folder, chopped_fastq )
#    print(porechop_step, file=fn)


# In[78]:

def print_porechop(barcode, folder, runid, fn):
    """
    Fill in something
    """
    fastq = '%s.barcode%s.fastq'% (runid, barcode)
    chopped_fastq = '%s.chopped.barcode%s.fastq'% (runid, barcode)
    porechop_step ="\n#now the porechop step\n"
    porechop_step += 'porechop -i %s/%s -o %s/%s --format fastq --middle_threshold 95\n'    %(folder,fastq, folder, chopped_fastq )
    porechop_step += "seqtk seq -a %s/%s > %s/%s\n"     %(folder, chopped_fastq, folder, chopped_fastq.replace('.fastq', '.fasta'))
    print(porechop_step, file=fn)


# In[79]:

def print_blastn_1(barcode, folder, runid, fn):
    """Fill in something here."""
    chopped_fasta = '%s.chopped.barcode%s.fasta'% (runid, barcode)
    hit_ids_fn = "%s/%s.%sblast.qseqid.barcode%s.txt" %(folder, runid, BLAST_DB_1, barcode)

    blastn_step = "\n#now the blast step against DB %s\n" % BLAST_DB_1
    blastn_step += "blastn -query %s/%s -db %s -evalue 0.01 -outfmt " %(folder,chopped_fasta, BLAST_DB_1 )
    blastn_step += "'6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms'"
    blastn_step += "-show_gis -num_threads %s | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge " % (threads)
    blastn_step += "> %s/%s.%s.%sblast_output\n\n" %(folder,chopped_fasta,date ,BLAST_DB_1)
    blastn_step += "cut -f 1 %s/%s.%s.%sblast_output" %(folder,chopped_fasta,date ,BLAST_DB_1)
    blastn_step += " > %s\n\n" % hit_ids_fn
    blastn_step += "filterbyname.sh in=%s/%s " % (folder, chopped_fasta)
    blastn_step += "out=%s/%s " % (folder, chopped_fasta.replace('.barcode%s.'%barcode,                                                              '.%shityes.barcode%s.'%(BLAST_DB_1, barcode)))
    blastn_step += "names=%s include=t\n" % hit_ids_fn
    blastn_step += "filterbyname.sh in=%s/%s " % (folder, chopped_fasta)
    blastn_step += "out=%s/%s " % (folder, chopped_fasta.replace('.barcode%s.'%barcode,                                                              '.%shitno.barcode%s.'%(BLAST_DB_1, barcode)))
    blastn_step += "names=%s include=f\n" % hit_ids_fn

    print(blastn_step, file =fn)


# In[80]:

#def print_blastn_2(barcode, folder, runid, fn):
#    """Fill in something here."""
#    nohit_fasta = '%s.chopped.%shitno.barcode%s.fasta'% (runid, BLAST_DB_1, barcode)
#    hit_ids_fn = "%s/%s.%shitno.%sblast.qseqid.barcode%s.txt" %(folder, runid, BLAST_DB_1, BLAST_DB_2, barcode)
#
#    blastn_step = "\n#now the blast step against database %s\n" % BLAST_DB_2
#    blastn_step += "blastn -query %s/%s -db %s -evalue 0.01 -outfmt " %(folder,nohit_fasta, BLAST_DB_2 )
#    blastn_step += "'6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms'"
#    blastn_step += "-show_gis -num_threads %s | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge " % (threads)
#    blastn_step += "> %s/%s.%s.%sblast_output\n\n" %(folder,nohit_fasta,date,BLAST_DB_2)
#    blastn_step += "cut -f 1 %s/%s.%s.%sblast_output" %(folder,nohit_fasta,date,BLAST_DB_2)
#    blastn_step += " > %s\n\n" % hit_ids_fn
#    blastn_step += "filterbyname.sh in=%s/%s " % (folder, nohit_fasta)
#    blastn_step += "out=%s/%s " % (folder, nohit_fasta.replace('.barcode%s.'%barcode,\
#                                                              '.%shityes.barcode%s.'%(BLAST_DB_2, barcode)))
#    blastn_step += "names=%s include=t\n" % hit_ids_fn
#    blastn_step += "filterbyname.sh in=%s/%s " % (folder, nohit_fasta)
#    blastn_step += "out=%s/%s " % (folder, nohit_fasta.replace('.barcode%s.'%barcode,\
#                                                              '.%shitno.barcode%s.'%( BLAST_DB_2, barcode)))
#    blastn_step += "names=%s include=f\n" % hit_ids_fn

#    print(blastn_step, file =fn)


# In[81]:

#add barcode 00
sep_barcodes_list.append('00')


# In[82]:

sep_barcodes_list


# In[85]:

for barcode in sep_barcodes_list:
    tmp_file = os.path.join(SCRIPT_FOLDER, 'barcode%s.sh' % barcode)
    with open(tmp_file, 'w') as fn:
        TMP_BARCODE_FOLDER = os.path.join(WORKSPACE, 'barcode%s' % barcode)

        print_begining_of_script(barcode, fn)
        print_nanolyse(barcode, runid, fn)
        print_porechop(barcode, TMP_BARCODE_FOLDER, runid, fn)
        print_blastn_1(barcode, TMP_BARCODE_FOLDER, runid, fn)
#        print_blastn_2(barcode, TMP_BARCODE_FOLDER, runid, fn)


# In[84]:

#now qsub all scripts
os.chdir(SCRIPT_FOLDER)
qsub_command = 'qsub %s'
for script in [os.path.join(SCRIPT_FOLDER, x) for x in os.listdir(SCRIPT_FOLDER) if x.endswith('.sh')]:
    qsub_command_stderr = subprocess.check_output(qsub_command % script, shell=True, stderr=subprocess.STDOUT)
    print(qsub_command_stderr)


# In[ ]:
