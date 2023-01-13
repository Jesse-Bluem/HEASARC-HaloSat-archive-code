#Update 10/9/2022 JKB - changed to Python 3.

import numpy as np
import glob, datetime, sys

def string_check(string,sub):
    if (string.find(sub) != -1): return 1
    else: return 0

#read directory from command line
hs_dir = sys.argv[1]

#set output file
file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H%M%S')  #for file name
log_name = hs_dir+'/searchdb/Checksum_'+file_time+'.txt'  #change to txt for easier opening
with open(log_name,'w') as log_file: 
  log_file.write('Checking the checksum file for correct file list. Only missing or incorrect files notices will be printed.')

#get list of files from file
cfile = glob.glob(hs_dir+'/v1/checksum_origin_*')
#print cfile[0]
check_file = open(cfile[0], 'r')
fdata = check_file.read()
data = fdata.split('\n')
c = []
for i in range(len(data)): c.append(data[i].split('  ')[-1])
c_arr = np.array(c, dtype=str)

#check that each target has correct files (regular targets)
for i in range(375):
    #set correct file names
    t = str('hs{:04d}01'.format(i+1))
    file_set = ['unfiltered/'+t+'.att.gz','unfiltered/'+t+'_s14.hk.gz','unfiltered/'+t+'_s38.hk.gz',\
                'unfiltered/'+t+'_s54.hk.gz','unfiltered/'+t+'_s14_uf.evt.gz','unfiltered/'+t+'_s38_uf.evt.gz',\
                'unfiltered/'+t+'_s54_uf.evt.gz','products/'+t+'.log.pdf.gz','products/'+t+'_s14.pi.gz',\
                'products/'+t+'_s38.pi.gz','products/'+t+'_s54.pi.gz','products/'+t+'_s14_cl.evt.gz',\
                'products/'+t+'_s38_cl.evt.gz','products/'+t+'_s54_cl.evt.gz']
    
    #get indices of all files for target
    idx = []
    for j in range(len(c_arr)):
        c = string_check(c_arr[j],t)
        #print c, j
        if (c==1): idx.append(j)
    t_files = c_arr[idx]
    #print t_files
    
    #check files for target for correct list, print files not found
    for fname in file_set:
        count=0
        for j in range(len(t_files)):
            c = string_check(t_files[j],fname)
            count+=c
        if (count==0): 
            with open(log_name,'a') as log_file: log_file.write('\nError! File not found: {}'.format(fname))
                
    #check for _mod files
    for j in range(len(t_files)):
        c = string_check(t_files[j],'mod')
        if (c==1):
            with open(log_name,'a') as log_file: log_file.write('\nError! Incorrect file found: {}'.format(t_files[j]))

#dark earth targets
for i in range(48):
    #set correct file names
    t = str('hs{:04d}01'.format(i+1001))
    file_set = ['unfiltered/'+t+'.att.gz','unfiltered/'+t+'_s14.hk.gz','unfiltered/'+t+'_s38.hk.gz',\
                'unfiltered/'+t+'_s54.hk.gz','unfiltered/'+t+'_s14_uf.evt.gz','unfiltered/'+t+'_s38_uf.evt.gz',\
                'unfiltered/'+t+'_s54_uf.evt.gz','products/'+t+'.log.pdf.gz','products/'+t+'_s14.pi.gz',\
                'products/'+t+'_s38.pi.gz','products/'+t+'_s54.pi.gz','products/'+t+'_s14_cl.evt.gz',\
                'products/'+t+'_s38_cl.evt.gz','products/'+t+'_s54_cl.evt.gz']
    
    #get indices of all files for target
    idx = []
    for j in range(len(c_arr)):
        c = string_check(c_arr[j],t)
        #print c, j
        if (c==1): idx.append(j)
    t_files = c_arr[idx]
    #print t_files
    
    #check files for target for correct list, print files not found
    for fname in file_set:
        count=0
        for j in range(len(t_files)):
            c = string_check(t_files[j],fname)
            count+=c
        if (count==0): 
            with open(log_name,'a') as log_file: log_file.write('\nError! File not found: {}'.format(fname))
                
    #check for _mod files
    for j in range(len(t_files)):
        c = string_check(t_files[j],'mod')
        if (c==1):
            with open(log_name,'a') as log_file: log_file.write('\nError! Incorrect file found: {}'.format(t_files[j]))
            
print('All targets checked.')
