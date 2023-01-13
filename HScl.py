#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:58:18 2019
@author: ringuette

master controller for uf and searchdb software

Updates Dec 18, 2019
-	Update code naming dependencies
-	Streamline code into functions 
-	Add command line inputs
-	Add error handling and reporting
-	Update directory structures and file names
-	Remove searchdb output after successful completion of Dohs6 per target

HSall.py Feb 5 2020 (and on)
- add lines to run cl software, but no gzip, checksum, etc functionality
- Note: UNIX times are 946684800-37 seconds larger than TAI times (tai_time=unix_time-epoch_diff+leap seconds)

Updates Oct 9, 2022 - JKB
- converted to Py3.
- Commented out unsed variables.
- removed 7 halo fields from special cuts on 12/7. Special cuts here were a mistake based on early contaminated data.
"""
import os, sys, datetime, time
import numpy as np

# directory containing HaloSat files such as response, target list, Python modules
ancillary_dir='/home/halosat/Server_Code/'
sys.path.append(ancillary_dir)
    
def single_source(hs,main_dir,data_dir,log_name,cl_time):
    #hsID = str('HS{:04d}'.format(hs))  #define string based on field number         UNUSED
    
    #setup cl data reduction
    time0 = time.time()
    hsname=str('{:04d}01'.format(hs))
    type1_list = [45,3,46,54,56,348] #non-default cuts for these sources JKB 12/7 - removed 293,279,84,275,112,102,233, - special cut is not required
    if hs==2: hard, vle = '5.00', '2.00' # 4.00 and 1.50 prior to 2020-09-17 JRR
    elif hs in type1_list: hard, vle = '0.50', '0.75'
    elif (hs==43) or (hs==44): hard, vle = '1.00', '0.75'
    else: hard, vle = '0.16', '0.75'  
    #print(hs, hard, vle)
    uf_dir = data_dir+'v1_local/'+hsname+'/unfiltered/'
    cl_dir = main_dir+'v1_local/'+hsname+'/products/'       
    if not os.path.exists(uf_dir): 
      print(str('No unfiltered data for target {}. Skipping target.'.format(hs)))
      return hs, cl_time  #don't bother if the uf data DNE 
    uf_dir+='hs'+hsname   
    if not os.path.exists(cl_dir): os.makedirs(cl_dir)
    cl_dir+='hs'+hsname
    if hs<1000: cl_string = '''"'''+uf_dir+'''" "'''+cl_dir+'''" "(hk2['LC_HARD'] <= '''+hard+\
        ''') & (hk2['LC_VLE'] <= '''+vle+''')"'''
    elif hs>1000: cl_string = '''"'''+uf_dir+'''" "'''+cl_dir+'''" "(hk2['LC_ALSI']/hk2['LC_UP'] > 1.50)"'''
        
    #perform cl data reduction
    os.system('python /home/halosat/Server_Code/clean_from_uf.py '+cl_string)    #WHY hscore? switched to py3 server code. JKB oct 12 2022
    stime = time.time()-time0
    print(str('Data cl reduction time: {:.3f}s\n'.format(stime)))
    with open(log_name,'a') as log_file: 
        log_file.write('\nData successfully filtered with hard < '+\
            hard+' ct/s and vle < '+vle+' ct/s cuts for '+str(hs))
        log_file.write('\nData filtering time: {:.3f}s\n'.format(stime))
    cl_time += stime
    
    return -1,cl_time  #return statement if all processes are successful

def all_sources(main_dir,data_dir,log_name,cl_time):  #code to run process on all HSIDs 372 and under and dark earth fields
    hs_failed = []
    with open(log_name,'a') as log_file:log_file.write('\nRunning for all sources.')  
    print('Running for all sources')
    for hs in range(375): #run science fields
        print('Processing '+str(hs+1)+' of 419 fields.')
        if hs==0: contents=str(hs+1)
        else: contents+=','+str(hs+1)
        test,cl_time=single_source(hs+1,main_dir,data_dir,log_name,cl_time)
        if test>0: hs_failed.append(test)  #if successful, length will be zero
        else: #remove searchdb output since processing successful
          hsID = str('HS{:04d}'.format(hs+1))  #define string based on field number        
                   
    for hs in range(46): #run dark earth fields, same as above
        print('Processing '+str(hs+375)+' of 419 fields.')
        contents+=','+str(hs+1001)
        test,cl_time=single_source(hs+1001,main_dir,data_dir,log_name,cl_time)
        if test>0: hs_failed.append(test)
        else: #remove searchdb output since processing successful
          hsID = str('HS{:04d}'.format(hs+1))    #define string based on field number      UNUSED
                   
    return hs_failed, contents,cl_time

'''
---------------------- Main Program -----------------------------------
'''    
    
#get inputs either from input argument or file
if len(sys.argv) == 4:
  contents=sys.argv[1] #list of fields to run
  main_dir=sys.argv[2]  #dir to put output
  data_dir=sys.argv[3]  #dir over searchdb, v1, and v1_local dirs with chosen uf data
  if not os.path.exists(main_dir): os.makedirs(main_dir)
else:
  print('You must specify the target list, the output and data directories.')
  raise Exception('failed')
#collect variables into standard format  
hs_list = contents.split(',') # list of sources to be run

#set log file
os.chdir(data_dir+'/searchdb/')
file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H%M%S')  #for file name
log_name = data_dir+'/searchdb/HSPROCcl_'+file_time+'.txt'  #change to txt for easier opening

#output file dirs to log file
cl_time = 0. #initialize time sums
with open(log_name,'w') as log_file: 
  log_file.write('Cleaned output will be stored in {} if successful.'.format(main_dir))

#run fields
hs_failed=[]
if hs_list[0]=='all': 
  hs_failed, contents,cl_time = all_sources(main_dir,data_dir,log_name,cl_time) #returns list of failed and list attempted
else:
  hs_list=np.array(hs_list,dtype=int)
  with open(log_name,'a') as log_file: log_file.write('\nRunning for sources '+contents)
  for i in range(len(hs_list)): #run list of fields
    print('Processing '+str(i+1)+' of '+str(len(hs_list))+' fields.')
    test,cl_time=single_source(hs_list[i],main_dir,data_dir,log_name,cl_time)  #if successful, length will be zero
    if (test>0): hs_failed.append(test)  #processing failed for source number

#output time stats
with open(log_name,'a') as log_file: 
  log_file.write('\nTotal data cl reduction time: {:.3f} s'.format(cl_time))
   
if len(hs_failed)>0: 
  failed=str(hs_failed[0])
  for f in range(len(hs_failed)-1): failed+=','+str(hs_failed[f+1])
  with open(log_name,'a') as log_file:log_file.write('\n\nClean processing failed for halosat ID numbers {}.'.format(failed))
else:
  with open(log_name,'a') as log_file:log_file.write('\n\nAll clean processing completed for all halosat ID numbers given.')
