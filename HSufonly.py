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

Updates Oct 9, 2022 (JKB)
-   Converted to Python 3.
-   Commented out unused numpy import.

#HEASARC archival times:
- winter 2020 tai times (start, stop): (0, 624499236)  into dir /home/halosat/HEASARC_Feb2020/

#ufonly: to be used for targets where searchdb was successful but the uf processing failed.
   the searchdb files must already exist in the normal directory
"""
import os, sys, datetime, glob, time
#import numpy as np
from numpy import array

# directory containing HaloSat files such as response, target list, Python modules
ancillary_dir='/home/halosat/Server_Code/'
sys.path.append(ancillary_dir)
from Dohs6_uf import Dohs_Single
    
def single_source(hs,main_dir,data_dir,log_name,dohs6_time):
    #initialize reduction timer and target name
    red_time0=time.time()  #start 'timer' for data reduction
    hsID = str('HS{:04d}'.format(hs))  #define string based on field number

    #make sure only one set of searchdb files exist for target
    HSlist = glob.glob(data_dir+hsID+'*.hk')  #check that only one file set remains
    if len(HSlist) > 1: 
        with open(log_name,'a') as log_file: log_file.write('\nERROR: More than one hk file remaining for '+hsID+'! Process failed. Skipping field.')
        print(str('\nERROR: More than one hk file remaining for '+hsID+'! Process failed. Skipping field.'))
        return hs,searchdb_time,dohs6_time #if too many files, returns to function to run next field. Never seen this happen with remove code. NO SEARCHDB_TIME?
       
    #perform uf data reduction
    try:
        if len(HSlist)==1: test=Dohs_Single(hs,main_dir=main_dir,data_dir=data_dir)      #use default uf cut
        else: 
          with open(log_name,'a') as log_file: log_file.write('\nNo data found with searchdb for '+str(hs))
          print('No data found with searchdb for '+str(hs))
          return hs,searchdb_time,dohs6_time       
        if len(test)==3: #check that data for all three DPUs were reduced 
          with open(log_name,'a') as log_file: log_file.write('\nData successfully reduced with unfiltered cuts for '+str(hs))
          print('Data successfully reduced with unfiltered cuts for '+str(hs))
          
          #remove searchdb output since processing successful       
          remove_list = glob.glob(data_dir+hsID+'*')
          for r in remove_list: os.remove(r)          
          
        else:  #break logic if data for all DPUs were not reduced      
          raise Exception('failed')
    except:# Exception as ex:
        with open(log_name,'a') as log_file: log_file.write('\nData reduction failed for '+str(hs))  
        print('Data reduction failed for '+str(hs))
        return hs,dohs6_time
    stime = time.time()-red_time0        
    print(str('Data reduction time: {:.3f}s'.format(stime)))
    with open(log_name,'a') as log_file: log_file.write('\nData reduction time: {:.3f}s'.format(stime))
    dohs6_time += stime
    
    return -1,dohs6_time  #return statement if all processes are successful
    
'''
---------------------- Main Program -----------------------------------
'''    

#update haloSat_targets.fits files for Dohs6_uf processes
os.chdir(ancillary_dir)
os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_targets.fits .")
    
#get inputs either from input argument or file
if len(sys.argv) == 3:
  contents=sys.argv[1] #list of fields to run
  main_dir=sys.argv[2]  #dir to put output
  if not os.path.exists(main_dir): os.makedirs(main_dir)
else:
  print('You must specify the target list and the directory.')
  raise Exception('failed')
  
#collect variables into standard format  
hs_list = contents.split(',') # list of sources to be run
data_dir=main_dir+'searchdb/'
if not os.path.exists(data_dir): os.makedirs(data_dir)
#print hs_list, run_par

#set log file
os.chdir(data_dir)
file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')  #for file name
log_name = data_dir+'HSPROCuf_'+file_time+'.txt'  #change to txt for easier opening

#output file dirs to log file
dohs6_time = 0. #initialize time sums
with open(log_name,'w') as log_file: 
  log_file.write('\nUnfiltered output will be stored in {} if successful.'.format(main_dir))

#run fields
hs_failed=[]
hs_list=array(hs_list,dtype=int)
with open(log_name,'a') as log_file: log_file.write('\nRunning for sources '+contents)
for i in range(len(hs_list)): #run list of fields
  print('Processing '+str(i+1)+' of '+str(len(hs_list))+' fields.')
  test,dohs6_time=single_source(hs_list[i],main_dir,data_dir,log_name,dohs6_time)  #if successful, length will be zero
  if (test>0): hs_failed.append(test)  #processing failed for source number

#output time stats
with open(log_name,'a') as log_file: 
  log_file.write('\nTotal data reduction time: {:.3f} s'.format(dohs6_time))
    
if len(hs_failed)>0: 
  failed=str(hs_failed[0])
  for f in range(len(hs_failed)-1): failed+=','+str(hs_failed[f+1])
  with open(log_name,'a') as log_file:log_file.write('\n\nProcessing failed for halosat ID numbers {}.'.format(failed))
else:
  with open(log_name,'a') as log_file:log_file.write('\n\nAll processing completed for all halosat ID numbers given.')
