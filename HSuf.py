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

Updates Oct 09, 2022
-	Updated to python 3 by JKB

#HEASARC archival times:
- winter 2020 tai times (start, stop): (0, 624499236)  into dir /home/halosat/HEASARC_Feb2020/
"""
import os, sys, datetime, glob, time
import numpy as np
from numpy import array

# directory containing HaloSat files such as response, target list, Python modules
ancillary_dir='/home/halosat/Server_Code/'
#ancillary_dir='E:/HaloSat/halosat_sourceanalysis/Dohs_Code/'
sys.path.append(ancillary_dir)
from Dohs6_uf import Dohs_Single
    
def single_source(hs,main_dir,data_dir,log_name,run_par,searchdb_time,dohs6_time):
    #pull variables from run_par
    gti_min=int(run_par[0])
    offset=float(run_par[1])
    start_time=run_par[2][0]
    end_time=run_par[2][1]
    run_str = ' with gti_min='+str(gti_min)+', offset<='+str(offset)+', start_time='+str(start_time)+', and end_time='+str(end_time)+'.'
    run_par_str = ' '+str(gti_min)+' '+str(offset)+' '+str(start_time)+' '+str(end_time)
    hsID = str('HS{:04d}'.format(hs))  #define string based on field number         

    #run searchdb for field with standard values
    print('\nRunning searchbd for HaloSat field ', str(hs), run_str)#, log_name
    time0 = time.time()
    os.system('python /home/halosat/Server_Code/searchdb.py '+str(hs)+run_par_str+' '+log_name) #JKB changed from /data to /server-cdoe OCT 8 2022
    stime = time.time()-time0
    searchdb_time+=stime
    print(stime)
    #check for success or failure in output, exit if failure
    with open(log_name,'r') as log_file: check=log_file.readlines()[-1]
    print('----------------------------')
    print(check, len(check.split('failed')))
    #stop
    if len(check.split('failed'))>1: return hs,searchdb_time,dohs6_time
    
    #print searchdb time, initialize reduction timer
    print(str('Data collection time: {:.3f}s'.format(stime))) 
    with open(log_name,'a') as log_file: log_file.write('\nData collection time: {:.3f}s'.format(stime))
    red_time0=time.time()  #start 'timer' for data reduction
    
    #get a list of file roots for this source, remove old ones
    HSlist = np.sort(glob.glob(data_dir+hsID+'*.hk'))
    for i in range(len(HSlist)-1): #save more recent file set, remove all earlier ones
        filename = HSlist[i].split('.')[0]
        remove_list = glob.glob(data_dir+filename+'*')
        for r in remove_list: os.remove(r)
    HSlist = glob.glob(data_dir+hsID+'*.hk')  #check that only one file set remains
    if len(HSlist) > 1: 
        with open(log_name,'a') as log_file: log_file.write('\nERROR: More than one hk file remaining for '+hsID+'! Process failed. Skipping field.')
        print(str('\nERROR: More than one hk file remaining for '+hsID+'! Process failed. Skipping field.'))
        return hs,searchdb_time,dohs6_time #if too many files, returns to function to run next field. Never seen this happen with remove code.
    
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
    except Exception as ex:
        print(ex) #JKB oct 2022 - lets give the user some sort of clue about what happened here, jeeze...
        with open(log_name,'a') as log_file: log_file.write('\nData reduction failed for '+str(hs))  
        print('Data reduction failed for '+str(hs)) 
        return hs,searchdb_time,dohs6_time
    stime = time.time()-red_time0        
    print(str('Data reduction time: {:.3f}s'.format(stime)))
    with open(log_name,'a') as log_file: log_file.write('\nData reduction time: {:.3f}s'.format(stime))
    dohs6_time += stime
    
    return -1,searchdb_time,dohs6_time  #return statement if all processes are successful
    
def all_sources(main_dir,data_dir,log_name,run_par,searchdb_time,dohs6_time):  #code to run process on all HSIDs 372 and under and dark earth fields
    hs_failed = []
    with open(log_name,'a') as log_file:log_file.write('\nRunning for all sources.')  
    print('Running for all sources')
    for hs in range(375): #run science fields
        print('Processing '+str(hs+1)+' of 419 fields.')
        if hs==0: contents=str(hs+1)
        else: contents+=','+str(hs+1)
        test,searchdb_time,dohs6_time=single_source(hs+1,main_dir,data_dir,log_name,run_par,searchdb_time,dohs6_time)
        if test>0: hs_failed.append(test)  #if successful, length will be zero
        else: #remove searchdb output since processing successful
          hsID = str('HS{:04d}'.format(hs+1))  #define string based on field number        
                   
    for hs in range(46): #run dark earth fields, same as above
        print('Processing '+str(hs+375)+' of 419 fields.')
        contents+=','+str(hs+1001)
        test,searchdb_time,dohs6_time=single_source(hs+1001,main_dir,data_dir,log_name,run_par,searchdb_time,dohs6_time)
        if test>0: hs_failed.append(test)
        else: #remove searchdb output since processing successful
          hsID = str('HS{:04d}'.format(hs+1))    #define string based on field number    UNUSED  
                   
    return hs_failed, contents,searchdb_time,dohs6_time

'''
---------------------- Main Program -----------------------------------
'''    

#update haloSat_targets.fits files for Dohs6_uf and searchdb processes
searchdb_dir='/home/halosat/Server_Code/' #server - JKB changed from data to halosat server_code oct 8 2022
os.chdir(searchdb_dir)
os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_targets.fits .")
os.chdir(ancillary_dir)
os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_targets.fits .")
    
#get inputs either from input argument or file
contents_check=0  #initialize check for running with automatic list 
if len(sys.argv) == 4:
  contents=sys.argv[1] #list of fields to run
  timecuts = array(sys.argv[2].split('[')[1].split(']')[0].split(','),dtype=int)  #timecuts
  main_dir=sys.argv[3]  #dir to put output
  if not os.path.exists(main_dir): os.makedirs(main_dir)
elif (len(sys.argv) == 2): #using default timecuts, all output should be directed to standard directories
  contents=sys.argv[1] #list of fields to run
  timecuts=array([0,10E15],dtype=int)   #default time cuts
  main_dir='/home/halosat/AllData/'       #server
elif (len(sys.argv) == 1):  #no arguments given, assume default time interval and get target list to run
  #copy over most recent list of fields to run
  os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_recent_targets.txt .")

  #read list of fields to run
  f = open(ancillary_dir+'haloSat_recent_targets.txt', 'r')
  contents = f.read()
  f.close()
  contents_check=1
  timecuts=array([0,10E15],dtype=int)
else:
  print('You must specify the directory structure when using timecuts.')
  raise Exception('failed')
  
#collect variables into standard format  
run_par=[64,0.25,timecuts]  #gti_min, offset, timecuts
hs_list = contents.split(',') # list of sources to be run
data_dir=main_dir+'searchdb/'
if not os.path.exists(data_dir): os.makedirs(data_dir)
#print hs_list, run_par

#set log file
os.chdir(data_dir)
file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H%M%S')  #for file name
log_name = data_dir+'HSPROC_'+file_time+'.txt'  #change to txt for easier opening
print(log_name)
#stop

#output file dirs to log file
searchdb_time, dohs6_time = 0., 0. #initialize time sums
with open(log_name,'w') as log_file: 
  log_file.write('Searchdb output will be stored in {} if errors occur.'.format(data_dir))
  log_file.write('\nUnfiltered output will be stored in {} if successful.'.format(main_dir))

#run fields
if hs_list[0]=='all': 
  hs_failed, contents,searchdb_time,dohs6_time = all_sources(main_dir,data_dir,log_name,run_par,searchdb_time,dohs6_time) #returns list of failed and list attempted
else:
  hs_failed=[]
  hs_list=array(hs_list,dtype=int)
  with open(log_name,'a') as log_file: log_file.write('\nRunning for sources '+contents)
  for i in range(len(hs_list)): #run list of fields
    print('Processing '+str(i+1)+' of '+str(len(hs_list))+' fields.')
    test,searchdb_time,dohs6_time=single_source(hs_list[i],main_dir,data_dir,log_name,run_par,searchdb_time,dohs6_time)  #if successful, length will be zero
    if (test>0): hs_failed.append(test)  #processing failed for source number

#output time stats
with open(log_name,'a') as log_file: 
  log_file.write('\n\nTotal data collection time: {:.3f} s'.format(searchdb_time))
  log_file.write('\nTotal data reduction time: {:.3f} s'.format(dohs6_time))
  log_file.write('\nTotal data processing time: {:.3f} s'.format(searchdb_time+dohs6_time))
    
if len(hs_failed)>0: 
  failed=str(hs_failed[0])
  for f in range(len(hs_failed)-1): failed+=','+str(hs_failed[f+1])
  with open(log_name,'a') as log_file:log_file.write('\n\nProcessing failed for halosat ID numbers {}.'.format(failed))
else:
  with open(log_name,'a') as log_file:log_file.write('\n\nAll processing completed for all halosat ID numbers given.')
