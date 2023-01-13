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

HScl_v1files.py Feb 6 2020
- runs ftchecksum, ftverify, and compresses all clean files
- compressed files are put into the v1 directory for HEASARC
- two types of calls:
  python /home/halosat/Server_Code'/HScl_v1files '/home/halosat/' '369'
      : first input if the top directory where v1_local and v1 directories exist
      : second input is the target number (or target numbers separated by commas)
  python /home/halosat/Server_Code'/HScl_v1files '/home/halosat/'
      : first input if the top directory where v1_local and v1 directories exist
      : lack of the second input runs all targets

Updates 10/9/2022 JKB
- Updated to python 3.
- Removed unused numpy import.

JKB NOTE: THIS CODE'S single_source DOESN'T MATCH HScl's FUNCTION
"""
import os, sys, datetime, glob
#import numpy as np
import halosat_analysis7_uf as hs7

# directory containing HaloSat files such as response, target list, Python modules
ancillary_dir='/home/halosat/Server_Code/'
#ancillary_dir='E:/HaloSat/halosat_sourceanalysis/Dohs_Code/'
sys.path.append(ancillary_dir)
   
def single_source(hs,main_dir,log_name):
    #set directory structure
    time0 = time.time()
    hsname=str('{:04d}01'.format(hs))
    cl_dir = main_dir+'v1_local/'+hsname+'/products/'  #starting dir
    gz_dir = main_dir+'v1/'+hsname+'/products/'   #final dir

    #check that cl dir and cl files exist for this target, return if no
    if not os.path.exists(cl_dir): 
      with open(log_name,'a') as log_file:log_file.write('\nNo files exist for target {}.'.format(hs))
      print('No files exist for target '+str(hs)+'.')
      return hs, 0
    pfnames = glob.glob(cl_dir+'hs'+hsname+'.*')
    if len(pfnames)==0: 
      with open(log_name,'a') as log_file:log_file.write('\nNo files exist for target {}.'.format(hs))
      print('No files exist for target '+str(hs)+'.')
      return hs, 0  #return if no files found

    #set list of file roots for this source
    if not os.path.exists(gz_dir): os.makedirs(gz_dir)   #make dir for compressed files
    ext_list, file_list = ['_s14.pi','_s38.pi','_s54.pi','_s14_cl.evt','_s38_cl.evt','_s54_cl.evt'], []  #only fits files
    for i in range(len(ext_list)):   #define list of files
        file_list.append('hs'+hsname+ext_list[i])
    #print cl_list, gz_list

    #check files with fverify and fchecksum
    with open(gz_dir+'hs'+hsname+'_filelist.txt','w') as report_file:  #write file list to a text file
      for f in file_list: report_file.write('{}{}\n'.format(cl_dir,f))
    os.system('ftchecksum @'+gz_dir+'hs'+hsname+'_filelist.txt update=yes datasum=yes chatter=1')  #checksum all files at once
    os.system('ftverify @'+gz_dir+'hs'+hsname+'_filelist.txt outfile='+gz_dir+'hs'+hsname+'_fverify.txt clobber=yes')  #fverify all files at once
    os.remove(gz_dir+'hs'+hsname+'_filelist.txt')  #remove file list
    
    #copy and gzip fits files to heasarc dir, add fchecksum output to report file
    with open(cl_dir+'hs'+hsname+'_Fresults.txt','w') as report_file:
      report_file.write('Ftchecksum Results:\n--------------------')        
    for f in file_list:
      print(cl_dir+f)
      hs7.print_checksum(cl_dir+f,cl_dir+'hs'+hsname+'_Fresults.txt')  #should error if checksum failed, prints checksum keywords to report file
      if os.path.exists(gz_dir+f+".gz"): os.remove(gz_dir+f+".gz")
      os.system("imcopy "+cl_dir+f+" "+gz_dir+f+".gz")          # compress and move fits files to output dir

    #save compressed version of pdf file in v1 dir  
    os.system('cp '+cl_dir+'hs'+hsname+'.pdf '+gz_dir+'hs'+hsname+'.log.pdf')  #copy txt file to new directory, will fail if DNE
    if os.path.exists(gz_dir+'hs'+hsname+'.log.pdf.gz'):  #remove txt.gz file if it already exists to avoid asking permission
            os.remove(gz_dir+'hs'+hsname+'.log.pdf.gz')     
    os.system('gzip '+gz_dir+'hs'+hsname+'.log.pdf')  #creates gz file and removes original. 7zip opens file just fine.
    
    #add fverify output to file
    f = open(gz_dir+'hs'+hsname+'_fverify.txt','r')
    contents=f.read()
    f.close
    os.remove(gz_dir+'hs'+hsname+'_fverify.txt') #remove file
    with open(cl_dir+'hs'+hsname+'_Fresults.txt','a') as report_file:  #add fverify output to report file
      report_file.write('\n\nFtverify Results:\n--------------------\n{}'.format(contents))    
      
    #add comment to log file
    stime = time.time()-time0
    with open(log_name,'a') as log_file:log_file.write('\nFiles successfully checked and compressed for target {} in {}s.'.format(hs,stime))
    print(str('\nFiles successfully checked and compressed for target {} in {}s.'.format(hs,stime)))
    
    return 0, stime  #return time
    
def all_sources(main_dir,log_name):  #code to run process on all HSIDs 372 and under and dark earth fields
    with open(log_name,'a') as log_file:log_file.write('\nRunning for all sources.')  
    hs_failed, gz_time = [], 0.
    print('Running for all sources')
    for hs in range(371): #run science fields
        print('Processing '+str(hs+1)+' of 418 fields.')
        if hs==0: contents=str(hs+1)  #build list of sources
        else: contents+=','+str(hs+1)
        status_check, single_time = single_source(hs+1,main_dir,log_name)  #run source
        gz_time += single_time   #total compression time
        if status_check>0: hs_failed.append(status_check)  #build failed list
                   
    for hs in range(45): #run dark earth fields, same as above
        print('Processing '+str(hs+372)+' of 418 fields.')
        contents+=','+str(hs+1001)
        status_check, single_time = single_source(hs+1001,main_dir,log_name)
        gz_time += single_time
        if status_check>0: hs_failed.append(status_check)
        
    return contents,hs_failed,gz_time

'''
---------------------- Main Program -----------------------------------
'''    
from numpy import array
import time

#get input
contents_check, contents = -1, -1  #initialize check for running with automatic list 
if len(sys.argv) == 3:
  contents=sys.argv[2] #list of fields to run
  hs_list = contents.split(',') # list of sources to be run
  main_dir=sys.argv[1]  #dir above v1 and v1_local
elif (len(sys.argv) == 2):  #only dir given, assume default time interval and running all targets
  contents_check=1
  main_dir=sys.argv[1]  #dir above v1 and v1_local

#set log file
os.chdir(main_dir)
file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H%M%S')  #for file name
log_name = main_dir+'v1_local/HSPROCclgz_'+file_time+'.txt'  #log file name
with open(log_name,'w') as log_file: #output file dir to log file
  log_file.write('Compressed output will be stored in {}.'.format(main_dir+'v1/'))

#run fields
if contents_check==1: 
  contents,hs_failed,gz_time = all_sources(main_dir,log_name) #returns list of failed and list attempted
else:
  hs_failed, gz_time = [], 0.
  hs_list=array(hs_list,dtype=int)
  with open(log_name,'a') as log_file: log_file.write('\nRunning for sources '+contents)
  for i in range(len(hs_list)): #run list of fields
    print('Checking and compressing '+str(i+1)+' of '+str(len(hs_list))+' fields.')
    status_check, single_time = single_source(hs_list[i],main_dir,log_name)  #if successful, length will be zero
    gz_time += single_time
    if status_check>0: hs_failed.append(status_check)    

#output time stats
with open(log_name,'a') as log_file: 
  log_file.write('\n\nTotal checking and compression time: {:.3f} s'.format(gz_time))
print(str('\nTotal checking and compression time: {:.3f} s'.format(gz_time)))

#output list of failed targets
if len(hs_failed)>0: 
  failed=str(hs_failed[0])
  for f in range(len(hs_failed)-1): failed+=','+str(hs_failed[f+1])
  with open(log_name,'a') as log_file:log_file.write('\nChecking and compression failed for halosat ID numbers {}.'.format(failed))
else:
  with open(log_name,'a') as log_file:log_file.write('\nAll checking and compression completed for all halosat ID numbers given.')
