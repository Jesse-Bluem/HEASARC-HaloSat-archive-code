# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:56:12 2020

@author: rringuette
Create directory listing for HEASARC data archival.
To be run on halosat2

"""
import os, datetime, sys
#import numpy as np

#read directory from command line
hs_dir = sys.argv[1]

#open directory list file with today's date in the filename
date = datetime.datetime.today().strftime('%Y%m%d')
dir_tree = hs_dir+'dirlist_'+date+'.txt'

#print directories containing files to output file
check=0
for path, dirs, files in os.walk(hs_dir):
    dirs.sort()  #sort directories in ascii order
    if len(files)>0 and path!=hs_dir: 
      if check==0: 
        with open(dir_tree,'w') as log_file: log_file.write('{}'.format(path))
        check=1
      else: 
        with open(dir_tree,'a') as log_file: log_file.write('\n{}'.format(path))
        
#run checksum script
os.system("perl /home/halosat/Server_Code/checksum_transfer.pl "+dir_tree+" O "+hs_dir+"checksum_origin_"+date+".txt")