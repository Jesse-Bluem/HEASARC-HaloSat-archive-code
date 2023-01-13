'''
- run all processes on targets with new data expected
- NOT for regenerating all data since some files error near midnight for searchdb process
- use only to update the uf and cl data base
- input is a comma-separated list of targets numbers
- if no input, run list for last two weeks
- if 'all' input, run all sources with ID<2000 in fits file

Code required:
- HSuf.py, which calls searchdb in /home/data/, and Dohs6_uf, filt_dictionary3, and halosat_analysis7_uf in Server_Code dir
- haloSat_as_flown.py in /home/larocca/asflown/
- HScl.py, which calls clean_from_uf in /home/halosat/hscore/
- HS_exptimes 
- HScl_v1files, which calls halosat_analysis7_uf
- HEASARC_dirlist, which calls checksum_transfer.pl 

Other important files: halosat_checksum.tar, HSufonly.py, and saa_boxes.fits. Do not remove these files. The others are routinely replaced.
Also, the /home/halosat/Server_Code/checksum/ directory has examples from Lorella.
'''
import sys, os
import astropy.io.fits as pyfits

os.system("export DISPLAY=localhost:10.0")
#set directories
main_dir='/home/halosat/AllData/'    #default directory for all data without timecuts
#main_dir='/home/halosat/HEASARC_test/' #test directory for final archive
code_dir='/home/halosat/Server_Code/'  #code directory
sys.path.append(code_dir)

if len(sys.argv)==1:  #if no list given, run list from last two weeks
  #copy over most recent list of fields to run
  os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_recent_targets.txt /home/halosat/Server_Code/.")

  #read list of fields to run
  f = open(code_dir+'haloSat_recent_targets.txt', 'r')
  contents = f.read()
  f.close()
else:  #list is given
  contents=sys.argv[1] #list of fields to run
  if contents=='all':  #if running all fields, get list from targets file
    os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_targets.fits .")  #update target list  
    t = pyfits.open(code_dir+'haloSat_targets.fits')  
    hsIDs = t[1].data['Target_ID']
    t.close()    
    contents=''
    for hs in hsIDs: 
      if hs < 2000: contents=contents+','+str(hs)
    contents=contents[1:]
  
#run processes
#os.system("python /home/halosat/Server_Code/asflown/haloSat_as_flown.py "+main_dir)  #regenerate as flown timeline, only for HEASARC
os.system("python "+code_dir+"HScl.py "+contents+" "+main_dir+" "+main_dir)  #run cl data
os.system("python "+code_dir+"HS_exptimes.py")  #generate fits file with exposure times for planning software
os.system("python "+code_dir+"HScl_v1files.py "+main_dir+" "+contents)  #compress cl files to v1 dir
#os.system("python "+code_dir+"HEASARC_dirlist.py "+main_dir+'v1/') #only for heasarc archival
#os.system("python "+code_dir+"Checksum_check.py "+main_dir) #only for heasarc archival

