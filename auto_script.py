'''
- run all processes on targets with new data expected (JKB: Irrelevent, mision completed)
- NOT for regenerating all data since some files error near midnight for searchdb process (JKB: Midnight bug fixed, this generates all the data now.)
- use only to update the uf and cl data base
- input is a comma-separated list of targets numbers, JKB: followed by a second command line argument of 'full' or 'run1' or 'run2'.
- JKB: run1 and run2 are subsets of the pipeline code that can be run separately to break up the process. run1 first does the searchdb call and other processing code without HEASARC archive stuff. run2 skips searchdb and does the full HEASARC archive.
- if no input, run list for last two weeks, JKB: with 'full', running all called pipeline code pieces
- if 'all' input, run all sources with ID<2000 in fits file, JKB: with 'full', running all called pipeline code pieces

Code required:
- HSuf.py, which calls searchdb in /home/data/, and Dohs6_uf, filt_dictionary3, and halosat_analysis7_uf in Server_Code dir
- haloSat_as_flown.py in /home/larocca/asflown/
- HScl.py, which calls clean_from_uf in /home/halosat/hscore/
- HS_exptimes 
- HScl_v1files, which calls halosat_analysis7_uf
- HEASARC_dirlist, which calls checksum_transfer.pl 

Other important files: halosat_checksum.tar, HSufonly.py, and saa_boxes.fits. Do not remove these files. The others are routinely replaced.
Also, the /home/halosat/Server_Code/checksum/ directory has examples from Lorella.

JKB: Example command line syntax: python /home/halosat/Server_Code/auto_script.py 'all' 'run1'

'''
import sys, os
import astropy.io.fits as pyfits

#set directories
main_dir='/home/halosat/AllData/'    #default directory for all data without timecuts
code_dir='/home/halosat/Server_Code/'  #code directory
sys.path.append(code_dir)

runset='full' #assume run all steps if not specified

if len(sys.argv)==1:  #if no list given, run list from last two weeks
  #copy over most recent list of fields to run
  os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_recent_targets.txt /home/halosat/Server_Code/.")

  #read list of fields to run
  f = open(code_dir+'haloSat_recent_targets.txt', 'r')
  contents = f.read()
  f.close()
  
else:  #list is given
  contents=sys.argv[1] #list of fields to run
  runset=sys.argv[2] #JKB. Add second variable to command line for which pieces of code to run.
  if contents=='all':  #if running all fields, get list from targets file
    os.system("sudo scp -i /home/halosat/.ssh/id_rsa -P53222 halosat@halosat.physics.uiowa.edu:/home/halosat/haloSat/haloSat_targets.fits .")  #update target list  
    t = pyfits.open(code_dir+'haloSat_targets.fits')  
    hsIDs = t[1].data['Target_ID']
    #print(hsIDs)
    t.close()    
    contents=''
    for hs in hsIDs: 
      if hs < 2000: contents=contents+','+str(hs)
    contents=contents[1:]
 
if runset == 'full':  
  #run processes - JKB 12/7: this is original full sequence. Split into two for easier tracking and commenting out. Or for personal preference, really.
  print('Running all file processing steps.')
  os.system("python "+code_dir+"HSuf.py "+contents)  #no timecuts to run all data, default directories
  os.system("python /home/halosat/Server_Code/asflown/haloSat_as_flown.py "+main_dir)  #regenerate as flown timeline, only for HEASARC
  os.system("python "+code_dir+"HScl.py "+contents+" "+main_dir+" "+main_dir)  #run cl data
  os.system("python "+code_dir+"HS_exptimes.py")  #generate fits file with exposure times for planning software
  os.system("python "+code_dir+"HScl_v1files.py "+main_dir+" "+contents)  #compress cl files to v1 dir
  os.system("python "+code_dir+"HEASARC_dirlist.py "+main_dir+'v1/') #only for heasarc archival
  os.system("python "+code_dir+"Checksum_check.py "+main_dir) #only for heasarc archival

elif runset == 'run1':
  #run1 - JKB 12/7: This is what is run in the first pass for archival
  print('Running run1 reduced file processing steps.')  
  os.system("python "+code_dir+"HSuf.py "+contents)  #no timecuts to run all data, default directories
  os.system("python "+code_dir+"HScl.py "+contents+" "+main_dir+" "+main_dir)  #run cl data
  os.system("python "+code_dir+"HS_exptimes.py")  #generate fits file with exposure times for planning software
  os.system("python "+code_dir+"HScl_v1files.py "+main_dir+" "+contents)  #compress cl files to v1 dir

elif runset == 'run2':
  #run2 - JKB 12/7:  This is the second pass that finishes the archival process.
  print('Running run2 reduced final file processing steps.')
  os.system("python /home/halosat/Server_Code/asflown/haloSat_as_flown.py "+main_dir)  #regenerate as flown timeline, only for HEASARC
  os.system("python "+code_dir+"HScl.py "+contents+" "+main_dir+" "+main_dir)  #run cl data
  os.system("python "+code_dir+"HS_exptimes.py")  #generate fits file with exposure times for planning software
  os.system("python "+code_dir+"HScl_v1files.py "+main_dir+" "+contents)  #compress cl files to v1 dir
  os.system("python "+code_dir+"HEASARC_dirlist.py "+main_dir+'v1/') #only for heasarc archival
  os.system("python "+code_dir+"Checksum_check.py "+main_dir) #only for heasarc archival
  
else:
  print("Second command line variable does not match any known setting: Options are 'full', 'run1', or 'run2'.")
