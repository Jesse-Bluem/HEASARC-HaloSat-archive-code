'''
Dohs6.py: Sep 6, 2019 Rebecca Ringuette
- Modified from previous versions from Phil Kaaret
- Standardized to produce HEASARC diretory structure and files,
    printed and plot outputs, and correct file headers

#main function: process files for list of source numbers
def Dohs_List(source_list,main_dir='/home/halosat/',data_dir='/home/halosat/v2p1/',\
                  ancillary_dir='/home/halosat/Server_Code/',searchdb_dir='/home/halosat/v2p1/',\
                  cut_note='uf',gti_min=64,timecuts=[],\
                  gti_max=2000,map_type=[0],lc_list=[0,3,4,6],loud=0,\
                  creator='Dohs6',filter_pars=[]):  
Dohs_List([source_num])   #use default uf cut
Dohs_List([hs],cut_note='cl')   #change to use clean cut specific to source # from library
Dohs_List([hs],cut_note='cl',filter_pars=[['LC_HARD','LC_VLE'],[[-1,0.25],[-1,2.]],[64,64]])
This cuts on LC_HARD band using 64s time bins to be greater than -1 and less than 0.25 ct/s
Also cuts on LC_VLE band using 64s time bins to be greater than -1 and less than 2.0 ct/s
Can add additional cuts on any time scale that is a multiple of 8.
e.g. filter_pars=[['LC_ALL'],[[20.,50.]],[8]] cuts on the total count rate using 
8 second bins to be greater than 20 and less than 50 counts/s.

   #use uf cuts with additional cuts in filter_pars

Updates Dec 18, 2019
-	Change Dohs_List to Dohs_Single, now runs one target per call
- Remove code that filters on longer timescales
-	Remove calls to code for all plots
-	Update file and directory naming conventions
-	More easily readable format of filt and light curve definitions in report and fits output
-	Update ratekey variable definition
-	UF processing report now in pdf form, compressed
-	Compress/rename att file instead of bus file
-	Update value of CREATOR keyword
-	Update names of python file dependencies

updates Jan 9, 2020
- update value of CREATOR keyword = db_HSuf
- remove 'hs' from directory structure
- change report_file writing method to append to be able to read where code stopped
- remove observation number feature (not needed)
- streamline output to report_file

updates Oct 9, 2022 - JKB
- conversion to Python 3.
- removed unused numpy import.
-time.clock has been deprecated as of py 3.3 - recommended change is to time.process_time, this change was made.
'''

import sys, os, time, fnmatch, glob
#import numpy as np
    
#find list of data files, remove ending for dohs
def datalist(pattern,data_dir):
    curpath = os.getcwd()  #get current working directory
    os.chdir(data_dir)    #change to dir where data is stored
    filelist=os.listdir('.')   #get a list of files
    filelist.sort()   
    specarr, counter = [], 0
    file_end = pattern.split('*')[1]
    for entry in filelist:  #find files that match the pattern given
        if fnmatch.fnmatch(entry,pattern):
            infile=entry.split(file_end)[0] #remove extension
            specarr.append(infile)
            counter+=1  #collect number of files
    os.chdir(curpath)  #change back to initial dir
    return specarr, counter   #return file list and number of files
'''
PARAMETER DESCRIPTIONS:
source_list = [199] # or whatever hsnames you have
main_dir = 'E:/HaloSat/halosat_sourceanalysis/' #dir above where tree goes  
data_dir = 'E:/HaloSat/halosat_sourceanalysis/AllSources_v2p1/' #dir where raw files are
ancillary_dir = 'E:/HaloSat/Server_Code/' #dir where code is stored  
cut_note = 'uf'    #default cut is unfiltered
loud = 0            #only the standard output
creator = 'db_HSuf'    #name of software for file headers
timecuts=[]         #default is no time cuts
obs = ''            #default is to not add observation group # to names of output files
'''

def Dohs_Single(h,main_dir='/home/halosat/',\
                  data_dir='/home/halosat/v2p1/',\
                  ancillary_dir='/home/halosat/Server_Code/',\
                  timecuts=[],loud=0,creator='db_hsuf'):  
    #'''
    # directory containing HaloSat files such as response, target list, Python modules
    sys.path.append(ancillary_dir)
    import filt_dictionary3 as fd
    import halosat_analysis7_uf as hs
    #from numpy import array     #UNUSED
    
    #start timer, initialize parameters
    start_time = time.process_time()
    lcpha = fd.lcdef()  #default definition of PHA ranges for light curves in 8s bins matching HK times  
    if (h > 1000): cut_note='uf_de'  #change nadir angle cut for dark earth fields
    else: cut_note='uf'      
    hsname = str('HS{:04d}'.format(h))  #set string corresponding to target number
    source_num = hsname.split('HS')[1]+'01'
    
    #check that searchdb output exists, get file root name
    file_check, n_files = datalist(hsname+'*.evt',data_dir)
    if n_files==0:
        print(str('{} file(s) found for field {}'.format(n_files,hsname)))
        return 0  #should skip to next field if file doesn't exist        
    if n_files>1: inroot=file_check[-1] #choose most recent file
    else: inroot=file_check[0]
    
    #set directory structure and output filename roots
    out_dir = main_dir+'v1/'+source_num+'/unfiltered/'  #heasarc dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    local_dir=main_dir+'v1_local/'+source_num+'/unfiltered/'  #local dir
    if not os.path.exists(local_dir): os.makedirs(local_dir)   
    outroot, hea_root = local_dir+'hs'+source_num, out_dir+'hs'+source_num  #add hs to beginning of file names HEA ROOT UNUSED
    rfile=outroot+'_uf.txt'
    
    #include searchdb output in report
    f = open(data_dir+inroot+'.txt','r')
    contents=f.read()
    f.close
    with open(rfile,'w') as report_file: #begin file
      report_file.write('Searchdb Output:\n----------------\n\n{}'.format(contents))
      report_file.write('\nUF Processing output:\n------------------')
    
    #set filter expression, print lcpha and formatted filter expressions to file
    lcpha_format='\nLight Curve band definitions [band name, low ADU, max ADU]:\n'
    for l in range(len(lcpha)): lcpha_format+='{0['+str(l)+']}\n'
    filt=fd.filt_def(cut_note).replace('  ','') #use uf cuts with no cuts added on 8s bins, remove double spaces
    f_arr=filt.split(' & ')    #split into separate filtering expressions
    with open(rfile,'a') as report_file: 
      report_file.write(lcpha_format.format(lcpha))    #write lcpha to log file
      report_file.write('\nCut Type = {0}\nCut expression (category restriction):\n'.format(cut_note)) 
    heasarc_filt = '"'
    for l in range(len(f_arr)):
      key_math=f_arr[l].split('(hk["')[1].split('"] ')  #remove hk wrapper from expression
      filt_el=[key_math[0],key_math[1].rstrip().rstrip(')')]  #remove trailing ) and store pieces
      with open(rfile,'a') as report_file: report_file.write('{0[0]} {0[1]}\n'.format(filt_el)) #write line to file
      heasarc_filt+=str('({0[0]} {0[1]}) & '.format(filt_el))  #exclude hk wrapping from expression      
    heasarc_filt=heasarc_filt[0:-3]+'"'  #remove last spaces and & and append a "
    
    # find pointing using target ID in file name
    # the file name must have the format HSnnnn where nnnn is the target number
    pointing, Target_ID, official_name, Object_type, Target_Name = hs.find_pointing(data_dir+inroot, ancillary_dir=ancillary_dir)
    with open(rfile,'a') as report_file:
      report_file.write('\nAnalyzing Target ID: {}, Offical Name: {}, Object Type: {}, HaloSat Name: {}'.\
          format(Target_ID, official_name, Object_type, Target_Name))
      report_file.write('\nTarget coordinates in [RA,DEC] are {}'.format(pointing))        
      report_file.write('\nInput file root: {}'.format(inroot))     

    #retrieve prefilter data, write obs groups to report file, load att data into dict
    pf_data = hs.prefilter_gen(pointing, data_dir+inroot, outroot, rfile) #prefilter dictionary
    bus = hs.load_att(data_dir+inroot)  #attitude dictionary

    # loop over all DPUs
    hk_dpu = {} #parent dictionary to store each hk dictionary
    for dpu in [14, 38, 54] :
      with open(rfile,'a') as report_file: 
        report_file.write('\n\nReducing data for DPU {}'.format(dpu))
        if timecuts!=[]: report_file.write('\nTAI Time interval(s) queried [start,stop]: {}'.format(timecuts))
      t1=time.process_time()

      # accumulate and filter event data, store important values in hk
      hk = hs.make_filter(data_dir+inroot, bus, dpu, filt, lcpha=lcpha, ancillary_dir=ancillary_dir, timecuts=timecuts,
                          pointing=pointing, loud=loud)          
      if ('TIME_DUR' not in hk.keys()): 
          hk_dpu[dpu]=hk
          continue   #skip dpu if no data
      hk['OBS_ID'], hk['POINTING'], hk['OBJECT'], hk['OBJTYPE'], hk['COMMENT']  = source_num, pointing, official_name, Object_type, Target_Name
      hk['FILT'], hk['HEASARC_FILT'] = filt, heasarc_filt
      hk['TOTAL_TIME']=sum(hk['TIME_DUR'])  #TIME_DUR now has a max of 8.025 s
                
      # print statistics
      with open(rfile,'a') as report_file:
        report_file.write('\nNumber of time bins = {}; number of bus records = {}'.format(len(hk['TIME']),hk['NBUS']))
        report_file.write('\nTotal exposure (s) = {:.3f}, in-band events = {}'.format(\
                hk['TOTAL_TIME'], sum(hk['LC_SCI']*hk['TIME_DUR'])))
        report_file.write('\n{:.3f} % time removed by {} cut.'.format(\
                (hk['TOTAL_TIME']-hk['ONTIME'])/hk['TOTAL_TIME']*100.,cut_note))  
      if hk['ONTIME']<8.:  #if no data left for dpu, move on
          print('No remaining GTIs after cuts for DPU', dpu)
          hk_dpu[dpu]=hk
          continue          
      with open(rfile,'a') as report_file:  #print statistics after filtering if data is left
        report_file.write('\nAfter filtering exposure (s) = {:.3f}, in-band events = {}'\
            .format(hk['ONTIME'], sum(hk['LC_SCI'][hk['GOOD_IDX']]*hk['TIME_DUR'][hk['GOOD_IDX']])))         
 
      #add prefilter data to hk, write evt and hk files, store hk dictionary for later return
      hk = hs.prefilter_hk(hk,pf_data,outroot)  #add prefilter data to hk dictionary
      # tempcor = True applies baseplate temperature correction for gain
      hk = hs.write_evt(hk, outroot, out_dir, tempcor=True, creator=creator, \
                        ratekey=['LC_SCI', 900, 2250], evt_type='_uf')  
      hs.write_hk(hk, outroot, out_dir, creator=creator, loud=0)  #only generated for uf files, so no evt_type
      hk_dpu[dpu]=hk      #store the hk for the current DPU in a parent dictionary

      #add processing time per detector to report
      with open(rfile,'a') as report_file: 
        report_file.write('\nProcessing time for DPU {}: {:.3f}s'.format(dpu,time.process_time()-t1))
    
    #remove prefilter files from local uf directory
    pfnames = glob.glob(outroot+'*.mkf')
    for pfname in pfnames : os.remove(pfname)
        
    #add results from fchecksum and fverify on all files into report (one report for all detectors, diff for each ID#
    ontime=[hk_dpu[14]['ONTIME'],hk_dpu[38]['ONTIME'],hk_dpu[54]['ONTIME']]
    if min(ontime)>0.:  #only proceed if data survived cuts for all three DPUs
        #"""comment out on windows machine
        #check files with fverify and fchecksum, add fverify output to report file 
        HSlist, n_files = datalist(hsname+'*.att',data_dir) #returns root of filename
        HSlist=HSlist[0]
        if os.path.exists(outroot+'.att'): os.remove(outroot+'.att')
        os.system('cp '+data_dir+HSlist+'.att '+outroot+'.att')
        gz_list = []
        with open(outroot+'_filelist.txt','w') as flist:  #create list of files for fchecksum and fverify
          flist.write('{}.att'.format(outroot))
          gz_list.append(str('{}.att'.format(outroot)))
          for dpu in [14,38,54]:
            flist.write('\n{}_s{}_uf.evt'.format(outroot,str(dpu)))
            gz_list.append(str('{}_s{}_uf.evt'.format(outroot,str(dpu))))
            flist.write('\n{}_s{}.hk'.format(outroot,str(dpu)))
            gz_list.append(str('{}_s{}.hk'.format(outroot,str(dpu))))
        os.system('ftchecksum @'+outroot+'_filelist.txt update=yes datasum=yes chatter=1')  #checksum all files at once
        os.system('ftverify @'+outroot+'_filelist.txt outfile='+outroot+'_fverify.txt clobber=yes')   
        os.remove(outroot+'_filelist.txt')
        
        #copy and gzip to heasarc dir, add fverify output to report file
        with open(rfile,'a') as report_file:
          report_file.write('\n\nFtchecksum Results:\n--------------------')        
        for f in gz_list:
          #print f
          filename = f.split('/')[-1]
          hs.print_checksum(f,rfile)  #should error if checksum failed, prints keywords to report file
          if os.path.exists(out_dir+filename+".gz"): os.remove(out_dir+filename+".gz")
          os.system("imcopy "+f+" "+out_dir+filename+".gz")          
        f = open(outroot+'_fverify.txt','r')
        contents=f.read()
        f.close
        os.remove(outroot+'_fverify.txt') #remove file
        with open(rfile,'a') as report_file:
          report_file.write('\n\nFtverify Results:\n--------------------\n{}'.format(contents))
    else: 
        print('No data surviving for at least one detector for this field. Please check your cuts.')
        with open(rfile,'a') as report_file:
          report_file.write('No data surviving for at least one detector for this field. Please check your cuts.')
    
    #add processing time to report
    red_time=time.process_time()-start_time
    with open(rfile,'a') as report_file:
      report_file.write('\nData reduction time: {:.3f}s'.format(red_time))
    
    #print str('Data reduction time: {:.3f}s'.format(red_time))
    if (n_files>0): return hk_dpu  #return parent hk dictionary if searchdb output was found
    else: return 0
