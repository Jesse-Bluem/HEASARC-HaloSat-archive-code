'''
# Routines for HaloSat data processing and filtering
# P. Kaaret 11/2/2019-5/2/2019, 5/31/2019
# R. Ringuette 1/22/2019, 4/25/2019, 5/1/2019
# this version returns info about duplicate removal from make_filter
#newfilter: changes lc cuts to be applied on 64s time scale
# v6: R. Ringuette 8/28/2019: implement index mapping, x3 speed improvement with 
#    same exposure time. Duplication no longer necessary. Also, no event times
#    with no hk times. Standardize output. Code altered has RR at the end of the line.
# v6_uf2: R. Ringuette Dec 17, 2019: More efficient way to calculate GTIs,
#    bug fixes (outputs now identical to Phil's fixed version), update keywords
#    output structure.
Updates Dec 17, 2019
-	Add print_checksum function for more efficient output of the CHECKSUM and DATASUM values to uf report
-	Remove plotting functions
-	Change calls to numpy package
-	Update boresight vector based on Crab slews
-	More efficient storage of arrays from hk and bus files
-	Remove gti_min requirement and filtering from code
-	Add calculation of position and velocity vectors from bus data
-	Remove extraneous code
-	Change TIME_DUR keyword default value to 8.0 s
-	Move PHA filtering to index map creation section
-	Correct index map creation error
-	More efficient GTI calculation method
-	Update keywords in file headers and extensions
-	Add SCREENING extension to hk and evt files

Updates Jan 6, 2020
- Add and update prefilter code from Phil
(all repeated data is within <1 degree of each other)
- Add prefilter outputs to hk file
- Update SCREENING extension format

halosat_analysis6_uf4.py (Jan 8 to Jan 20ish, 2020)
- remove SCREENING extension for evt and hk files
- update keywords as noted in Jan8 meeting summary
- remove position and velocity calculation from bus data
- add checksum, datasum, and fverify output to print_checksum routine
- remove mfk files after prefilter data is saved
- use att files instead of bus files
- remove X and Y calculations and all functions used for this
- test ANTISUN_DLON = angular difference between ecliptic longitude and anti-sun direction 

halosat_analysis7_uf.py (Jan 31, 2020)
- add SUN_RA and SUN_DEC to hk output
- all headers checked with HEASARC (and Lorella)
- better IAU option 2 naming code
- add COMMENT keyword for numerical offset targets
- adjust pi value calc and limits (lines 532, 588, 589)
- NOTE: better to establish data columns line by line rather than in a loop.
    It seems the fits headers associated with the data are arranged by order of creation,
    so the units and null values created after the columns are not placed properly unless 
    it is established with the data column. Must use evthdu.header.insert(hdu#,(label,value,comment)).

Updates Oct 9, 2022 - JKB
- Updated to Python 3.
- Commented out unused imports and variables.
- Line 561 -1 to channel (gain fix).
- DATE-OBS and DATE-END have 32 seconds added to them.
- time.clock changed to time.process_time (Python 3 deprecates time.clock)
- (Jan 2023) added axis=0 to np.append calls. Without axis=o the 3D position and velocity vectors were flattened into a list, resulting in bad pos/vel data in the end hk FITS files.
'''

#from numpy import ndarray UNUSED
import numpy as np
import astropy.io.fits as pyfits
import time, datetime, os, glob
from astropy import units as u
from astropy.coordinates import SkyCoord
#from astropy.time import Time UNUSED

#getDate gets a datetime format from a TAI seconds from 2000/01/01
def getDate(time):  #from Daniel
  import pandas as pd
  epoch_date = pd.to_datetime('2000/01/01 00:00:00')
  time1 = time - 37.0
  date = epoch_date + datetime.timedelta(seconds=time1)
  return date

#getTime gets a TAI seconds from 2000/01/01 from a datetime
def getTime(date):  #from Daniel
  import pandas as pd
  epoch_date = pd.to_datetime('2000/01/01 00:00:00')
  date1 = pd.to_datetime(date)
  time = (date1-epoch_date).total_seconds()
  time += 37
  return time


#return header values for CHECKSUM and DATASUM calculation results, R. Ringuette
def print_checksum(check_file,rfile):
    #collect checksum and datasum values from primary header of file
    t = pyfits.open(check_file)
    h0 = t[0].header
    t.close()
    check_string = str('File = {},\nCHECKSUM = {}, DATASUM = {}.'.format(check_file,h0['CHECKSUM'], h0['DATASUM']))
    with open(rfile,'a') as report_file:  #print to file
      report_file.write('\n{0}'.format(check_string))    
    return  

#used on in OFFSET calculation
def distance(p1, p2):  #p1 = RA1, DEC1; p2=RA2, DEC2
# computes the angular distance in arcsec between two objects
# great circle angular distance between two coordinates    
# p1 and p2 are each two component vectors of sky coordinates
# [ra, dec] both in degrees
# returns angular distance in degrees
  # unpack vectors of sky coordinates
  (r1, d1) = p1[0], p1[1]
  (r2, d2) = p2[0], p2[1]
  # convert first position to 3-vector
  theta = (np.pi/180)*(90.0-d1)
  phi = (np.pi/180)*r1
  x1 = np.cos(phi)*np.sin(theta)
  y1 = np.sin(phi)*np.sin(theta)
  z1 = np.cos(theta)
  # convert second position to 3-vector
  theta = (np.pi/180)*(90.0-d2)
  phi = (np.pi/180)*r2
  x2 = np.cos(phi)*np.sin(theta)
  y2 = np.sin(phi)*np.sin(theta)
  z2 = np.cos(theta)
  # find dot product
  t = x1*x2 + y1*y2 + z1*z2
  dist = np.arccos(t)*(180/np.pi)
  return dist

#from Anna Zajcyk, modified from bct_time_to_local
def tai_time_conv(event_time, reftime='UTC', refout='DB'):
    from datetime import datetime, timedelta
    from pytz import utc

    # convert BCT times to unix times then to Central US time
    # BCT EDU time is seconds since Epoch 2000-01-01 TAI time (GPS time)
    # epoch_linux = datetime(1970, 1, 1, 0, 0, 0)   # reference for POSIX epoch
    # epoch_linux = utc.localize(epoch_linux)   # POSIX localized to UTC
    epoch_bct = datetime(2000, 1, 1, 0, 0, 0)    # reference epoch yr = 2000.0 for BCT EDU
    epoch_bct = utc.localize(epoch_bct)   # epoch localized to UTC
    if (reftime == 'UTC'): #default is utc
        leapseconds = 37.0 # GPS time is offset from UTC time because of leap seconds
    elif (reftime == 'tai'): 
        leapseconds = 0.0   # want time to stay in TAI
    utc_time=[]
    for i in range(len(event_time)):
      time_delta = timedelta(seconds=event_time[i]-leapseconds) # create datetime object from event_time [sec]
      true_time = epoch_bct + time_delta # get event time in UTC date:time
      utc_true_time = true_time.astimezone(utc) # true time should already be in UTC, but just making sure...
      if (refout == ''):  #default is 'hea'
          out_time = utc_true_time.strftime('%Y-%m-%d %H:%M:%S')
      elif (refout == 'tru'):
          out_time = utc_true_time.strftime('%m/%d/%Y %H:%M:%S')
      elif (refout == 'DB'):
          out_time = utc_true_time.strftime('%Y-%m-%dT%H:%M:%S')
      utc_time.append(out_time)
    utc_time=np.array(utc_time,dtype=str)
    return utc_time

def find_pointing(inroot, ancillary_dir='') : # output matches file values for HS0079
# try to figure out the target ID from the input file name
# and then get the coordinates from the HaloSat targets file
  # get target ID from inroot file name
  # assumes name is of the form ...HSnnnn... where nnnn is target number
  try : 
    i = inroot.find('HS') # find HS in inroot file name
    Target_ID = int(inroot[i+2:i+6]) # extract following digits
  except :
  	print('Cannot find Target_ID from file name ', inroot)
  	return []
  # read Halosat target library
  td = pyfits.getdata(ancillary_dir+'haloSat_targets.fits', 1)
  q = np.where(td['Target_ID'] == Target_ID)[0]
  Target_Name=td['Target_Name'][q][0].replace(' ','_') #replace spaces with underscores 
  object_type = td['Type'][q][0].upper()  #all caps for HEASARC
  if len(q) == 0 :
    print('Target ID ', Target_ID, ' not found in ', ancillary_dir+'haloSat_targets.fits')
    return [],'',''
  pointing = [td['RA'][q][0], td['DEC'][q][0]]
  
  #change Target_Name to IAU convention option 2 if not a common name
  #all common names are IDs 1-12, Lockman Hole is 63, common name offsets are 344 to 354
  if (Target_ID>12) and (Target_ID!=63) and ((Target_ID<344) or (Target_ID>354)):  
    crd = SkyCoord(ra=pointing[0]*u.degree, dec=pointing[1]*u.degree, frame='icrs',distance=1*u.Mpc) #new distance
    crd2 = crd.to_string('hmsdms')
    crd3=crd2.split("+")  #from Jesse (thru except statement)
    try:
        crd4=crd3[1] #plus dec
        crd5='+'+crd4[0:2]+crd4[3:5]#+crd4[6:8]
    except:
        crd4=crd3[0].split("-")[1] #minus dec
        crd5='-'+crd4[0:2]+crd4[3:5]#+crd4[6:8]
    official_name = 'HaloSat J'+crd2[0:2]+crd2[3:5]+crd5  #+crd2[6:8]    
    print(official_name)
    #Target_Name = 'HaloSat J'+crd2[0:2]+crd2[3:5]+crd2[6:8]+crd2[15:18]+crd2[19:21]+crd2[22:24]
  else: official_name = Target_Name 
  
  return pointing, Target_ID, official_name, object_type, Target_Name

def load_att(inroot,timecuts=[]):
  f = pyfits.open(inroot+'.att') # read in the FITS file with bus data
  d = f[1].data # put all of the data into the dictionary d
  f.close() # close the FITS file
  bus = {} # make empty dictionary
  #bus_dpu = np.array(d['DPU_ID']) # find the DPU ID for each bus record
  if len(timecuts) == 0 : # no time cuts
    #  q1 = (bus_dpu == dpu) # pick out only the HK records for the selected DPU
    for k in d.names : # loop over all columns in the FITS record array
      bus[k] = np.array(d[k]) 
  else : # have time cuts
    bus_time0 = np.array(d['TIME']) # find the time for each bus record
    q0 = np.zeros(len(bus_time0), dtype=bool) # fill array with false
    for t0, t1 in timecuts : # loop over selected time intervals
      q0 = q0 | ((t0 <= bus_time0) & (bus_time0 <= t1)) # set true for records in this time interval
    # pick out only the HK records for the selected DPU and within a time interval
    q2 = np.unique(np.array(d['TIME'][q0]), return_index=True)[1]
    for k in d.names : # loop over all columns in the FITS record array
      bus[k] = np.array(d[k][q0][q2])
  return bus

def prefilter_gen(pointing, inroot, outroot, rfile):
  """
  ->inroot has data_dir, outroot has local dir
  ->run separately per det in case one fails
  Read the HaloSat spacecraft data (.bus), 
  group sets of pointing into observations, and generate a filter
  file for each observation using the ftool prefilter.
  Need to copy the following files to your hs_dir or local directory:
  halosat_tle.txt rigidity.date leapsec.fits
  """
  hs_dir = '/home/halosat/hscore/'
  # load spacecraft att gti data into arrays
  f = pyfits.open(inroot+'.att') # read in the FITS file with bus data
  gti = f['GTI'].data # put the GTI data into the dictionary gti
  f.close() # close the FITS file

  # find contiguous sets of pointings and group them into observations
  obs_start = [gti['START'][0]]
  obs_end = [gti['STOP'][0]]
  obs_gap = 86400.0 # maximum gap between pointings in an observation
  for s, e in np.transpose([gti['START'][1:], gti['STOP'][1:]]) :
    if (s - obs_end[-1]) < obs_gap : # if new GTI is soon enough after end of current observation
      obs_end[-1] = e # move the observation end to the end of the GTI
    else :
      obs_start.append(s) # otherwise start a new observation
      obs_end.append(e)
  with open(rfile,'a') as report_file: 
    report_file.write('\n\n{} observation groups found.'.format(len(obs_start)))
    report_file.write('\nObs#  Start (TAI)   Stop (TAI)   Start (UTC)          Stop (UTC)')
  for i in range(len(obs_start)) : 
    #print 'Observation ', i, ': ', obs_start[i], ' to ', obs_end[i]
    start_date, end_date = tai_time_conv([obs_start[i]])[0], tai_time_conv([obs_end[i]])[0] #convert to UTC date and time
    with open(rfile,'a') as report_file: 
      report_file.write('\n{:02d}    {}   {}  {}  {}'.format(i+1, obs_start[i],obs_end[i],start_date,end_date))
    # generate command string to run prefilter
    # start early so that SAA_TIME is valid at start
    s = 'prefilter attname='+inroot+'.att'+ \
        ' start='+str(obs_start[i]-2700.0)+' end='+str(obs_end[i])+ \
        ' interval=8.0 ranom='+str(pointing[0])+' decnom='+str(pointing[1])+ \
        ' outname='+outroot+'_'+str(i).zfill(2)+'.mkf columns=ALL orbmode=TLE_TEXT2'+ \
        ' orbname='+hs_dir+'halosat_tle.txt alignfile=NONE'+ \
        ' leapname='+hs_dir+'leapsec.fits'+ \
        ' rigname='+hs_dir+'rigidity.data missepoch=2000-01-01T00:00:00'+ \
        " timeadj='CONST:-37.0'"+ \
        ' origin=NASA/GSFC clobber=yes chatter=0'
    #print s
    os.system(s)
  
  #store chosen prefilter data in dictionary  
  prefilter_list=np.array(['SAT_ALT','ELV','BR_EARTH','FOV_FLAG','SUNSHINE','SUN_ANGLE','MOON_ANGLE',\
    'RAM_ANGLE','ANG_DIST','COR_ASCA','COR_SAX','MCILWAIN_L','SAA','SAA_TIME','RA','DEC',\
    'SAT_LAT','SAT_LON','POSITION','VELOCITY','ROLL','SUN_RA','SUN_DEC'],dtype=str)  #prefilter does not have NADIR_ANGLE  
  pf_data={}
  pfnames = glob.glob(outroot+'*.mkf') #outroot has the directory tree in it
  for pfname in pfnames : #store data from each file
    f = pyfits.open(pfname) # read in the prefilter file
    d = f[1].data # put all of the data into the FITS record array d
    f.close() # close the FITS file
    
    # fill arrays for the various quantities of interest with data
    if len(pf_data)>1: 
      for k in prefilter_list: pf_data[k] = np.append(pf_data[k],np.array(d[k]),axis=0) #JKB axis=0 fix for position and velocity vector bug
      pf_data['TIME'] = np.append(pf_data['TIME'],np.array(d['TIME']))
    else: 
      for k in prefilter_list: pf_data[k] = np.array(d[k])  #store chosen data in arrays
      pf_data['TIME'] = np.array(d['TIME']) # time      
  return pf_data  
  
def prefilter_hk(hk, pf_data, outroot):   # load chosen data from prefilter file(s) into hk dictionary
  #initialize dictionary arrays
  dt = 6.0 # maximum allowed time offset between HK and prefilter records
  n = len(hk['TIME']) # number of time bins in HK
  prefilter_list=np.array(['SAT_ALT','ELV','BR_EARTH','FOV_FLAG','SUNSHINE','SUN_ANGLE','MOON_ANGLE',\
    'RAM_ANGLE','ANG_DIST','COR_ASCA','COR_SAX','MCILWAIN_L','SAA','SAA_TIME','RA','DEC',\
    'SAT_LAT','SAT_LON','POSITION','VELOCITY','ROLL','SUN_RA','SUN_DEC'],dtype=str)  #prefilter does not have NADIR_ANGLE
  integer_index=[3,4,12]  #indices for integer arrays
  decimal_index=[0,1,2,5,6,7,8,9,10,11,13,14,15,16,17,20,21,22]  #indices for decimal arrays
  index_3D = [18,19]  #indices for 3D arrays
  hk_pf={}  #initialize prefilter data storage
  for k in prefilter_list[decimal_index]: hk_pf[k]=np.zeros(n,dtype=float)-200.  #all decimal arrays set to -200 default
  for k in prefilter_list[integer_index]: hk_pf[k]=np.zeros(n,dtype=int)-999   #all integer arrays set to -999 default
  for k in prefilter_list[index_3D]: hk_pf[k]=np.zeros((n,3),dtype=float)-200.   #position and velocity arrays

  #translate prefilter data to hk time 
  pf_time = pf_data['TIME'] # time
  for i in range(n) :  # translate prefilter data to HK time bins
    q = np.argmin(abs(hk['TIME'][i] - pf_time))
    if abs(hk['TIME'][i] - pf_time[q]) <= dt :
      for k in prefilter_list:
        hk_pf[k][i]=pf_data[k][q]    
      
  # copy prefilter data into hk dictionary
  for k in prefilter_list:
    if k in hk: hk[k+'P']=hk_pf[k]  #if key repeated, add 'P' to key
    else: hk[k]=hk_pf[k]      #else, just store in dictionary
     
  return hk

def make_filter(inroot, bus, dpu, filt, timecuts=[], lcpha=[], loud=0, \
                ancillary_dir='', pointing=[], evtfiles=[]) :
  """
  Read HaloSat data from housekeeping (.hk), spacecraft (.bus), and
  event (.evt) files.  Convert into filter data with uniform time bins.
  Returns dictionary of data for filtering, including light curves.
  """

  hk = {} # make empty dictionary for housekeeping and pointing/orbit data
  hk['GOOD_DPU']=dpu  #store dpu number
  # load instrument housekeeping data into arrays
  f = pyfits.open(inroot+'.hk') # read in the FITS file with HK data
  d1 = f[1].data # put all of the data into the FITS record array d
  f.close() # close the FITS file
  hk_dpu = np.array(d1['DPU_ID']) # find the DPU ID for each HK record
  #only keep times within time range given, altered from bus file section, RR
  if len(timecuts) == 0 : # no time cuts
    q1 = (hk_dpu == dpu) # pick out only the HK records for the selected DPU
  else : # have time cuts
    hk_time0 = np.array(d1['TIME']) # find the time for each hk record
    q0 = np.zeros(len(hk_time0), dtype=bool) # fill array with false
    for t0, t1 in timecuts : # loop over selected time intervals
      q0 = q0 | ((t0 <= hk_time0) & (hk_time0 <= t1)) # set true for records in this time interval
    # pick out only the HK records for the selected DPU and within a time interval
    q1 = (hk_dpu == dpu) & q0  
  q2 = np.unique(np.array(d1['TIME'][q1]), return_index=True)[1] # find unique times of HK records for the selected DPU
  #q2 = t[1] # indices of the first occurrence of each time
  # for each key in d, extract only entries for this DPU and with unique times
  for k in d1.names : # loop over all columns in the FITS record array
    #t = d1[k][q1]
    hk[k] = np.array(d1[k][q1][q2])
    #hk[k] = np.array(t[q2])  
  if (loud > 2) : print('HK time range = ', min(hk['TIME']), max(hk['TIME']))
  #output meaningful error messages, RR
  if len(hk['TIME']) < 1 : # no time bins
    print('No time bins found in housekeeping data')
    hk['ONTIME']=0.0
    return hk
  time_dur = np.diff(hk['TIME'])  #calculate difference between hk times, RR
  hk['TIME_DIFF'] = np.insert(time_dur,0,0.)  #save time duration between each hk_time, insert 0 to match GTI array structure, RR

  # generate time bins based on Hk times
  dt = 8.0 # time between HK records !!! update if we change the interval
  jitter = 0.1 # account for jitter in times

  # translate bus data to HK time bins
  hk['NBUS'] = len(bus['TIME'])  #store for printing to log file
  #hk['BUS_TIME'] = bus['TIME']  #testing consecutive att times
  hk_time = hk['TIME'] # put time into an array for convenvience
  n = len(hk_time) # number of time bins
  #print 'Number of time bins = ', n, ' number of bus records = ', len(bus_time)
  ra, dec = np.zeros(n), np.zeros(n) # arrays to hold ra, dec
  sat_lat, sat_lon = np.zeros(n), np.zeros(n) # arrays to hold latitude, longitude
  nadir_angle = np.zeros(n) # array to hold nadir angle
  #x, y = np.zeros(n), np.zeros(n) # arrays to hold x, y UNUSED
  for i in range(n) :
    #q = np.where( (hk_time[i] <= bus_time) & (bus_time < (hk_time[i]+dt)) )[0]  #changed Dec9 to be [0]
    q = np.where(abs(hk_time[i] - bus['TIME']) <= 0.5* dt)[0]  #bus['TIME']=bus_time
    if len(q) > 0 :
      ra[i], dec[i]  = np.mean(bus['RA'][q]), np.mean(bus['DEC'][q])
      sat_lat[i], sat_lon[i] = np.mean(bus['SAT_LAT'][q]), np.mean(bus['SAT_LON'][q])
      nadir_angle[i] = np.mean(bus['NADIR_ANGLE'][q])
    else :
      ra[i], dec[i] = -200.0, -200.0
      sat_lat[i], sat_lon[i] = -200.0, -200.0
      nadir_angle[i] = -200.0

  # copy RA, DEC into main dictionary !!! add new variables
  hk['RA'], hk['DEC'] = np.array(ra), np.array(dec)
  hk['SAT_LAT'], hk['SAT_LON'] = np.array(sat_lat), np.array(sat_lon)
  hk['NADIR_ANGLE'] = np.array(nadir_angle)
  # only find pointing offset if pointing is specified
  if len(pointing) == 2 :
    # calculate the pointing offset
    hk['OFFSET'] = np.zeros(n)
    for i in range(n) :
      p = [ra[i], dec[i]] # pointing of instrument bore sight
      hk['OFFSET'][i] = distance(pointing, p) # offset of bore sight from target

  # find times when we are in SAA boxes
  f = pyfits.open(ancillary_dir+'saa_boxes.fits') # read in FITS file with SAA boxes
  d = f[1].data # put all of the data into the variable d
  f.close() # close the FITS file
  # move the data from d to numpy arrays
  saa_valid = np.array(d['SAA_VALID']) # TAI Times of when each set of SAA boxes is valid
  sb1 = np.array(d['SAA_BOX1']) # longitude [min, max], latitude [min, max]
  sb2 = np.array(d['SAA_BOX2']) # longitude [min, max], latitude [min, max]
  saa_start = saa_valid[:,0] # start time of each SAA box valid interval
  saa_stop = saa_valid[:,1] # end time of each SAA box valid interval Dec9
  in_saa = np.zeros(n, dtype=bool) 
  for i in range(n) : 
    lon, lat = hk['SAT_LON'][i], hk['SAT_LAT'][i]
    w = np.where((saa_start <= hk_time[i])&(saa_stop >= hk_time[i]))[0] # which SAA box to use for this HK record Dec9
    in_saa[i] = ((sb1[w,0] <= lon) & (lon <= sb1[w,1]) & (sb1[w,2] <= lat) & (lat <= sb1[w,3])) | \
                ((sb2[w,0] <= lon) & (lon <= sb2[w,1]) & (sb2[w,2] <= lat) & (lat <= sb2[w,3]))
  hk['IN_SAA'] = in_saa # write the result to the dictionary

  # use a single event file from the database
  # load instrument event data into arrays
  f = pyfits.open(inroot+'.evt') # read in the FITS file with event data
  d = f[1].data # put all of the data into the dictionary d
  f.close() # close the FITS file
  q = (np.array(d['DPU_ID']) == dpu) # pick out only the events for the selected DPU
  hk['EVT_TIME'] = np.array(d['TIME'][q]) # event time
  hk['EVT_PHA'] = np.array(d['PHA'][q]) # event pulse height amplitude    
  t0 = time.process_time()
  if len(lcpha) > 0 :
    # find light curve(s) of event rate in HK time bins in PHA bands set by user in lcpha
    # we break the event list into pieces to make this go faster
    n = len(hk_time) # number of time bins
    time_chunk = 1500 # [s] seconds of data to process in each chunk (orbit length) Dec9
    if (loud > 1) : print('Making light curves, # time bins = ', n,)
    nlc = len(lcpha)
    lc = np.zeros((nlc,n)) # array to hold counts in time bins
    q = np.where((hk_time[0] <= hk['EVT_TIME']) & (hk['EVT_TIME'] < (hk_time[0]+time_chunk)))[0]
    et = hk['EVT_TIME'][q]
    ep = hk['EVT_PHA'][q]
    if len(et) > 0 : et_max = max(et) 
    else : et_max = min(hk['EVT_TIME'])
    idx_map=[]  #initialize map of indices, each position will correspond to the hk time at the same index, RR
    time_durp1 = np.append(hk['TIME_DIFF'],(dt))#+jitter))  #skip first (0.) value, append 8.1 at end to include end events, RR
    new_time_dur=[]
    for i in range(n) :
      if hk_time[i]+dt >= et_max :
        q = np.where((hk_time[i] <= hk['EVT_TIME']) & (hk['EVT_TIME'] <= (hk_time[i]+time_chunk)))[0]
        et = hk['EVT_TIME'][q]
        ep = hk['EVT_PHA'][q]
        if len(et) > 0 : et_max = max(et)
        else : et_max = hk_time[i]+time_chunk
      if (time_durp1[i+1]>(dt+jitter)): new_time_dur.append(dt)  #if break, use 8.0 seconds
      else: new_time_dur.append(time_durp1[i+1])  #otherwise, use time duration
      time_dur=new_time_dur[i]
      if (time_dur < jitter): #skip hk times within 0.1s (jitter) of another
          idx_map.append([[]])  #if length of time bin too short, append an empty list for this hk time, RR
          if (loud>1):  #print info on short hk times for trouble shooting, RR
              print(str('\nHK time {}s at index {} skipped with duration {:.3f}s.'.format(hk_time[i], i, time_dur)))
              print(str('BPL_TEMP = {:.3f}, DPU_TEMP = {:.3f}, TEC_PWM = {}.'.format(hk['BPL_TEMP'][i],hk['DPU_TEMP'][i],hk['TEC_PWM'][i])))
              print(str('Next HK time {}s at index {}'.format(hk_time[i+1], i+1)))
              print(str('with BPL_TEMP = {:.3f}, DPU_TEMP = {:.3f}, TEC_PWM = {}.'.format(hk['BPL_TEMP'][i+1],hk['DPU_TEMP'][i+1],hk['TEC_PWM'][i+1])))          
          continue   #go to next hk time, RR
      idx=np.where((hk_time[i] <= et) & (et < (hk_time[i]+time_dur)))[0]  #avoids double counting or missing events by using the actual time duration of the bin, RR          
      filt_idx=np.where((ep[idx]>=600) & (ep[idx]<=12900))[0]  # new Dec9
      for j in range(nlc) : 
        lc[j,i] = np.sum((lcpha[j][1] <= ep[idx]) & (ep[idx] < lcpha[j][2]) ) / time_dur  #more accurate count rate based on actual length of hk bin, RR
      if len(q)>0: idx_map.append([q[idx[filt_idx]]])  #if events are found in this hk time bin, add list of evt array indices to this position, RR Dec9
      else: #if no events found for this hk time bin, append an empty list, RR
          idx_map.append([[]])
    # copy the light curves into the main dictionary
    for j in range(nlc) : 
      hk[lcpha[j][0]] = lc[j,:]
    hk['IDX_MAP']=idx_map   #store map of event indices for each hk time for later use. RR
    hk['TIME_DUR']=np.array(new_time_dur)

  # evaluate filter expression to find times when data is 'good'
  good = eval(filt)
  good_idx=np.linspace(0,len(good)-1,len(good),dtype=int)
  good_idx=good_idx[good]  #save hk['GOOD'] definition in index form instead of boolean, RR (faster)
  # add good to data dictionary
  #hk['GOOD'] = good
  #hk['FILT_EXP'] = filt # save filter expression
  hk['GOOD_IDX'] = good_idx
  hk['ONTIME']=np.sum(hk['TIME_DUR'][hk['GOOD_IDX']])               #updated Dec13
  
  #find GTIs, no filtering on GTI length at this stage. new Dec13
  good_diff = np.diff(hk['TIME'][hk['GOOD_IDX']])   #difference between consecutive good times
  idx = np.where(good_diff>8.1)[0]    #find indices after which there are breaks in time
  hk['GTI_START'] = np.insert(hk['TIME'][hk['GOOD_IDX']][idx+1],0,hk['TIME'][hk['GOOD_IDX']][0]) #first good time, then times after breaks
  hk['GTI_STOP'] = np.append(hk['TIME'][hk['GOOD_IDX']][idx]+hk['TIME_DUR'][hk['GOOD_IDX']][idx],\
    hk['TIME'][hk['GOOD_IDX']][-1]+hk['TIME_DUR'][hk['GOOD_IDX']][-1]) 
  #times of breaks + duration of hk time bin, then last good time + duration of last good hk time bin  
  #identical exposure time to method a few lines up (hk['ONTIME'])
  #could filter by gti length by using the idx variable

  return hk

#write evt files separately from pi files
def write_evt(hk, outroot, main_dir, creator='', evt_type='_uf',loud=0, tempcor=True, ratekey=[]) :
  """
  Convert event pulse heights to pulse invariant (energy) values and write spectrum.
  select is a list of time intervals to accept
  """

  # copy event data from dictionary in to arrays
  dpu = hk['GOOD_DPU'] # DPU ID for each event, same for all
  evt_time = hk['EVT_TIME'] # event time
  pha = hk['EVT_PHA'] # event pulse height amplitude
  # make list of good events
  good_time = [] # time of events in GTIs
  good_pha = [] # pha of events in GTIs
  bt = [] # baseplate temp for each event, RR
  #use index mapping saved from earlier to match surviving hk times with events, RR
  for i in hk['GOOD_IDX']:  #this increases speed by 3x to 5x from previous method, RR
      q = hk['IDX_MAP'][i][0]  #collect event indices for current hk time, RR
      good_time.extend(evt_time[q]) # add event times to list for time, RR
      good_pha.extend(pha[q]) # add PHA of events in this GTI to pha list, RR
      bt.extend(np.repeat(hk['BPL_TEMP'][i],len(hk['IDX_MAP'][i][0])))   #repeat BPL_TEMP value for this hk time for all events in this time bin, RR
  good_time = np.array(good_time)
  good_pha = np.array(good_pha)
  bt = np.array(bt)
  hk['GOOD_TIME'] = good_time # time for each event
  hk['GOOD_PHA'] = good_pha # pulse height amplitude for each event

  # set range and step of pi values, 
  # pi starts at energy = 0, increases in steps of ebin, so energy = pi*ebin
  # this must match program that wrote FITS file
  # these definitions must match with the arf and rmf
  emin = 0.1 # [keV] minimum energy value
  #emax = 9.2 # [keV] maximum energy value
  ebin = 0.02 # [keV] energy bin size for pi
  #nbin = int((emax-emin)/ebin) # number of bins in histogram, 455
  #print 'Number of bins in spectrum = ', nbin
  
  # convert PHA to energy, then to pi value
  #pi_unit = str(ebin)+ 'keV' UNUSED
  # !!! should read in coefficients from CALDB
  # find the index into the conversion arrays for this DPU
  dpu_list = np.array([14, 38, 54]) # Note Anna's DPU order!!!
  j = np.where(dpu_list == dpu)[0][0] 
  if tempcor == False : # do not do correction for baseplate temperature
    # coefficients are temperature-averaged values from Anna's report of 6/8/2019
    # Group 1 coefficients for conversion of ADC channels to eV
    c1 = np.array([  0.54384,   0.55626,   0.54930])
    c2 = np.array([-30.78253, -33.78921, -36.55629]) 
    ## Group 2 coefficients for conversion of ADC channels to eV
    #c1 = np.array([0.54420, 0.55715, 0.54943])
    #c2 = np.array([ -31.24,  -33.91,  -37.90]) 
    #c3 = np.array([1.000, 1.020, 0.986]) # E = (C1*PHA + C2)/C3
    # convert PHA to energy, then to pi value
    good_keV = (c1[j]*good_pha+c2[j])/1000.0 # find energy of this event in keV
  else : # do correction for baseplate temperature
    #if (loud > 0) : print 'Finding HK record for each event. ',
    # coefficients are temperature-dependent values from Anna's report of 6/4/2018
    dpu_list = np.array([14, 38, 54])
    c1a1 = np.array([  0.0000018615,   0.0000023380,  0.0000006617])
    c1a2 = np.array([  0.0000058429,  -0.0000472094, -0.0000147105])
    c1a3 = np.array([  0.5437512278,   0.5564678190,  0.5493653560])
    c2a1 = np.array([ -0.1426338885,  -0.1355792938, -0.1809170943])
    c2a2 = np.array([-29.9819221014, -32.9267861288,-35.5069756003])
    c1 = c1a1[j]*bt**2 + c1a2[j]*bt + c1a3[j] # calculate energy calibration coefficient C1
    c2 = c2a1[j]*bt + c2a2[j] # calculate energy calibration coefficient C2
    #print 'DPU ', dpu, ' - ', 'temp=', np.mean(bt), '  c1=', np.mean(c1) , '  c2=', np.mean(c2)
    # convert PHA to energy, then to pi value
    good_keV = (c1*good_pha+c2)/1000.0 # find energy of this event in keV
  # convert energy to pi value    
  #good_pi = np.floor((good_keV - emin)/ebin).astype(int)
  good_pi = np.floor(1.+(good_keV - emin)/ebin).astype(int)-1 #(JKB) fix channel being offset by 1.
  hk['GOOD_PI'] = np.array(good_pi)

  #create time strings
  date_obs = tai_time_conv([min(hk['GTI_START'])+32])[0]  #UTC time
  date_end = tai_time_conv([max(hk['GTI_STOP'])+32])[0]   #UTC time
  file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')

  # write out filtered event file
  # make the primary header (standard in all files)
  prihdr = pyfits.Header()
  prihdr['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
  prihdr['INSTRUME'] = ('SDD'+str(hk['GOOD_DPU']), 'Instrument name') 
  prihdr['OBS_ID'] = (hk['OBS_ID'], 'Observation ID')
  prihdr['OBJECT'] = (hk['OBJECT'], 'Object/Target Name')
  prihdr['OBJTYPE'] = (hk['OBJTYPE'], 'Object/Target type')  
  prihdr['DATE-OBS'] = (date_obs, 'Start date of observations')
  prihdr['DATE-END'] = (date_end, 'End date of observations')
  prihdr['DATE'] = (file_time, 'File creation date')  #checksum and datasum are added with fverify etc
  prihdu = pyfits.PrimaryHDU(header=prihdr) # make the primary header

  # define the columns for events
  col1 = pyfits.Column(name='TIME', format='1D', unit='s', array=good_time) # event time
  col3 = pyfits.Column(name='PHA', format='1I', unit='chan', array=good_pha) # event PHA 
  col4 = pyfits.Column(name='PI', format='1I', unit='chan', array=good_pi) # event pi
  cols = pyfits.ColDefs([col1, col3, col4]) # create a ColDefs (column-definitions) object for all columns:
  evthdu = pyfits.BinTableHDU.from_columns(cols) # create a new binary table HDU object 
  
  #set data type keywords for event extension of evt file
  evthdu.header.comments['TTYPE1'] = 'Time of events'  #8th header
  evthdu.header.comments['TFORM1'] = 'data format of field'
  evthdu.header.comments['TUNIT1'] = 'physical unit of field'
  evthdu.header.comments['TTYPE2'] = 'Pulse Height Analyzer'
  evthdu.header.comments['TFORM2'] = 'data format of fields'
  evthdu.header.comments['TUNIT2'] = 'physical unit of field'
  evthdu.header.insert(14,('TLMIN2', 600, 'minimum legal value of column'))  #add column 2 headers
  evthdu.header.insert(15,('TLMAX2', 12900, 'maximum legal value of the column'))   
  evthdu.header.comments['TLMIN2'] = 'minimum legal value of column'
  evthdu.header.comments['TLMAX2'] = 'maximum legal value of the column'
  evthdu.header.comments['TTYPE3'] = 'Pulse Invariant'
  evthdu.header.comments['TFORM3'] = 'data format of field'
  evthdu.header.comments['TUNIT3'] = 'physical unit of field'
  evthdu.header.insert(21,('TLMIN3',1, 'minimum legal value of the column'))
  evthdu.header.insert(22,('TLMAX3',455, 'maximum legal value of the column'))
  evthdu.header.insert(23,('TNULL3',-1, 'null value'))
  
  #standard extension keys
  evthdu.header['PI2ENE'] = (0.02, 'PI conversion from chan to energy keV')
  evthdu.header['EXTNAME'] = ('EVENTS', 'Binary table extension name')  
  evthdu.header['HDUCLASS'] = ('OGIP', 'format conforms to OGIP/GSFC standards')
  evthdu.header['HDUCLAS1'] = ('EVENTS', 'First class level') 
  evthdu.header['HDUCLAS2'] = ('ALL', 'Second class level')  
  evthdu.header['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
  evthdu.header['INSTRUME'] = ('SDD'+str(hk['GOOD_DPU']), 'Instrument name') 
  evthdu.header['DATAMODE'] = ('PHOTON', 'Instrument datamode')
  evthdu.header['OBSERVER'] = ('PHILIP KAARET', 'Principal Investigator')  
  evthdu.header['OBS_ID'] = (hk['OBS_ID'], 'Observation ID')
  evthdu.header['OBJECT'] = (hk['OBJECT'], 'Object/Target Name')
  evthdu.header['OBJTYPE'] = (hk['OBJTYPE'], 'Object/Target type')  
  #if (float(hk['OBS_ID']) >=35501) & (float(hk['OBS_ID'])<37301):
  #  evthdu.header['COMMENT'] = (hk['COMMENT'], 'Offset target association')  
  evthdu.header['EQUINOX'] = (2000, '[yr] Equinox of celestial coord system')
  evthdu.header['RADECSYS'] = ('FK5', 'Celestial coord system')
  evthdu.header['RA_NOM'] = (hk['POINTING'][0], '[deg] R.A. of nominal aspect point [J2000]')  
  evthdu.header['DEC_NOM'] = (hk['POINTING'][1], '[deg] Dec. of nominal aspect point [J2000]') 
  evthdu.header['RA_OBJ'] = (hk['POINTING'][0], '[deg] Object Right ascension [J2000]')
  evthdu.header['DEC_OBJ'] = (hk['POINTING'][1], '[deg] Object Declination [J2000]')
  evthdu.header['TIMESYS'] = ('TT', 'Reference Time System')
  evthdu.header['MJDREFI'] = (51544, '[d] MJD reference day (2000-01-01T00:00:00)')
  evthdu.header['MJDREFF'] = (0.00074287037037037, '[d] MJD reference (fraction of day)')  #MJD
  evthdu.header['TIMEREF'] = ('LOCAL', 'Reference Frame')
  evthdu.header['TASSIGN'] = ('SATELLITE', 'Time assigned by clock')
  evthdu.header['TIMEUNIT'] = ('s', 'Time unit for timing header keyword')
  evthdu.header['TIMEDEL'] = (0.05,'[s] Data time resolution.')   
  evthdu.header['TIMEZERO'] = (0.0, '[s] Time Zero')
  evthdu.header['TIMEPIXR'] = (1, 'Bin time beginning=1 middle=0.5 end=1')
  evthdu.header['TIERRELA'] = (4.0E-06, '[s/s] relative errors expressed as rate')
  evthdu.header['TIERABSO'] = (1.0, '[s] timing precision in seconds')
  evthdu.header['TSTART'] = (min(hk['GTI_START']), '[s] Observation start time')
  evthdu.header['TSTOP'] = (max(hk['GTI_STOP']), '[s] Observation stop time')
  evthdu.header['TELAPSE'] = (max(hk['GTI_STOP'])-min(hk['GTI_START']), '[s] Stop-Start')  
  evthdu.header['ONTIME'] = (hk['ONTIME'], '[s] Observation time of source')
  evthdu.header['EXPOSURE'] = (hk['ONTIME'], '[s] exposure')      
  evthdu.header['DATE-OBS'] = (date_obs, 'Start date of observations')
  evthdu.header['DATE-END'] = (date_end, 'End date of observations')
  evthdu.header['CLOCKAPP'] = (True, 'Clock correction applied? (F/T)')
  evthdu.header['DEADAPP'] = (False, 'Has deadtime been applied to data? (F/T)')
  evthdu.header['ORIGIN'] = ('UNIVERSITY OF IOWA', 'Origin of fits file')
  evthdu.header['PROCVER'] = ('hsuf_20221026', 'Processing script version number')   #prev. hsuf_20200131_hscl_20200131
  evthdu.header['SOFTVER'] = ('Hea_11apr2022_V6.30.1', 'Software version')     
  evthdu.header['CALDBVER'] = ('hs20200129', 'CALDB index version used')   
  evthdu.header['TLM2FITS'] = ('db_20221002', 'Telemetry converter FITS number')   #location?!    
  evthdu.header['CREATOR'] = (creator, 'Software creator of the file')  
  evthdu.header['DATE'] = (file_time, 'File creation date')  
  #checksum and datasum are added with fverify and fchecksum later

  # define the columns for GTI extension
  col1 = pyfits.Column(name='START', format='1D', unit='s',array=hk['GTI_START'])
  col2 = pyfits.Column(name='STOP', format='1D', unit='s', array=hk['GTI_STOP'])
  cols = pyfits.ColDefs([col1, col2]) # create a ColDefs (column-definitions) object for all columns:
  gtihdu = pyfits.BinTableHDU.from_columns(cols) # create a new binary table HDU object 
  
  # add data type keywords for gti extension of evt file
  gtihdu.header.comments['TTYPE1'] = 'Start time'
  gtihdu.header.comments['TFORM1'] = 'data format of field'
  gtihdu.header.comments['TUNIT1'] = 'physical unit of field'
  gtihdu.header.comments['TTYPE2'] = 'Stop time'
  gtihdu.header.comments['TFORM2'] = 'data format of field'
  gtihdu.header.comments['TUNIT2'] = 'physical unit of field'

  #standard extension keys
  gtihdu.header['EXTNAME'] = ('GTI', 'Binary table extension name')  
  gtihdu.header['HDUCLASS'] = ('OGIP', 'format conforms to OGIP/GSFC standards')
  gtihdu.header['HDUCLAS1'] = ('GTI', 'First class level')  
  gtihdu.header['HDUCLAS2'] = ('ALL', 'Second class level')  
  gtihdu.header['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
  gtihdu.header['INSTRUME'] = ('SDD'+str(hk['GOOD_DPU']), 'Instrument name') 
  gtihdu.header['OBSERVER'] = ('PHILIP KAARET', 'Principal Investigator')   
  gtihdu.header['OBS_ID'] = (hk['OBS_ID'], 'Observation ID')
  gtihdu.header['OBJECT'] = (hk['OBJECT'], 'Object/Target Name')
  gtihdu.header['OBJTYPE'] = (hk['OBJTYPE'], 'Object/Target type')
  #if (float(hk['OBS_ID']) >=35501) & (float(hk['OBS_ID'])<37301):
  #  gtihdu.header['COMMENT'] = (hk['COMMENT'], 'Offset target association')  
  gtihdu.header['EQUINOX'] = (2000, '[yr] Equinox of celestial coord system')
  gtihdu.header['RADECSYS'] = ('FK5', 'Celestial coord system')
  gtihdu.header['RA_NOM'] = (hk['POINTING'][0], '[deg] R.A. of nominal aspect point [J2000]')  
  gtihdu.header['DEC_NOM'] = (hk['POINTING'][1], '[deg] Dec. of nominal aspect point [J2000]') 
  gtihdu.header['RA_OBJ'] = (hk['POINTING'][0], '[deg] Object Right ascension [J2000]')
  gtihdu.header['DEC_OBJ'] = (hk['POINTING'][1], '[deg] Object Declination [J2000]')
  gtihdu.header['TIMESYS'] = ('TT', 'Reference Time System')
  gtihdu.header['MJDREFI'] = (51544, '[d] MJD reference day (2000-01-01T00:00:00)')
  gtihdu.header['MJDREFF'] = (0.00074287037037037, '[d] MJD reference (fraction of day)')
  gtihdu.header['TIMEREF'] = ('LOCAL', 'Reference Frame')
  gtihdu.header['TASSIGN'] = ('SATELLITE', 'Time assigned by clock')
  gtihdu.header['TIMEUNIT'] = ('s', 'Time unit for timing header keyword')
  gtihdu.header['TIMEZERO'] = (0.0, '[s] Time Zero')
  gtihdu.header['TSTART'] = (min(hk['GTI_START']), '[s] Observation start time')
  gtihdu.header['TSTOP'] = (max(hk['GTI_STOP']), '[s] Observation stop time')
  gtihdu.header['DATE-OBS'] = (date_obs, 'Start date of observations')
  gtihdu.header['DATE-END'] = (date_end, 'End date of observations')
  gtihdu.header['CLOCKAPP'] = (True, 'Clock correction applied? (F/T)')
  gtihdu.header['ORIGIN'] = ('UNIVERSITY OF IOWA', 'Origin of fits file') 
  gtihdu.header['PROCVER'] = ('hsuf_20221026', 'Processing script version number') #prev. hsuf_20200131_hscl_20200131
  gtihdu.header['SOFTVER'] = ('Hea_11apr2022_V6.30.1', 'Software version')    
  gtihdu.header['CALDBVER'] = ('hs20200129', 'CALDB index version used')   
  gtihdu.header['TLM2FITS'] = ('db_20221002', 'Telemetry converter FITS version')   #location?!      
  gtihdu.header['CREATOR'] = (creator, 'Software creator of the file')
  gtihdu.header['DATE'] = (file_time, 'File creation date')  
  #checksum and datasum are added with fverify and fchecksum later  
  
  # join primary with extension and write
  hdulist = pyfits.HDUList([prihdu, evthdu, gtihdu])
  hdulist.writeto(outroot+'_s'+str(hk['GOOD_DPU'])+evt_type+'.evt', overwrite=True)
  
  return hk

#write a custom hk file with 8s light curve data for chosen LC bands, not the final version, RR
def write_hk(hk, outroot, main_dir, creator='', loud=0):  
    
  #create time strings
  date_obs = tai_time_conv([min(hk['GTI_START'])+32])[0]  #UTC time
  date_end = tai_time_conv([max(hk['GTI_STOP'])+32])[0]   #UTC time
  file_time = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
    
  # make the primary header (standard in all files)
  prihdr = pyfits.Header()
  prihdr['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
  prihdr['INSTRUME'] = ('SDD'+str(hk['GOOD_DPU']), 'Instrument name') 
  prihdr['OBS_ID'] = (hk['OBS_ID'], 'Observation ID')
  prihdr['OBJECT'] = (hk['OBJECT'], 'Object/Target Name')
  prihdr['OBJTYPE'] = (hk['OBJTYPE'], 'Object/Target type')     
  prihdr['DATE-OBS'] = (date_obs, 'Start date of observations')
  prihdr['DATE-END'] = (date_end, 'End date of observations')
  prihdr['DATE'] = (file_time, 'File creation date')  #checksum and datasum are added with fverify etc
  prihdu = pyfits.PrimaryHDU(header=prihdr) # make the primary header

  # define the columns for the FITS hk file, RR
  key_list=np.array(['TIME','SDD_TEMP','TEC_PWM','ADC_TEMP','LRS_PCNT','SDDHVSET','FLEXICNT',\
            'SDD1','MON_3P3V','MON_M5V','MON_P5V','SDD0','DAC_TEMP','BPL_TEMP',\
            'HSK_PCNT','SDDHVMON','DPU_TEMP','RA','DEC','SAT_LAT','SAT_LON',\
            'NADIR_ANGLE','OFFSET','LC_GAP1','LC_SCI','LC_ALL','LC_HARD','LC_VLE',\
            'LC_OXYGEN','LC_RESET','LC_ALSI','LC_UP'],dtype=str)
  key_unitlist=np.array(['TIME','SDDHVSET','SDD1','MON_3P3V','MON_M5V','MON_P5V','SDD0','SDDHVMON',\
                'RA','DEC','SAT_LAT','SAT_LON','NADIR_ANGLE','OFFSET','LC_GAP1','LC_SCI',\
                'LC_ALL','LC_HARD','LC_VLE','LC_OXYGEN','LC_RESET','LC_ALSI','LC_UP'],dtype=str)
  unit_list = np.array(['s','V','V','V','V','V','V','V','deg','deg','deg','deg','deg','deg',\
                'count/s','count/s','count/s','count/s','count/s','count/s','count/s','count/s','count/s'],dtype=str)
  hk_col=[]  
  for k in key_list: 
    kidx = np.where(key_unitlist==k)[0]  #check if data type has a unit, if so, attach in append call
    if len(kidx)==1: hk_col.append(pyfits.Column(name=k, format='1D', unit=unit_list[kidx][0], array=hk[k][hk['GOOD_IDX']]))
    else: hk_col.append(pyfits.Column(name=k, format='1D', array=hk[k][hk['GOOD_IDX']]))
  hk_col.append(pyfits.Column(name='IN_SAA', format='L', array=hk['IN_SAA'][hk['GOOD_IDX']]))    
  
  #define columns from prefilter for the FITS hk file, RR
  prefilter_list = np.array(['SAT_ALT','ELV','BR_EARTH','FOV_FLAG','SUNSHINE','SUN_ANGLE','MOON_ANGLE',\
    'RAM_ANGLE','ANG_DIST','COR_ASCA','COR_SAX','MCILWAIN_L','SAA','SAA_TIME','RAP','DECP',\
    'SAT_LATP','SAT_LONP'],dtype=str)
  key_unitlist = np.array(['SAT_ALT','ELV','BR_EARTH','SUN_ANGLE','MOON_ANGLE','RAM_ANGLE','ANG_DIST',\
                  'COR_ASCA','COR_SAX','SAA_TIME','RAP','DECP','SAT_LATP','SAT_LONP'],dtype=str)
  unit_list = np.array(['km','deg','deg','deg','deg','deg','deg','GeV/c','GeV/c','s','deg','deg','deg','deg'],dtype=str)
  for k in prefilter_list: 
      if (k=='FOV_FLAG') or (k=='SUNSHINE') or (k=='SAA'): hk_col.append(pyfits.Column(name=k, format='1I', array=hk[k][hk['GOOD_IDX']]))
      else: 
        kidx = np.where(key_unitlist==k)[0]  #check if data type has a unit, if so, attach in append call
        if len(kidx)==1: hk_col.append(pyfits.Column(name=k, format='1D', unit=unit_list[kidx][0], array=hk[k][hk['GOOD_IDX']]))
        else: hk_col.append(pyfits.Column(name=k, format='1D', array=hk[k][hk['GOOD_IDX']]))
  hk_col.append(pyfits.Column(name='POSITION', format='3E', unit='km', array=hk['POSITION'][hk['GOOD_IDX']]))   
  hk_col.append(pyfits.Column(name='VELOCITY', format='3E', unit='km/s', array=hk['VELOCITY'][hk['GOOD_IDX']]))   
  hk_col.append(pyfits.Column(name='ROLL', format='1D', unit='deg', array=hk['ROLL'][hk['GOOD_IDX']]))    
  hk_col.append(pyfits.Column(name='SUN_RA', format='1D', unit='deg', array=hk['SUN_RA'][hk['GOOD_IDX']]))
  hk_col.append(pyfits.Column(name='SUN_DEC', format='1D', unit='deg', array=hk['SUN_DEC'][hk['GOOD_IDX']]))  
    
  #put data in table
  cols = pyfits.ColDefs(hk_col) # create a ColDefs (column-definitions) object for all columns:
  hkhdu = pyfits.BinTableHDU.from_columns(cols) #, [], 'SPECTRUM') # create a new binary table HDU object     

  # add column keywords to header, RR
  hkhdu.header.comments['TTYPE1'] = 'Time'   #8th position
  hkhdu.header.comments['TFORM1'] = 'data format of field'
  hkhdu.header.comments['TUNIT1'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE2'] = '[deg C] Temperature readout for Si chip '
  hkhdu.header.comments['TFORM2'] = 'data format of field'
  hkhdu.header.comments['TTYPE3'] = 'Value duty cycle of the TEC'
  hkhdu.header.comments['TFORM3'] = 'data format of field'
  hkhdu.header.comments['TTYPE4'] = '[deg C] Temperature readout ADC'
  hkhdu.header.comments['TFORM4'] = 'data format of field'
  hkhdu.header.comments['TTYPE5'] = 'Low rate science packet'
  hkhdu.header.comments['TFORM5'] = 'data format of field'
  hkhdu.header.comments['TTYPE6'] = 'Voltage DAC on set line'
  hkhdu.header.comments['TFORM6'] = 'data format of field'
  hkhdu.header.comments['TUNIT6'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE7'] = 'FLEXI counter'
  hkhdu.header.comments['TFORM7'] = 'data format of field'
  hkhdu.header.comments['TTYPE8'] = 'reset pulse threshold '
  hkhdu.header.comments['TFORM8'] = 'data format of field'
  hkhdu.header.comments['TUNIT8'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE9'] = 'Voltage Monitor on 3.3V line'
  hkhdu.header.comments['TFORM9'] = 'data format of field'
  hkhdu.header.comments['TUNIT9'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE10'] = 'Voltage Monitor on -5V line'
  hkhdu.header.comments['TFORM10'] = 'data format of field'
  hkhdu.header.comments['TUNIT10'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE11'] = 'Voltage Monitor on +5V line'
  hkhdu.header.comments['TFORM11'] = 'data format of field'
  hkhdu.header.comments['TUNIT11'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE12'] = 'event detection threshold'
  hkhdu.header.comments['TFORM12'] = 'data format of field'
  hkhdu.header.comments['TUNIT12'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE13'] = '[deg C] Temperature of the DAC'
  hkhdu.header.comments['TFORM13'] = 'data format of field'
  hkhdu.header.comments['TTYPE14'] = '[deg C] Temperature of the baseplate sensor'
  hkhdu.header.comments['TFORM14'] = 'data format of field'
  hkhdu.header.comments['TTYPE15'] = 'Housekeeping packet counter'
  hkhdu.header.comments['TFORM15'] = 'data format of field'
  hkhdu.header.comments['TTYPE16'] = 'Voltage on DAC monitor line'
  hkhdu.header.comments['TFORM16'] = 'data format of field'
  hkhdu.header.comments['TUNIT16'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE17'] = '[deg C] Temperature of DPU'
  hkhdu.header.comments['TFORM17'] = 'data format of field'
  
  #bus data
  hkhdu.header.comments['TTYPE18'] = 'Right Ascension of the pointing'  #50th header
  hkhdu.header.comments['TFORM18'] = 'data format of field'
  hkhdu.header.comments['TUNIT18'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE19'] = 'Declination of the pointing'
  hkhdu.header.comments['TFORM19'] = 'data format of field'
  hkhdu.header.comments['TUNIT19'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE20'] = 'Satellite Latitude'
  hkhdu.header.comments['TFORM20'] = 'data format of field'
  hkhdu.header.comments['TUNIT20'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE21'] = 'Satellite Longitude'
  hkhdu.header.comments['TFORM21'] = 'data format of field'
  hkhdu.header.comments['TUNIT21'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE22'] = 'Angular distance from pointing to nadir'
  hkhdu.header.comments['TFORM22'] = 'data format of field'
  hkhdu.header.comments['TUNIT22'] = 'physical unit of field'  
  hkhdu.header.comments['TTYPE23'] = 'Offset'
  hkhdu.header.comments['TFORM23'] = 'data format of field'
  hkhdu.header.comments['TUNIT23'] = 'physical unit of field'
  
  #count rate data
  hkhdu.header.comments['TTYPE24'] = 'Averaged 8s Count Rate in LC_GAP1 band'  #68th header
  hkhdu.header.comments['TFORM24'] = 'data format of field'
  hkhdu.header.comments['TUNIT24'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE25'] = 'Averaged 8s Count Rate in LC_SCI band'
  hkhdu.header.comments['TFORM25'] = 'data format of field'
  hkhdu.header.comments['TUNIT25'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE26'] = 'Averaged 8s Total Count Rate'
  hkhdu.header.comments['TFORM26'] = 'data format of field'
  hkhdu.header.comments['TUNIT26'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE27'] = 'Averaged 8s Count Rate in LC_HARD band'
  hkhdu.header.comments['TFORM27'] = 'data format of field'
  hkhdu.header.comments['TUNIT27'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE28'] = 'Averaged 8s Count Rate in LC_VLE band'
  hkhdu.header.comments['TFORM28'] = 'data format of field'
  hkhdu.header.comments['TUNIT28'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE29'] = 'Averaged 8s Count Rate in LC_OXYGEN band'
  hkhdu.header.comments['TFORM29'] = 'data format of field'
  hkhdu.header.comments['TUNIT29'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE30'] = 'Averaged 8s Count Rate in LC_RESET band'
  hkhdu.header.comments['TFORM30'] = 'data format of field'
  hkhdu.header.comments['TUNIT30'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE31'] = 'Averaged 8s Count Rate in LC_ALSI band'
  hkhdu.header.comments['TFORM31'] = 'data format of field'
  hkhdu.header.comments['TUNIT31'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE32'] = 'Averaged 8s Count Rate in LC_UP band'
  hkhdu.header.comments['TFORM32'] = 'data format of field'
  hkhdu.header.comments['TUNIT32'] = 'physical unit of field'

  #add columns from prefilter
  hkhdu.header.comments['TTYPE33'] = 'Satellite in SAA?'  #95th header
  hkhdu.header.comments['TFORM33'] = 'data format of field'
  hkhdu.header.comments['TTYPE34'] = 'distance from earth center'
  hkhdu.header.comments['TFORM34'] = 'data format of field'
  hkhdu.header.comments['TUNIT34'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE35'] = 'angle between pointing and earth limb'  #100th header
  hkhdu.header.comments['TFORM35'] = 'data format of field'
  hkhdu.header.comments['TUNIT35'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE36'] = 'angle between pointing and bright earth'
  hkhdu.header.comments['TFORM36'] = 'data format of field'
  hkhdu.header.comments['TUNIT36'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE37'] = '0=sky; 1=dark earth; 2=bright earth'
  hkhdu.header.comments['TFORM37'] = 'data format of field'
  hkhdu.header.insert(108,('TNULL37',-999, 'null value'))   
  hkhdu.header.comments['TTYPE38'] = '1=in sunshine; 0=not'
  hkhdu.header.comments['TFORM38'] = 'data format of field'
  hkhdu.header.insert(111,('TNULL38',-999, 'null value'))    
  hkhdu.header.comments['TTYPE39'] = 'angle between pointing vector and sun vector'
  hkhdu.header.comments['TFORM39'] = 'data format of field'
  hkhdu.header.comments['TUNIT39'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE40'] = 'angle between pointing vector and moon vector'
  hkhdu.header.comments['TFORM40'] = 'data format of field'
  hkhdu.header.comments['TUNIT40'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE41'] = 'angle between pointing and velocity vectors'
  hkhdu.header.comments['TFORM41'] = 'data format of field'
  hkhdu.header.comments['TUNIT41'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE42'] = 'angular distance of pointing from nominal'
  hkhdu.header.comments['TFORM42'] = 'data format of field'
  hkhdu.header.comments['TUNIT42'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE43'] = 'magnetic cut off rigidity (ASCA map)'
  hkhdu.header.comments['TFORM43'] = 'data format of field'
  hkhdu.header.comments['TUNIT43'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE44'] = 'magnetic cut off rigidity (IGRF map)'
  hkhdu.header.comments['TFORM44'] = 'data format of field'
  hkhdu.header.comments['TUNIT44'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE45'] = 'McIlwain L parameter (SAX)'
  hkhdu.header.comments['TFORM45'] = 'data format of field'
  hkhdu.header.comments['TTYPE46'] = '1=in SAA; 0=not'
  hkhdu.header.comments['TFORM46'] = 'data format of field'
  hkhdu.header.insert(134,('TNULL46',-999, 'null value'))    
  hkhdu.header.comments['TTYPE47'] = 'time since entering/exiting SAA'
  hkhdu.header.comments['TFORM47'] = 'data format of field'
  hkhdu.header.comments['TUNIT47'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE48'] = 'Right Ascension by prefilter'
  hkhdu.header.comments['TFORM48'] = 'data format of field'
  hkhdu.header.comments['TUNIT48'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE49'] = 'Declination by prefilter'
  hkhdu.header.comments['TFORM49'] = 'data format of field'
  hkhdu.header.comments['TUNIT49'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE50'] = 'Satellite Latitude by prefilter'
  hkhdu.header.comments['TFORM50'] = 'data format of field'
  hkhdu.header.comments['TUNIT50'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE51'] = 'Satellite Longitude by prefilter'
  hkhdu.header.comments['TFORM51'] = 'data format of field'
  hkhdu.header.comments['TUNIT51'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE52'] = 'ECI Position satellite (X, Y, Z)'
  hkhdu.header.comments['TFORM52'] = 'data format of field'
  hkhdu.header.comments['TUNIT52'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE53'] = 'ECI velocity satellite (vX,vY,vZ)'
  hkhdu.header.comments['TFORM53'] = 'data format of field'
  hkhdu.header.comments['TUNIT53'] = 'physical unit of field'
  hkhdu.header.comments['TTYPE54'] = 'pointing axis roll'
  hkhdu.header.comments['TFORM54'] = 'data format of field'
  hkhdu.header.comments['TUNIT54'] = 'physical unit of field'   
  hkhdu.header.comments['TTYPE55'] = 'Right Ascension of Sun'
  hkhdu.header.comments['TFORM55'] = 'data format of field'
  hkhdu.header.comments['TUNIT55'] = 'physical unit of field'   
  hkhdu.header.comments['TTYPE56'] = 'Declination of Sun'
  hkhdu.header.comments['TFORM56'] = 'data format of field'
  hkhdu.header.comments['TUNIT56'] = 'physical unit of field' 

  #set standard extension header keys, updated Dec12 2019 RR
  hkhdu.header['EXTNAME'] = ('HK', 'Binary table extension name')  
  hkhdu.header['HDUCLASS'] = ('OGIP', 'format conforms to OGIP/GSFC standards')
  hkhdu.header['HDUCLAS1'] = ('TEMPORALDATA', 'First class level')  
  hkhdu.header['HDUCLAS2'] = ('HK', 'Second class level')
  hkhdu.header['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
  hkhdu.header['INSTRUME'] = ('SDD'+str(hk['GOOD_DPU']), 'Instrument name')   
  hkhdu.header['OBSERVER'] = ('PHILIP KAARET', 'Principal Investigator')
  hkhdu.header['OBS_ID'] = (hk['OBS_ID'], 'Observation ID')
  hkhdu.header['OBJECT'] = (hk['OBJECT'], 'Object/Target Name')
  hkhdu.header['OBJTYPE'] = (hk['OBJTYPE'], 'Object/Target type')   
  #if (float(hk['OBS_ID']) >=35501) & (float(hk['OBS_ID'])<37301):
  #  hkhdu.header['COMMENT'] = (hk['COMMENT'], 'Offset target association')  
  hkhdu.header['EQUINOX'] = (2000, '[yr] Equinox of celestial coord system')
  hkhdu.header['RADECSYS'] = ('FK5', 'Celestial coord system')
  hkhdu.header['RA_NOM'] = (hk['POINTING'][0], '[deg] R.A. of nominal aspect point [J2000]')  
  hkhdu.header['DEC_NOM'] = (hk['POINTING'][1], '[deg] Dec. of nominal aspect point [J2000]')  
  hkhdu.header['RA_OBJ'] = (hk['POINTING'][0], '[deg] Object Right ascension [J2000]')
  hkhdu.header['DEC_OBJ'] = (hk['POINTING'][1], '[deg] Object Declination [J2000]')
  hkhdu.header['TIMESYS'] = ('TT', 'Reference Time System')
  hkhdu.header['MJDREFI'] = (51544, '[d] MJD reference day (2000-01-01T00:00:00)')
  hkhdu.header['MJDREFF'] = (0.00074287037037037, '[d] MJD reference (fraction of day)')
  hkhdu.header['TIMEREF'] = ('LOCAL', 'Reference Frame')
  hkhdu.header['TASSIGN'] = ('SATELLITE', 'Time assigned by clock')
  hkhdu.header['TIMEUNIT'] = ('s', 'Time unit for timing header keyword')
  hkhdu.header['TIMEDEL'] = (8.0,'[s] Data time resolution.')   
  hkhdu.header['TIMEZERO'] = (0.0, '[s] Time Zero')
  hkhdu.header['TIMEPIXR'] = (1, 'Bin time beginning=1 middle=0.5 end=1')
  hkhdu.header['TIERRELA'] = (4.0E-06, '[s/s] relative errors expressed as rate')
  hkhdu.header['TIERABSO'] = (1.0, '[s] timing precision in seconds')
  hkhdu.header['TSTART'] = (min(hk['GTI_START']), '[s] Observation start time')
  hkhdu.header['TSTOP'] = (max(hk['GTI_STOP']), '[s] Observation stop time')
  hkhdu.header['DATE-OBS'] = (date_obs, 'Start date of observations')
  hkhdu.header['DATE-END'] = (date_end, 'End date of observations') 
  hkhdu.header['CLOCKAPP'] = (True, 'Clock correction applied? (F/T)')
  hkhdu.header['ORIGIN'] = ('UNIVERSITY OF IOWA', 'Origin of fits file')
  hkhdu.header['PROCVER'] = ('hsuf_20221026', 'Processing script version number')   #prev. hsuf_20200131_hscl_20200131
  hkhdu.header['SOFTVER'] = ('Hea_11apr2022_V6.30.1', 'Software version')   
  hkhdu.header['CALDBVER'] = ('hs20200129', 'CALDB index version used')   
  hkhdu.header['TLM2FITS'] = ('db_20221002', 'Telemetry converter FITS number')   #location?!  
  hkhdu.header['CREATOR'] = (creator, 'Software creator of the file')
  hkhdu.header['DATE'] = (file_time, 'File creation date')  
  #checksum and datasum are added with fverify and fchecksum later
    
  # join primary with extension, write, and check
  hdulist = pyfits.HDUList([prihdu, hkhdu])
  hdulist.writeto(outroot+'_s'+str(hk['GOOD_DPU'])+'.hk', overwrite=True)

  return
