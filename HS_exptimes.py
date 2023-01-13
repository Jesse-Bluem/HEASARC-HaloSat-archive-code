'''
HS_exptimes.py
Feb 28, 2020
Rebecca Ringuette

- parse through all directories under main dir given, default is AllData on halosat2 server
- collect exposure times from .pi files for each detector
- keep total of three detectors
- put into a fits file for planning use

Updates 10/9/2022 JKB
Converted to Python 3 and removed unused OS import.
Changed a dict.keys() usage with a loop that was broken in Py3.
Removed iteritems usage which has been removed in Python 3.
NOTE: running subset of fields resulted in data for HS1031 being written to fits under 133, which was not run in the subset. All other fields were fine.
NOTE 2: Running additional fields moved the data for 1031 to a different field (138). Other fields after 134 also broke and wrote to the wrong fields, starting with 319.
Fix on line 52: change trigger condition for sequence jump to dark earth to be 400 instead of 100, as 400 > than # of regular fields.
'''

import astropy.io.fits as pyfits
import sys, datetime, glob, operator
import numpy as np

#initialize top directory and exp_time
if len(sys.argv)==2: main_dir=sys.argv[1]
else: main_dir='/home/halosat/AllData/v1_local/'
code_dir = '/home/halosat/Server_Code/'
out_dir = '/home/halosat/halosat_planning/'
exp_times, end_times = {}, {}
todays_date=str(datetime.date.today())

#get list of pi files with paths, initialize data storage
pi_files = glob.glob('/home/halosat/AllData/v1_local/*/products/*.pi')
for pi in pi_files:
    t = pyfits.open(pi)
    hsID = t['SPECTRUM'].header['OBS_ID']
    if hsID[0:4] not in exp_times.keys(): 
        exp_times[hsID[0:4]] = [t['SPECTRUM'].header['EXPOSURE']]
        end_times[hsID[0:4]] = [t['SPECTRUM'].header['DATE-END']]
    else: 
        exp_times[hsID[0:4]].append(t['SPECTRUM'].header['EXPOSURE']) 
        end_times[hsID[0:4]].append(t['SPECTRUM'].header['DATE-END'])
   
#insert zero exposure time for missing key values
keyser=[]#orig py2 method doesnt work for py3, just iterate dict to list JKB oct 2022
#note: This originally used dict.keys(), which still works in Py3 and all other instances in the code except for this one.
#keys() format was changed, which seemed to only break this instance, where it looped over .keys, in .keys seems to be unaffected.
for key in exp_times:
	keyser.append(key)
test = np.sort(keyser)
diff = np.diff(np.array(test,dtype=int))
idx = np.where((diff>1) & (diff<400))[0]  #ignore value jump for dark earth fields - JKB 10/2022 THIS BREAKS RUNNING THE CODE FOR A SUBSECTION OF FIELDS IF THE GAP IS 100 OR MORE. CHANGED TO 400
for i in idx:
    for k in range(int(diff[i]-1)): 
        exp_times[str('{:04d}'.format(int(test[i])+1+k))]=[0.]    
        end_times[str('{:04d}'.format(int(test[i])+1+k))]=['2000-01-01T00:00:00']
try: #Until Cusp 1 data comes in.
    print(exp_times['0374'])
except:
    exp_times['0374']=[0.]
    end_times['0374']=['2000-01-01T00:00:00']
try: #Until Cusp 2 data comes in.
    print(exp_times['0375'])
except:
    exp_times['0375']=[0.]
    end_times['0375']=['2000-01-01T00:00:00']

keyser=[]
for key in exp_times:
	keyser.append(key)
#collect hsIDs and total exposure times into arrays
hsIDs = np.sort(keyser)
total_time = np.zeros(len(hsIDs),dtype=float)
for i in range(len(hsIDs)): 
    total_time[i]=sum(exp_times[hsIDs[i]])
#end_time = max(max(end_times.iteritems(), key=operator.itemgetter(1))[1]) #defunct with Python 3, iteritems has been removed.
ender=[]
for key in end_times:
	ender.append(end_times[key])
end_time=max(max(ender)) #JKB oct 2022 is this equivalent? just want max GTI value. This replaced iteritems. After testing output it looks good.

#read in fits targets file
t = pyfits.open(code_dir+'haloSat_targets.fits')
HSIDs = t[1].data['Target_ID']
d = t[1].data
h0 = t[0].header
h1 = t[1].header
t.close()

#extend hsIDs with zeros for non-science source to compare
total_time = np.append(total_time,np.repeat(0.,len(HSIDs)-len(hsIDs)))
hsIDs = np.append(hsIDs,np.repeat('0000',len(HSIDs)-len(hsIDs)))    

#compare hsIDs, print if different
diff = abs(np.array(hsIDs,dtype=int)-np.array(HSIDs,dtype=int))
idx = np.where(diff>0)[0]
if len(idx)>0: 
    idx2 = np.where(hsIDs[idx]!='0000')[0]  #ignore sightline and avoid sources
    if len(idx2)>0:
      print('ID number order does not match!')
      print(hsIDs[idx[idx2]])
      #shouldn't output anything if it worked correctly

#add column to fits object
cols=t[1].columns
c1=pyfits.Column(name='Observed_Seconds',array=total_time,unit='s',format='E')
new_cols=d.columns+c1
hdu=pyfits.BinTableHDU.from_columns(new_cols)
hdu.header['FILEDATE']=(todays_date,'Date file was generated.')
hdu.header['ENDTIME']=(end_time,'End of latest GTI.')
hdu.header['EXTNAME']=('TARGETS','Extension name.')
prihdu = pyfits.PrimaryHDU(header=h0)
finalhdu=pyfits.HDUList([prihdu,hdu])
finalhdu.writeto(out_dir+'haloSat_Targets.fits',overwrite=True)
