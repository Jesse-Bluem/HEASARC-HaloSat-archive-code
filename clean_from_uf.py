#UPDATES 10/9/2022 JKB - removed unused imports. Fixed a variable typo (probably never caught because that print statement was only for error conditiions). Updated to Python 3. Added a check during PDF file write to make sure GTI survived cuts before checking for sun angle, as a GTI of zero length after cuts crashed the code. It now returns an automatic zero sun angle if the GTI length was reduced to zero, rather than querying and breaking.
#JKB 12/7 - PDF y-axis scaling searched for a zero, which returned a 0 value for fields with cuts >= 1. I changed to look for an = sign instead, so any number works. Added similar axis scaling to VLE, which was fixed at 0 to 2. Moved tight layout on PDF page 1 to before text writing code block, as it was not applying tight layout in Python 3 due to the text, for some reason.

import os, sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numpy import array, where, transpose
from numpy import zeros, linspace, concatenate, mean, percentile, vstack
from datetime import datetime, timedelta
import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord


def equatorial2galactic(ra, dec):
# Converts a list of equatorial coordinates (ra,dec) to galactic (l,b).
  # Insert equatorial coordinates into Astropy SkyCoord
  equ = SkyCoord(ra=ra, dec=dec, unit='degree', frame='icrs')
  # Get the galactic coordinates
  gal = equ.galactic
  l = gal.l.value
  b = gal.b.value
  # Return (l,b)
  return (l, b)

def getDate(time) :
  epoch = datetime(2000, 1, 1, 0, 0, 0) # reference epoch yr = 2000.0 for BCT EDU
  # epoch = utc.localize(epoch) # epoch localized to UTC
  time1 = time - 37.0
  date = epoch + timedelta(seconds=time1)
  return date

def find_obs(filename) :
  if os.path.isfile(filename) :
    h = pyfits.open(filename) # read file
    # header = h['GTI'].header # get header of GTI extension
    d = h['GTI'].data # get data in GTI extension
    tstart = d.field('START') 
    tstop = d.field('STOP')
    h.close() # close the file
    # find contiguous sets of GTIs and group them into observations
    obs_start = [tstart[0]]
    obs_end = [tstop[0]]
    obs_total = [tstop[0]-tstart[0]]
    obs_gap = 86400.0 # maximum gap between pointings in an observation
    for s, e in transpose([tstart[1:], tstop[1:]]) :
      if (s - obs_end[-1]) < obs_gap : # if new GTI is soon enough after end of current observation
        obs_end[-1] = e # move the observation end to the end of the GTI
        obs_total[-1] += e-s # add GTI to total time
      else :
        obs_start.append(s) # otherwise start a new observation
        obs_end.append(e)
        obs_total.append(e-s)
    # return arrays with start and end of each observation
    return vstack((array(obs_start), array(obs_end), array(obs_total)))
  else :
    print('Cannot find file ', filename)
    return (0.0, 0.0, 0.0)


# directory containing HaloSat files such as response, target list, Python modules
# is set in environment variable HALOSOFT
if 'HALOSOFT' in os.environ.keys() :
  hs_dir = os.environ['HALOSOFT']
else : # complain if not set
  print('Set the environment variable HALOSOFT to point at the HaloSat software directory.')
  print("In bash, use, e.g., export HALOSOFT=/home/halosat/hscore/")
  sys.exit(1)
sys.path.append(hs_dir) # append directory to Python path
import hs_clean_p3 as hs

# default values for various parameters
# default is to do no timecuts, by not setting timecuts variable
group = False # default is to do no grouping
gti_min = 64.0 # [seconds] minimum GTI duration
gti_max = 64.0 # [seconds] maximum GTI duration 
tempcor = True # do temperature-dependent gain correction
pdf = True # write all plots and text output to a pdf

if len(sys.argv) == 4 : # have command line arguments
  # plt.ioff() # turn off interactive plotting
  inroot = sys.argv[1] # assume first is input root
  outroot = sys.argv[2] # and second is output root
  filt = sys.argv[3] # filter expression
  # redirect output to outroot+'.log'
  orig_stdout = sys.stdout
  f = open(outroot+'.log', 'w')
  sys.stdout = f
elif len(sys.argv) >= 5 : # have command line arguments
  # plt.ioff() # turn off interactive plotting
  inroot = sys.argv[1] # assume first is input root
  outroot = sys.argv[2] # and second is output root
  filt = sys.argv[3] # filter expression
  # interpret all remaining arguments as Python expressions
  for i in range(4, len(sys.argv)) :
    exec(sys.argv[i]) # evaluate argument as Python expression
else :
  print('Need at least three command line arguments:')
  print('1) input root filename with directory')
  print('2) output root filename with directory')
  print('3) filter expression')
  print('Extra arguments are interpreted as Python expressions.')
  print('Use to set variables: timecuts, group, gti_min, gti_max, tempcor, pdf')

# redirect output to outroot+'.log'
orig_stdout = sys.stdout
f = open(outroot+'.log', 'w')
sys.stdout = f

def cut_value(filt, cut) :
# extract a cut level from the filter expression
# assume the cut level starts with '0' and ends with ')'
# return 0 if cut if not found
  i = filt.find(cut)
  print(filt)
  if i >= 0 :
    j = filt[i:-1].find('=') #JKB changed from searching for an assumed 0 start to the equals sign, since some cuts don't start with a zero.
    k = filt[i+j:].find(')')
    #print(filt, cut, filt[i+j+1:i+j+k], i, j+1, k)
    return float(filt[i+j+1:i+j+k])
  else :
    return 0.0

# attempt to extract cut levels from filter expression
hard_cut = cut_value(filt, 'LC_HARD')
vle_cut = cut_value(filt, 'LC_VLE')
#sun_cut = cut_value(filt, 'SUN_ANGLE')
#print 'hard_cut = ', hard_cut, '  vle_cut = ', vle_cut
#, ' sun_cut = ', sun_cut

if pdf :
  # the plot output file
  pdff = PdfPages(outroot+'.pdf')
  print('Will save output to ', outroot+'.pdf')
else : plt.ion() # interactive plotting

print('Input file root =', inroot)
print('Output file root =', outroot)
if 'timecuts' in locals() : print('Time cuts =', timecuts) #no timecuts?
else : print('No time cuts')
print('Time bin size, min =', gti_min, ' max =', gti_max)
print('Temperature correction =', tempcor)
print('Group =', group)
print('Filter expression =', filt)
print

# find pointing from attitude file
# pointing = hs.find_pointing(inroot)
## generate filter files using prefilter
#print 'Generating prefilter files'
#hs.run_prefilter(inroot, pointing, hs_dir=hs_dir)
#print


pcolor = {14:'r', 54:'g', 38:'b'} # plot color dictionary based on DPU number

spec = {} # dictionary to hold spectra
ontime = {} # dictionary to hold on times
exposure = {} # dictionary to hold exposure times after cuts
events = {} # dictionary to hold number of events in SCI band after cuts
hkp = {} # merged HK dictionary for plotting

firstdpu = True
for dpu in [14, 54, 38] :
  print('Processing data for DPU', dpu)
  # read in unfiltered data
  if 'timecuts' in locals() : # use time cuts if specified
    hk, evt = hs.read_uf(inroot, dpu, timecuts=timecuts, hs_dir=hs_dir, loud=1) #no timecuts?
  else :
    hk, evt = hs.read_uf(inroot, dpu, hs_dir=hs_dir, loud=1) 

  # print some statistics
  print('Total exposure (s) = {0:.1f}, science band events = {1:.3f}'\
        .format(8*float(sum(hk['IN_SAA'] == False)), 8*float(sum(hk['LC_SCI']))))
  ontime[dpu] = 8*float(sum(hk['IN_SAA'] == False))

  # rebin the HK data with larger time bins
  hk2 = hs.time_rebin(hk, gti_min=gti_min, gti_max=gti_max, loud=1)
  # add dictionary entry = VLE rate * rigidity
  hk2['VLE_COR'] = hk2['LC_VLE']*hk2['COR_SAX']

  # filter on new time bins
  print('Filtering with ', filt)
  # filter expression
  g = eval(filt)
  hk2['GOOD'] = g # save which are good intervals in hk2
  if (sum(g) == 0) :
    print('*** No good time intervals found.  Please check your cuts or your data files.')
    if pdf :
      pdff.close() # write the pdf
      os.remove(outroot+'.pdf') # delete the file, since it is probably empty
      sys.stdout = orig_stdout # set stdout back the way it started
      f.close() # close the output file
    sys.exit(1) # exit with error code 1 

  gti0 = hk2['TIME'][g]
  gti1 = hk2['TIME_END'][g]
  # calculate some statistics
  gtime = sum(gti1 - gti0) # total time in GTI

  # total events in first LC band in GTI
  gevents = sum(hk2['LC_SCI'][g]*(gti1-gti0))
  print('After filtering, exposure (s) = {0:.1f}, science band events = {1:.1f}'\
        .format(gtime, gevents))

  # find continuous sets of GTI and merge them
  gti = [[gti0[0], gti1[0]]] # list for merged GTIs, start with first unmerged GTI
  gti_gap = 0.3 # maximum gap between GTI to be considered continuous
  for s, e in transpose([gti0[1:], gti1[1:]]) :
    if (s - gti[-1][1]) < gti_gap : # if new GTI is soon enough after the current GTI end
      gti[-1][1] = e # move the GTI end
    else :
      gti.append([s, e]) # otherwise start a new observation
  # change GTI to have list of start times, then list of end times
  gti = transpose(gti) 
  print('Total time in GTIs = ', sum(gti[1]-gti[0]))

  # write files with the spectrum and cleaned events
  # filter string in Lorella's format
  filt_nohk = filt.replace('hk2[', '').replace(']', '')
  print('Filter = '+filt_nohk)
  try : groupi = int(group) # convert group string to an integer
  except ValueError : groupi = 0 # if fail, then set to 0 = False
  (s, t) = hs.write_evt_pi(evt, dpu, gti, inroot, outroot, filt_nohk, \
                           group=groupi, loud=1, ratekey=True)
  spec[dpu] = s
  exposure[dpu] = t
  events[dpu] = gevents
  print

  # merge the hk2 dictionaries into hkp for plotting
  if firstdpu :
    for k, v in hk2.items(): # k = key, v = value of key
      hkp[k] = v
    firstdpu = False
  else : 
    for k, v in hk2.items(): # k = key, v = value of key
      hkp[k] = concatenate((hkp[k], v))

# make plots
if pdf :
  # create a figure that is 8 by 10.5 inches
  fig1 = plt.figure(outroot, figsize=(10.0, 7.5), dpi=100)
  plt.clf()
  #plt.title(inroot+' '+filt)
  grid_size = (3, 3)

# generate a list of strings with observation info
# find observation start and stop dates and exposure times for
# attitude, unfiltered, and clean files 
# find observation times from attitude file
t = find_obs(inroot+'.att')
atts, atte, atttot = t[0], t[1], t[2]
ftimes = zeros((len(atttot), 6)) # make array to hold times after filtering

# find observation times from unfiltered files
ot = [find_obs(inroot+'_s14_uf.evt')]
ot.append(find_obs(inroot+'_s54_uf.evt'))
ot.append(find_obs(inroot+'_s38_uf.evt'))
# find observation times from clean files
ot.append(find_obs(outroot+'_s14_cl.evt'))
ot.append(find_obs(outroot+'_s54_cl.evt'))
ot.append(find_obs(outroot+'_s38_cl.evt'))

# match filtered observations with attitude observations
for i in range(len(atts)) :
  for j in range(6) : 
    q = (atts[i] <= ot[j][0]) & (ot[j][0] <= atte[i]) # find matching observation
    if sum(q) == 1 : ftimes[i][j] = ot[j][2][q]

#print ftimes
# get the nominal pointing coordinates from the RA_NOM and DEC_NOM keywords
# in the attitude file
h = pyfits.open(inroot+'.att') # read attitude file
attheader = h['GTI'].header # get header of GTI extension
h.close() # close the file
ra_nom = attheader['RA_NOM'] # get the coordinates
dec_nom = attheader['DEC_NOM']
obs_id = attheader['OBS_ID'] # get the Observation ID
objname = attheader['OBJECT'] # get the target name
objtype = attheader['OBJTYPE'] # get the target name

# get the old HaloSat name from the targets file
# in the attitude file
h = pyfits.open(hs_dir+'haloSat_targets.fits')
d = h[1].data # get data in extension 1
h.close() # close the file
target = int(obs_id[0:4]) # get target ID from Obs ID
target_id = d['Target_ID'] # get the target IDs
hsnames = d['Target_Name'] # get the old target names
q = where(target == target_id)[0]
if len(q) == 0 :
  print('Target ID ', target_id, ' not found in ', hs_dir+'haloSat_targets.fits') #JKB, fix capitalization error on Target_ID 
  sys.exit(1)
hsname = hsnames[q][0]
#print 'Old HaloSat name = ', hsname

# list of strings describing observation
if inroot.find('/') != -1 : # there is a / in the inroot name, so strip out the directory
  inroot_nodir = inroot[inroot.find('/')+1:]
else :
  inroot_nodir = inroot
s = ['ObsID='+obs_id+'  Name='+objname]
# find Galactic coordinates
(l_nom, b_nom) = equatorial2galactic(ra_nom, dec_nom)
s.append('RA={0:7.3f}'.format(ra_nom)+' DEC={0:+7.3f}'.format(dec_nom)+ \
         '  l={0:7.3f}'.format(l_nom)+' b={0:+7.3f}'.format(b_nom))
s.append('AltName='+hsname+'  Type='+objtype)
filt_nohk = filt.replace('hk2[', '').replace(']', '')
s.append('Filter = '+filt_nohk)
#t = 'Science band events:'
#for dpu in [14, 54, 38] :
#  t += ' D{0:d}={1:d}'.format(dpu, int(events[dpu]))
#s.append(t)
for dpu in [14, 54, 38] :
  s.append('S{0:d} unfilt={1:.1f}ks clean={2:.1f}ks events={3:d}'.format(dpu, \
           ontime[dpu]/1E3, exposure[dpu]/1E3, int(events[dpu])))
# find range of baseplate temperature
t = hkp['BPL_TEMP']
s.append('BPL Temp = {0:.1f}C ({1:.1f}-{2:.1f})'.\
           format(mean(t), percentile(t, 5), percentile(t, 95)))

if pdf : plt.subplot2grid(grid_size, (0, 2), rowspan=1, colspan=1)
# plot the hard rate versus Sun angle after filtering
hs.plotxy(hkp, 'SUN_ANGLE', 'LC_HARD', title='Hard rate vs Sun angle',
          onlygood=True, color=pcolor, pdf=pdf,
          xbin=linspace(min(hk2['SUN_ANGLE'][g]), max(hk2['SUN_ANGLE'][g]), 16),
          dt=64.0, minbin=6)
## plot the angle between the target and the Sun
##hs.plotts(hkp, 'SUN_ANGLE', title='Sun angle', pdf=pdf, psym='.', color=pcolor)

if pdf : plt.subplot2grid(grid_size, (1, 0), rowspan=1, colspan=2)
# plot the rate of hard band events versus time bin
hs.plotts(hkp, 'LC_HARD', title='Hard rate', psym='.', color=pcolor, \
	      collapse=True, dt=64.0, ylim=[0.0, max(0.5, 2*hard_cut)], 
        plot_bar=hard_cut, pdf=pdf)
          #, hist=True, nbin=32)
# print hard_cut, [0.0, max(0.5, 2*hard_cut)]

if pdf : plt.subplot2grid(grid_size, (1, 2), rowspan=1, colspan=1)
# plot the VLE rate versus rigidity
hs.plotxy(hkp, 'PFLUX', 'LC_HARD', title='Hard rate vs proton flux', pdf=pdf,
        ylim=[0.0, max(0.5, 2*hard_cut)], psym='.', color=pcolor)

if pdf : plt.subplot2grid(grid_size, (2, 0), rowspan=1, colspan=2)
# plot the rate of VLE band events versus time bin. JKB 12/7 edit to scale VLE plot Y axis like hard rate axis
hs.plotts(hkp, 'LC_VLE', title='VLE rate', psym='.', color=pcolor, legend='best', \
	      collapse=True, dt=64.0, ylim=[0.0, max(2.0, 2*vle_cut)], plot_bar=vle_cut, pdf=pdf)
          # hist=True, nbin=32,

if pdf : plt.subplot2grid(grid_size, (2, 2), rowspan=1, colspan=1)
# plot the proton flux versus time bin
hs.plotts(hkp, 'PFLUX', title='Proton Flux', pdf=pdf, psym='.', color=pcolor)

# plot the VLE rate versus rigidity
#hs.plotxy(hkp, 'COR_SAX', 'LC_VLE', title='VLE rate vs rigidity', pdf=pdf,
#	      ylim=[0.0, vle_cut], psym='.', color=pcolor)
## array that runs over range of rigidity
#cor = linspace(min(hkp['COR_SAX']), max(hkp['COR_SAX']), 100)
#hs.plotts(hkp, 'ANTISUN_DLON', title='ELon-ELon(AntiSun)', pdf=pdf, psym='.', color=pcolor)

if pdf :   
  plt.tight_layout() #JKB 12/7 - apply tight layout before adding text in Python 3 or it messes up PDF for some reason.
  hp = 0.05
  monofont = {'fontname':'DejaVu Sans Mono'}
  for i in range(len(s)) :
    plt.text(hp, 0.96-0.03*i, s[i], transform=fig1.transFigure, **monofont)
else :
  print(s)

if pdf :    
  #plt.tight_layout() JKB 12/7 moved this up, was erroring out in Py3 and not applying tight layout. Text block was issue.
  pdff.savefig(fig1) # save the figure to the pdf file
  # create a figure that is 8 by 10.5 inches
  fig2 = plt.figure(outroot, figsize=(10.0, 7.5), dpi=100)
  plt.clf()
  #plt.title(inroot+' '+filt)
  grid_size = (3, 3)

#  # plot the rate of reset events
#  hs.plotts(hk2, 'LC_RESET', title='Reset rate', 
#            collapse=True, orb_break=False,
#            add=True, dpu=dpu, psym='+'+pcolor[dpu], legend='best')

if pdf : plt.subplot2grid(grid_size, (0, 0), rowspan=1, colspan=2)
# plot the rate of oxygen band events after filtering
hs.plotts(hkp, 'LC_OXYGEN', title='Oxygen rate', onlygood=True, 
        collapse=True, psym='.', color=pcolor, pdf=pdf)

if pdf : plt.subplot2grid(grid_size, (0, 2), rowspan=1, colspan=1)
# plot the oxygen rate versus Sun angle after filtering
hs.plotxy(hkp, 'SUN_ANGLE', 'LC_SCI', title='Science rate vs Sun angle',
          onlygood=True, color=pcolor, pdf=pdf,
          xbin=linspace(min(hk2['SUN_ANGLE'][g]), max(hk2['SUN_ANGLE'][g]), 16),
          dt=64.0, minbin=6)

if pdf : plt.subplot2grid(grid_size, (1, 0), rowspan=1, colspan=1)
# plot the oxygen rate versus proton flux after filtering
hs.plotxy(hkp, 'LC_HARD', 'LC_OXYGEN', title='Oxygen rate vs hard flux',
          onlygood=True, color=pcolor, pdf=pdf,
          xbin=linspace(0.0, hard_cut, 16), dt=64.0, minbin=6)

if pdf : plt.subplot2grid(grid_size, (1, 1), rowspan=1, colspan=1)
# plot the oxygen rate versus proton flux after filtering
hs.plotxy(hkp, 'LC_VLE', 'LC_OXYGEN', title='Oxygen rate vs VLE flux',
          onlygood=True, color=pcolor, pdf=pdf,
          xbin=linspace(0.0, vle_cut, 16), dt=64.0, minbin=6)

if pdf : plt.subplot2grid(grid_size, (1, 2), rowspan=1, colspan=1)
# plot the oxygen rate versus proton flux after filtering
hs.plotxy(hkp, 'PFLUX', 'LC_OXYGEN', title='Oxygen rate vs proton flux',
          onlygood=True, color=pcolor, pdf=pdf,
          xbin=linspace(0.0, max(hk2['PFLUX'][g]), 16), dt=64.0, minbin=6)

if pdf : plt.subplot2grid(grid_size, (2, 0), rowspan=1, colspan=1)
# plot the soft rate versus rigidity flux after filtering
hs.plotxy(hkp, 'COR_SAX', 'LC_SCI', title='Science rate vs rigidity', 
          onlygood=True, psym='.', color=pcolor, pdf=pdf,
          xbin=linspace(0.0, 7.0, 35), dt=64.0, minbin=6) 

## plot the oxygen rate versus rigidity flux after filtering
#hs.plotxy(hkp, 'COR_SAX', 'LC_OXYGEN', title='Oxygen rate vs rigidity', 
#          onlygood=True, psym='.', color=pcolor, pdf=pdf,
#          xbin=linspace(0.0, 6.0, 30), dt=64.0, minbin=6) 

#if pdf : plt.subplot2grid(grid_size, (1, 1), rowspan=1, colspan=1)
# plot the soft rate versus VLE rate after filtering
#hs.plotxy(hkp, 'LC_VLE', 'LC_SCI', title='Soft rate vs VLE rate', 
#          onlygood=True, psym='.', color=pcolor, pdf=pdf,
#          xbin=linspace(0.0, vle_cut, 16), dt=64.0, minbin=12) 
## plot the oxygen rate versus VLE rate * rigidity flux after filtering
#hs.plotxy(hkp, 'VLE_COR', 'LC_SCI', title='Soft rate vs VLE*rigidity', 
#          onlygood=True, psym='.', color=pcolor, pdf=pdf,
#          xbin=linspace(0.0, 15.0, 60), dt=64.0, minbin=12) 

if pdf : plt.subplot2grid(grid_size, (2, 1), rowspan=1, colspan=1)
# plot the science rate versus time since SAA after filtering
hs.plotxy(hkp, 'SAA_TIME', 'LC_SCI', title='Soft rate vs time since SAA',
          onlygood=True, xlim=[0.0,3000.0], pdf=pdf,
          psym='.', color=pcolor)

if pdf : plt.subplot2grid(grid_size, (2, 2), rowspan=1, colspan=1)
# plot the rate of oxygen band events versus nadir angle after filtering
hs.plotxy(hkp, 'NADIR_ANGLE', 'LC_OXYGEN', title='Oxygen rate vs nadir angle',
          onlygood=True, color=pcolor, pdf=pdf,
          xbin=linspace(min(hk2['NADIR_ANGLE'][g]), max(hk2['NADIR_ANGLE'][g]), 16), 
          dt=64.0, minbin=6)
#          onlygood=True, pdf=pdf, psym='.', color=pcolor)

## plot the science rate versus time since SAA after filtering
#hs.plotxy(hkp, 'SAA_TIME', 'LC_ALSI', title='Inst line rate vs time since SAA',
#          onlygood=True, xlim=[0.0,3000.0], pdf=pdf,
#          psym='+', color=pcolor)
#            xbin=linspace(0.0, 3300.0, 33), dt=64.0, minbin=12) 

#  # plot SAT_LAT vs SAT_LON after filtering
#  hs.plotxy(hk2, 'SAT_LON', 'SAT_LAT', title='Spacecraft position', 
#            onlygood=True, 
#            add=True, dpu=dpu, psym='.', color=pcolor[dpu], legend='best')

if pdf :    
  plt.tight_layout()
  pdff.savefig(fig2) # save the figure to the pdf file
  # create a figure that is 8 by 10.5 inches
  fig3 = plt.figure(outroot, figsize=(10.0, 7.5), dpi=100)
  plt.clf()
  plt.axis('off')


s =     ['Observation blocks                         Exposure (s)']
s.append('Start               Stop                   Att   uf14   uf38   uf54   cl14   cl38   cl54   Sun angle')
for i in range(len(atts)) :
  #print(atts[i])
  # find sun angles in observations  
  tc = "& (hk2['TIME'] > "+str(atts[i])+") & (hk2['TIME'] <= "+str(atte[i])+")" 
  #print('tc = ', filt + tc)
  gi = eval(filt + tc)
  #print(gi)
  if True in gi: #True means GTI survived cuts
  	sunmin = min(hk2['SUN_ANGLE'][gi])
  	sunmax = max(hk2['SUN_ANGLE'][gi])
  	t = str(getDate(atts[i]))[0:19]+' '+str(getDate(atte[i]))[0:19]+' '+'{0:6d}'.format(int(atttot[i]))
  	for j in range(6) : t += ' {0:6d}'.format(int(ftimes[i][j]))
  	t += '  {0:5.1f}:{1:5.1f}'.format(sunmin, sunmax)
  	s.append(t)
  else: #no data survived cuts, include a null value for sun angle to keep track of GTI existing. JKB 10/2022
  	sunmin=''
  	sunmax=''
  	t = str(getDate(atts[i]))[0:19]+' '+str(getDate(atte[i]))[0:19]+' '+'{0:6d}'.format(int(atttot[i]))
  	for j in range(6) : t += ' {0:6d}'.format(int(ftimes[i][j]))
  	t += ''
  	s.append(t)
  #print sunmin, sunmax
  #t = str(getDate(atts[i]))[0:19]+' '+str(getDate(atte[i]))[0:19]+' '+ \
   # '{0:6d}'.format(int(atttot[i]))
  #for j in range(6) : t += ' {0:6d}'.format(int(ftimes[i][j]))
  #t += '  {0:5.1f}:{1:5.1f}'.format(sunmin, sunmax)
  #s.append(t)

if pdf :    
  hp = 0.05
  monofont = {'fontname':'DejaVu Sans Mono'}
  for i in range(len(s)) :
    plt.text(hp, 0.96-0.03*i, s[i], transform=fig1.transFigure, **monofont)
else :
  print(s)

if pdf :
  plt.tight_layout() # fit plots together
  pdff.savefig(fig3) # save the figure to the pdf file
  pdff.close() # write the pdf

sys.stdout = orig_stdout # set stdout back the way it started
f.close() # close the output file
