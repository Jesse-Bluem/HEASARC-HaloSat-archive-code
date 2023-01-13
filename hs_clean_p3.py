# Routines for HaloSat data processing and filtering
# P. Kaaret 11/2/2019-5/2/2019, 5/31/2019-7/30/2019, 9/14-20/2019
# R. Ringuette 1/22/2019, 4/25/2019, 5/1/2019, 7/26/2019
# J. Bluem 10/16/2022 - conversion for Python 3. This version is called by the HEASARC file processing code in the Server_Code directory. Removed fmt='+', which threw a warning message every time. Added 32 seconds to DATE-OBS and DATE-END to correct for time issues. Quality flag moved by 1 channel to match edits for gain shift.

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.path as path
from matplotlib.lines import Line2D
from quaternion import Quaternion
import time, glob, os
from astropy.coordinates import SkyCoord
import astropy.units
from datetime import datetime, timedelta


def getDate(time) :
  epoch = datetime(2000, 1, 1, 0, 0, 0) # reference epoch yr = 2000.0 for BCT EDU
  # epoch = utc.localize(epoch) # epoch localized to UTC
  time1 = time - 37.0
  date = epoch + timedelta(seconds=time1)
  return date


def equatorial2ecliptic(ra, dec):
# Converts a list of equatorial coordinates (ra,dec) to galactic (l,b).
  # Insert equatorial coordinates into Astropy SkyCoord
  # set distance to 1 Mpc to avoid any affects from shifting reference point
  equ = SkyCoord(ra=ra*astropy.units.degree, dec=dec*astropy.units.degree, \
                 frame='icrs', distance=1*astropy.units.Mpc)
  # Get the ecliptic coordinates
  ecl = equ.heliocentrictrueecliptic
  elon = ecl.lon.value # longitude
  elat = ecl.lat.value # latitude
  return (elon, elat)


# plotting functions
def plothist(hk, v, onlygood=False, add=False, dpu=0, \
             nbin=512, prange=[], color='', plot_bar=False, \
             title='', xlabel='', ylabel='', legend='', pdf=False) :
# plot a histogram of HK or light curve data
  # if user did not specify a figure title then use variable v
  if len(title) == 0 : title = 'Histogram of '+v
  if not pdf : plt.figure(title) # move to the appropriate plot window
  if (not add) or (not plt.fignum_exists(title)) : plt.clf()
  # ylabel = 'Number/bin' # !!!rr this overwrite any user ylabel
  # label for data is either variable name, variable name 'for GTI', or DPU number
  if onlygood :
    x = np.array(hk[v][hk['GOOD']])
    vlabel = v+' for GTI'
    # ylabel+=' for GTI' # !!!rr why add this to y? seems redundant
  else :
    x = np.array(hk[v])
    vlabel = v
  # add DPU to label if plotting multiple DPUs
  if dpu != 0 : vlabel += ' DPU '+str(dpu)
  # set range for x axis
  if len(prange) != 2 : prange = [min(x), max(x)]
  # do the plotting
  if len(color) > 0 :
    plt.hist(x, bins=nbin, range=prange, label=vlabel, color=color, histtype='step')
  else :
    plt.hist(x, bins=nbin, range=prange, label=vlabel, histtype='step')
  if not pdf : plt.title(title) # set plot title
  if len(xlabel) == 0 : plt.xlabel(v) # set the X-axis label
  else : plt.xlabel(xlabel)
  if len(ylabel) == 0 : plt.ylabel('Number/bin') # set the X-axis label
  else : plt.ylabel(ylabel) 
  if plot_bar: plt.axvline(x=plot_bar,color='black',linestyle='--')
  if len(legend) > 0 : plt.legend(loc=legend) # make a legend if requested
  if not pdf : plt.show()


def orbbreak(t, y, dpu):
  sort_idx = np.argsort(t) # sort data in order of time
  t = t[sort_idx]
  break_idx = np.where(np.diff(t) > 2000)[0]
  ocolor = {14:'lightblue', 54:'palegreen', 38:'pink'} # !!!rr match user colors?
  for i in range(len(break_idx)):
    plt.plot([break_idx[i],break_idx[i]], [min(y),max(y)], '--', color=ocolor[dpu])

def plotts(hk, v, collapse=True, orb_break=False, dt=64.0, onlygood=False, 
           psym='.', color={14:'r', 54:'g', 38:'b'}, title='', ylabel='', ylim=None, 
           plot_bar=None, legend='', hist=False, nbin=512, pdf=False) :
# plot a time series of HK or light curve data
  # collapse=True ignore gaps in light curve, =False plots against time
  # orb_break=True plots lines between orbits

  # figure out what we should call the plot
  # add to title and figure name if plotting only GTI
  if onlygood : gs = ' for GTI'
  else : gs = ''
  # if user did not specify a figure title then use variable name
  if len(title) == 0 : title = v+gs
  else : title = title+gs
  if not pdf : 
    plt.figure(title) # move to the appropriate plot window
    plt.clf()
  # should we plot on the the 'good' data?
  if onlygood : 
    y = hk[v][hk['GOOD']] # plot only data in GTI
    t = hk['TIME'][hk['GOOD']]
    d = hk['DPU_ID'][hk['GOOD']]
  else : 
    y = hk[v]  # plot all data
    t = hk['TIME']          
    d = hk['DPU_ID']
  # make array with color set by DPU for each point
  c = np.vectorize(color.get)(d)
  # for i in range(len(d)) : c[i] = color[d[i]]
  # print 'c = ', c[0:10]
  # did the user specify limits on y?
  if ylim is not None : 
    plt.ylim(ylim) # set the plot limits in y
    # set points with y > upper limit to the upper limit
    y = y.copy() # make a copy so we don't overwrite data in the dictionary
    q = np.where(y > ylim[1])
    if len(q) > 0 : y[q] = ylim[1]
  # should we collapse the time axis?
  if collapse :
    # sort arrays in time order
    q = np.argsort(t)
    y, t, c = y[q], t[q], c[q]
    # find gaps in time, include one at start and one past end
    g = np.concatenate(([t[0]], t[np.where(np.diff(t) > 2000.0)[0]+1], [np.max(t)+dt]))
    tn = 0.0*t # new time array with gaps removed
    gi, toff, tint = 0, 0.0, 0.0
    for i in range(len(t)) :
      if t[i] >= g[gi+1] :
        gi += 1
        toff += tint + 2*dt
        tint = 0.0
      tint = max(tint, t[i] - g[gi]) # update interval duration, use of max is probably redundant
      tn[i] = t[i] - g[gi] + toff
    td = np.diff(tn)
    plt.scatter(tn/1E3, y, s=5.0, c=c, marker=psym) # plot the data
    plt.xlabel('Time (ks)') # set the X-axis label
    # if orb_break: orbbreak(t,y,dpu) # only plot the orbit breaks if collapsed and if requested
    title += ' vs elapsed time'
  else:
    # plt.plot(t, y, psym, label=dlabel, ms=2.0) # plot the data
    plt.scatter(t/1E3-min(t)/1E3, y, s=2.0, c=c, marker=psym) # plot the data
    plt.xlabel('Time (ks)') # set the X-axis label
    # obvious without orbbreak to tell where the orbital breaks are
    title += ' vs time'
  # should we plot a horizontal histograms?
  if hist :
    print('nbin = ', nbin)
    # !!! need to scale x-axis to make this visible
    if ylim is None : ylim = [min(y), max(y)] # set range for histograms
    for dpu in np.unique(d) :  # make a histogram for each dpu
      plt.hist(y[np.where(d == dpu)], bins=nbin, range=ylim, color=color[dpu],
               histtype='step', orientation='horizontal', linewidth=1)
  # pretty up the plot with labels, etc.
  if not pdf : plt.title(title) # set plot title
  if len(ylabel) == 0 : plt.ylabel(v) # set the Y-axis label
  else : plt.ylabel(ylabel)
  if len(legend) > 0 :
    # !!! should make this work if not all three DPUs being used
    s = [Line2D([0], [0], color=color[14], lw=1),
         Line2D([0], [0], color=color[54], lw=1),
         Line2D([0], [0], color=color[38], lw=1)]
    plt.legend(s, ['14', '54', '38'], loc=legend)
    #plt.legend([str(k) for k in color.keys()], loc=legend) # make a legend if requested
  if plot_bar is not None : plt.axhline(y=plot_bar,color='black',linestyle='--')
  if not pdf : plt.show()


def linear_fit(x, y, yerr) :
  """
  Find the best slope (a) and intercept (b) that minimize the chi-sqaured in fitting
  a line to a set of data points (x,y) with error on y.
  Note y = a*x + b
  Returns a, error on a, b, error on b, chi squared, degrees of freedom
  """

  # Sums we will use to find slope and intercept that minimize chi-squared
  sumy = np.sum(y/(yerr**2))
  sum1 = np.sum(1/(yerr**2))
  sumx = np.sum(x/(yerr**2))
  sumxy = np.sum((x*y)/(yerr**2))
  sumx2 = np.sum((x**2)/(yerr**2))
  delta = ((sum1*sumx2) - (sumx**2))
  # find slope
  a = ((sum1*sumxy) - (sumx*sumy))*(1/delta)
  # find intercept
  b = ((sumx2*sumy) - (sumx*sumxy))*(1/delta)
  # error on slope and intercept
  aerr = np.sqrt((sum1)*(1/delta))
  berr = np.sqrt((sumx)*(1/delta))
  # chi squared deviation of data from fitted line
  chi2 = np.sum(((y-(a*x+b))/yerr)**2)
  # calculate the number of data points
  num_points = np.size(x)
  # degrees of freedom  is 2 fewer, since in the linear fit we have two parameters
  dof = num_points - 2
  return a, aerr, b, berr, chi2, dof

def plotxy(hk, x, y, onlygood=False, \
           psym='.', color='k', title='', xlabel='', ylabel='', legend='', \
           xlim=None, ylim=None, xbin=None, dt=8.0, minbin=6, pdf=False) :
# plot one HK or light curve data set versus another
  # add to title and figure name if plotting only GTI
  if onlygood : gs = ' for GTI'
  else : gs = ''
  # if user did not specify a figure title then use variable names
  if len(title) == 0 : title = y+' vs '+x+gs
  else : title = title+gs
  if not pdf : # if this plot is in its own window then
    fig = plt.figure(title) # move to the appropriate plot window
    plt.clf() # clear the plot window

  xd, yd = hk[x], hk[y] # selected X and Y data
  d = hk['DPU_ID'] # DPU for each data point
  if onlygood : # use only the data passing cuts
    p = hk['GOOD'] # find the good data
    xd, yd, d = xd[p], yd[p], d[p]
    #xd, yd, d = hk[x][p], hk[y][p], d[p]
  #else : # use all the data
  #  xd, yd = hk[x], hk[y]
  # make array with color for each point set by DPU
  c = np.vectorize(color.get)(d)

  if xbin is not None : # use the specified x bins to bin the data
    n = len(xbin)-1 # number of bins
    xb, cb, nb = [], [], [] # lists to hold the binned data
    for i in range(n) : # loop over the bins
      # find data points in this bin
      q = (xbin[i] <= xd) & (xd < xbin[i+1])
      if np.sum(q) >= minbin : # make a bin if at least minbin data points
        xb.append(np.mean(xd[q])) # mean of x
        cb.append(np.sum(dt*yd[q])) # total counts, dt is conversion from rate to counts
        nb.append(len(xd[q])) # number of data point in this bin
    if len(xb) > 0 : # are there any data points?
      # convert lists to arrays
      xb = np.array(xb)
      yb = np.array(cb)/(dt*np.array(nb)) # rate in bin
      yberr = (np.sqrt(np.array(cb)+1.0)+1.5)/(dt*np.array(nb)) # Gehrels error on rate
      # plot the binned data
      plt.errorbar(xb, yb, yberr, marker='+', markersize=3, 
                   color='k', linestyle='None') #JKB removed fmt='+' Python warns its already defined by marker for every loop, lets clean it up.
      # find the weighted average and error
      b = np.sum(yb/(yberr**2))/np.sum(1/(yberr**2))
      berr = np.sqrt(1/np.sum(1/(yberr**2)))
      chi2 = sum(((yb - b)/yberr)**2)
      dof = len(yb)-1 # degrees of freedom
      # print 'Weighted average = ', b, berr, chi2, dof
      # plot the weighted average
      xf = np.linspace(min(xb), max(xb))
      plt.plot(xf, 0*xf+b, ':', color='b')
      if dof > 1 : # avoid divide by zero error
        plt.legend([r'$\chi^2/\nu$'+' = {0:.2f}'.format(chi2/dof)], loc='upper left')
      ## do a linear fit to the binned data
      #(a, aerr, b, berr, chi2, dof) = linear_fit(xb, yb, yberr)
      #print a, aerr, b, berr, chi2, dof
      # plot the fit
      #xf = np.linspace(min(xb), max(xb))
      #plt.plot(xf, a*xf+b, ':', color='b')
      #n = -int(np.floor(np.log10(abs(a)))) # negative of exponent base 10 fo slope for nicer formatting
      #plt.text(0.1, 0.9, 'Slope*10^{0:d} = {1:.2f}+/-{2:.2f}'.format(n, a*10**n, aerr*10**n), transform=ax.transAxes)
      #plt.text(0.1, 0.8, 'Chi^2/DoF = {0:.2f}'.format(chi2/dof), transform=ax.transAxes)
    else :
      xbin = None # didn't have enough data to bin, so try plotting without bins
  if xbin is None : # plot individual data points
    # plt.plot(xd, yd, marker=psym, color=color, label=dlabel, linestyle='None')
    plt.scatter(xd, yd, s=2.0, c=c, marker=psym) # plot the data
    if len(legend) > 0 :
      # !!! should make this work if not all three DPUs being used
      s = [Line2D([0], [0], color=color[14], lw=1),
           Line2D([0], [0], color=color[54], lw=1),
           Line2D([0], [0], color=color[38], lw=1)]
      plt.legend(s, ['14', '54', '38'], loc=legend)
  
  # set x and y axis ranges if specified by the user
  if xlim is not None : plt.xlim(xlim)
  if ylim is not None : plt.ylim(ylim)
  if len(xlabel) == 0 : plt.xlabel(x) # set the X-axis label
  else : plt.xlabel(xlabel)
  if len(ylabel) == 0 : plt.ylabel(y) # set the Y-axis label
  else : plt.ylabel(ylabel)
  if not pdf : 
    plt.title(title) # set plot title
    plt.show()


def find_pointing(inroot, hs_dir='') :
  # get the nominal pointing coordinates from the RA_NOM and DEC_NOM keywords
  # in the attitude file

  hdulist = pyfits.open(inroot+'.att') # read attitude file
  header = hdulist['ATTITUDE'].header # get header of ATTITUDE extension
  hdulist.close() # close the file
  ra_nom = header['RA_NOM'] # get the coordinates
  dec_nom = header['DEC_NOM']
  print('Target coordinates = ', ra_nom, dec_nom)
  return [ra_nom, dec_nom]


def read_omni(omnifile, hs_dir='') :
  """
  Read solar wind data and return lists of: time, speed, density.
  Solar wind data are described at https://omniweb.gsfc.nasa.gov/ow_min.html
  Download omni_5min2018.asc and omni_5min2018.asc from
  https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/
  wind speed [km/s] and proton density [#/cm^3] are tabulated versus time
  proton flux [#/cm^2/s] = 1E5*speed*density
  Henley and Shelton (2010, ApJS 187, 388) removed flux > 2E8 cm^-2 s^-1.
  """
  t0 = 568080000.0+37 # TAI seconds at 2018-01-01 following BCT's conventions
  f = open(hs_dir+omnifile) # open the solar wind data file which is in text format
  time, speed, density = [], [], [] # lists to hold data
  for s in f : # read each line in the file in turn
    # data are in fixed columns
    # print s[0:14], s[124:131], s[156:162]
    year = int(s[0:4]) # find date and time of observation
    day = int(s[5:8]) # day of year (1-365)
    hour = int(s[9:11]) # hour of day (0-23)
    minute = int(s[12:14]) # minute of hour (0-59)
    # convert to seconds since 2018-01-01 and then into TAI seconds since 2000-01-01
    tai = t0 + (year-2018)*31536000 + (day-1)*86400 + hour*3600 + minute*60
    time.append(tai)
    speed.append(float(s[124:131])) # flow speed [km/s]  124:131
    density.append(float(s[156:162])) # proton density, n/cc  156:162
    # print year, day, hour, minute, time[-1], speed[-1], density[-1]
  f.close() # close the text data file
  return (time, speed, density)


def read_uf(inroot, dpu, timecuts=[], hs_dir='', loud=0) :
  """
  Read HaloSat data from housekeeping (.hk), spacecraft (.bus), and
  event (.evt) files.  Convert into filter data with uniform time bins.
  Returns dictionary of data for filtering, including light curves.
  """

  # load housekeeping data into hk dictionary
  hk = {} # make empty dictionary for housekeeping and pointing/orbit data
  hkfile = inroot+'_s'+str(dpu)+'.hk'
  print('Reading ', hkfile)
  f = pyfits.open(hkfile) # read in the FITS file with HK data
  dhk = f[1].data # put all of the data into the FITS record array d
  f.close() # close the FITS file
  # copy HK data into hk dictionary
  for k in dhk.names : # loop over all columns in the FITS record array
    hk[k] = np.array(dhk[k]) # enter the column into the dictionary
  if len(hk['TIME']) < 1 : # no time bins
    print('No time bins found in housekeeping data')
    return hk
  if (loud > 2) : print('HK time range = ', min(hk['TIME']), max(hk['TIME']))
  # add DPU_ID to dictionary
  hk['DPU_ID'] = dpu + 0*hk['IN_SAA'] # use IN_SAA because it is integer array

  # Read solar wind data into hk dictionary
  # read the solar wind data for 2018 and 2019
  # !!! will need to add 2020
  swtime1, speed1, density1 = read_omni('omni_5min2018.asc', hs_dir=hs_dir)
  swtime2, speed2, density2 = read_omni('omni_5min2019.asc', hs_dir=hs_dir)
  # convert to arrays
  swtime = np.array(swtime1 + swtime2)  # make numpy arrays from the lists
  speed = np.array(speed1 + speed2)
  density = np.array(density1 + density2)
  # select times with valid data (speed = 99999.9 indicates bad data)
  q = np.where(speed < 10000.0)
  swtime = swtime[q]
  # calculate proton flux [#/cm^2/s] = 1E5*speed*density
  # use units of 10^8 cm^-2 s^-1
  protonflux = 1E-3*speed[q]*density[q]
  # interpolate to get proton flux at times of HK records
  pflux = np.interp(hk['TIME'], swtime, protonflux, right=0) # min(protonflux))
  # copy into hk dictionary
  hk['PFLUX'] = pflux # proton flux in units of 10^8 cm^-2 s^-1

  # load event data into evt dictionary
  uffile = inroot+'_s'+str(dpu)+'_uf.evt'
  print('Reading ', uffile)
  f = pyfits.open(uffile) # read in the FITS file with event data
  d = f[1].data # put all of the data into the dictionary d
  f.close() # close the FITS file
  # copy event data into evt dictionary
  evt = {} # make dictionary
  evt['TIME'] = np.array(d['TIME']) # event time
  evt['PHA'] = np.array(d['PHA']) # event pulse height amplitude
  evt['PI'] = np.array(d['PI']) # event pulse height invariant

  # if timecuts are set, return only data within selected time ranges
  if len(timecuts) == 0 : # no time cuts
    # return the HK and event data
    return hk, evt
  else : # have time cuts
    # select HK data in given time intervals
    time0 = np.array(hk['TIME']) # find the time for each bus record
    q = np.zeros(len(time0), dtype=bool) # fill array with false
    for t0, t1 in timecuts : # loop over selected time intervals
      # set q0 true for records in this time interval
      q = q | ((t0 <= time0) & (time0 <= t1))
    # write new dictionary with HK data from selected time intervals
    hkc = {} # make empty dictionary
    # loop over items and copy entries in selected time range to new dictionary
    for k, v in hk.items() : hkc[k] = v[q]
    # select event data in given time intervals
    time0 = np.array(evt['TIME']) # find the time for each bus record
    q = np.zeros(len(time0), dtype=bool) # fill array with false
    for t0, t1 in timecuts : # loop over selected time intervals
      # set q0 true for records in this time interval
      q = q | ((t0 <= time0) & (time0 <= t1))
    # write new dictionary with evt data from selected time intervals
    evtc = {} # make empty dictionary
    # loop over items and copy entries in selected time range to new dictionary
    for k, v in evt.items() : evtc[k] = v[q]
    # return the seelcted HK and event data
    return hkc, evtc
  

def time_rebin(hk, gti_min=0, gti_max=5400, loud=0) :
  """
  Find good time intervals (GTIs).
  """
  dhk = 1.5*8.0 # maximum allowed time between succesive HK records
  # !!! try decreasing after Jesse re-runs database with new time procedure

  hk_time = hk['TIME'] # array with times of HK intervals

  # find continuous good intervals
  gti_start, gti_stop = [hk_time[0]], [hk_time[0]]
  n = len(hk_time) # number of time bins
  for i in range(0, n) : # !!! start at 1?
    # print i, time[i], gti_start[-1], gti_stop[-1]
    # check if next time bin is continuous with current GTI and is good
    if (hk_time[i] - gti_stop[-1]) < dhk :
      gti_stop[-1] = hk_time[i] # then extend current GTI
      # if GTI is longer than gti_max, then start a new GTI
      if (gti_stop[-1]-gti_start[-1]) >= gti_max :
        gti_start.append(hk_time[i])
        gti_stop.append(hk_time[i])    
    else : # if we're not extending this GTI, then start a new GTI
      gti_start.append(hk_time[i])
      gti_stop.append(hk_time[i])
  # change lists to arrays
  gti_start = np.array(gti_start)
  gti_stop = np.array(gti_stop)

  # remove GTI's that are too short
  dur = gti_stop-gti_start
  q = dur >= gti_min
  gti_start = gti_start[q]
  gti_stop = gti_stop[q]

  if len(gti_start) == 0 :
    print('*** No good time intervals found.  Please check your cuts or your data files.')
    return hk

  if (loud > 1) :
    # print GTIs
    for i in range(len(gti_start)) :
      print('GTI = ', gti_start[i], gti_stop[i], gti_stop[i]-gti_start[i])

  hk2 = {} # new hk dictionary with different time binning
  # loop over all of the items in the HK dictionary
  for k, v in hk.items(): # k = key, v = value of key
    a = [] # list to hold average value for each new time bin
    for i in range(len(gti_start)) : # loop over GTIs
      # find HK bins in this GTI
      q = (gti_start[i] <= hk_time) & (hk_time < gti_stop[i])
      # find average of hk values
      a.append(np.mean(v[q]))
    hk2[k] = np.array(a) # add array of average values
  # add start and stop times for each interval
  hk2['TIME'] = gti_start
  hk2['TIME_END'] = gti_stop
  
  # return the new hk2 dictionay and the gti info
  return hk2


def write_evt_pi(evt, dpu, gti, inroot, outroot, filt, loud=0, group=False, ratekey=False) :
  """
  Write clean event and spectrum files.
  evt = dictionary of event data
  dpu = DPU number (integer)
  gti = 2*N array, gti[0] = good interval start times, gti[1] = stop times
  inroot = root of input (unfiltered) file names
  outroot = root of output (clean) file names
  loud - adjusts amount of output, 0=minimum, 
  group = False for no group, otherwise minimum number of counts per bin
  if ratekey is set, then a keyword is added with the rate in the band, 
    e.g. ratekey=['LC_HARD', 5600, 12900]
  """

  # check that we have at least one GTI
  if len(gti[0]) < 1 :
    print('No GTI passed to write_spectrum.')
    print('You may want to check your cuts.')
    exit(1)

  # copy event data from dictionary to arrays
  evt_time = evt['TIME'] # event time
  evt_pha = evt['PHA'] # event pulse height amplitude
  evt_pi = evt['PI'] # event pulse height amplitude
  #print 'Number of events = ', len(evt_time)

  # filter events
  q = np.zeros((len(evt_time)), dtype=bool) # q = True means accept event
  ontime = 0.0
  for t in np.transpose(gti) :
    # if the event time is in this filter interval, then set q to True
    q = q | ((t[0] <= evt_time) & (evt_time <= t[1]))
    ontime += t[1]-t[0] # have to recalculate ontime since filtering
    #print 't, ontime, sum(q) = ', t[0], t[1], ontime, np.sum(q)
  # filter events
  evt_time = evt_time[q]
  evt_pha = evt_pha[q]
  evt_pi = evt_pi[q]
  #print 'Number of events after filtering = ', len(good_pha)

  # read unfiltered event file to get header
  uffile = inroot+'_s'+str(dpu)+'_uf.evt'
  print('Reading ', uffile)
  f = pyfits.open(uffile) # open the unfiltered file
  uf_prihdr = f[0].header # get primary header
  uf_evthdr = f['EVENTS'].header # get event header
  uf_gtihdr = f['GTI'].header # get GTI header
  f.close() # close the FITS file

  # write out clean event file
  # make the FITS header
  # make the primary header and copy header keywords from unfiltered file
  prihdu = pyfits.PrimaryHDU(header=uf_prihdr)
  # modify keywords that need to change
  prihdu.header['DATE-OBS'] = (str(getDate(min(gti[0])+32)).replace(' ', 'T')[0:19], 'Start date of observations') #QQQ
  prihdu.header['DATE-END'] = (str(getDate(max(gti[0])+32)).replace(' ', 'T')[0:19], 'End date of observations') #QQQ
  prihdu.header['DATE'] = (str(datetime.now()).replace(' ', 'T')[0:19], 'File creation date')

  # define the columns for events
  col1 = pyfits.Column(name='TIME', format='1D', array=evt_time) # event time
  col2 = pyfits.Column(name='PHA', format='1I', array=evt_pha) # event PHA
  col3 = pyfits.Column(name='PI', format='1I', array=evt_pi) # event PI
  cols = pyfits.ColDefs([col1, col2, col3]) # create a ColDefs (column-definitions) object for all columns:
  # create a binary table HDU for the event data, copy the header from the unfiltered file
  evthdu = pyfits.BinTableHDU.from_columns(cols, header=uf_evthdr)
  # modify keywords that need to change content or location
  evthdu.header.remove('TTYPE1') # remove keyword in wrong place
  evthdu.header.insert(8, ('TTYPE1', 'TIME', 'Time of event')) # insert in right place
  evthdu.header.remove('TFORM1') # remove keyword in wrong place
  evthdu.header.insert(9, ('TFORM1', '1D', 'data format of field'))
  #evthdu.header.remove('TUNIT1') # remove keyword in wrong place
  evthdu.header.insert(10, ('TUNIT1', 's', 'physical unit of field'))
  evthdu.header.remove('TTYPE2') # remove keyword in wrong place
  evthdu.header.insert(11, ('TTYPE2', 'PHA', 'Pulse Height Analyzer'))
  evthdu.header.remove('TFORM2') # remove keyword in wrong place
  evthdu.header.insert(12, ('TFORM2', '1I', 'data format of field'))
  #evthdu.header.remove('TUNIT2') # remove keyword in wrong place
  evthdu.header.insert(13, ('TUNIT2', 'chan', 'physical unit of field'))
  evthdu.header.remove('TTYPE3') # remove keyword in wrong place
  evthdu.header.insert(16, ('TTYPE3', 'PI', 'Pulse Invariant'))
  evthdu.header.remove('TFORM3') # remove keyword in wrong place
  evthdu.header.insert(17, ('TFORM3', '1I', 'data format of field'))
  #evthdu.header.remove('TUNIT3') # remove keyword in wrong place
  evthdu.header.insert(18, ('TUNIT3', 'chan', 'physical unit of field'))
  evthdu.header.insert(21, ('TNULL3', -1))
  evthdu.header['HDUCLAS2'] = ('ACCEPTED', 'Second Class level')
  evthdu.header['DATE-OBS'] = (str(getDate(min(gti[0])+32)).replace(' ', 'T')[0:19], \
  	                           'Start date of observations')
  evthdu.header['TSTART'] = (min(gti[0]), '[s] Observation start time')
  evthdu.header['DATE-END'] = (str(getDate(max(gti[1])+32)).replace(' ', 'T')[0:19], \
  	                           'End date of observations')
  evthdu.header['TSTOP'] = (max(gti[1]), '[s] Observation stop time')
  evthdu.header['TELAPSE'] = (max(gti[1])-min(gti[0]), '[s] Stop - Start')
  evthdu.header['ONTIME'] = (ontime, '[s] Observation time on target')
  evthdu.header['EXPOSURE'] = (ontime, '[s] Exposure')
  evthdu.header['PROCVER'] = (evthdu.header['PROCVER']+'_hscl_20221026', \
                              'Processing script version number')
  evthdu.header['CREATOR'] = ('db_hsuf_hscl', 'Software creator of the file')
  evthdu.header['DATE'] = (str(datetime.now()).replace(' ', 'T')[0:19], \
  	                       'File creation date')

  # define the columns for GTI
  col1 = pyfits.Column(name='START', format='1D', array=gti[0])
  col2 = pyfits.Column(name='STOP', format='1D', array=gti[1])
  cols = pyfits.ColDefs([col1, col2]) # create a ColDefs (column-definitions) object for all columns:
  # create a binary table HDU for the GTI data, copy the header from the unfiltered file
  gtihdu = pyfits.BinTableHDU.from_columns(cols, header=uf_gtihdr)
  # modify keywords that need to change content or location
  gtihdu.header.remove('TTYPE1') # remove keyword in wrong place
  gtihdu.header.insert(8, ('TTYPE1', 'START', 'Start time')) # insert in right place
  gtihdu.header.remove('TFORM1') # remove keyword in wrong place
  gtihdu.header.insert(9, ('TFORM1', '1D', 'data format of field'))
  #evthdu.header.remove('TUNIT1') # remove keyword in wrong place
  gtihdu.header.insert(10, ('TUNIT1', 's', 'physical unit of field'))
  gtihdu.header.remove('TTYPE2') # remove keyword in wrong place
  gtihdu.header.insert(11, ('TTYPE2', 'STOP', 'Stop time'))
  gtihdu.header.remove('TFORM2') # remove keyword in wrong place
  gtihdu.header.insert(12, ('TFORM2', '1D', 'data format of field'))
  #evthdu.header.remove('TUNIT2') # remove keyword in wrong place
  gtihdu.header.insert(13, ('TUNIT2', 's', 'physical unit of field'))
  gtihdu.header['HDUCLAS2'] = ('STANDARD', 'Second Class level')
  gtihdu.header['DATE-OBS'] = (str(getDate(min(gti[0])+32)).replace(' ', 'T')[0:19], \
  	                           'Start date of observations')
  gtihdu.header['TSTART'] = (min(gti[0]), '[s] Observation start time')
  gtihdu.header['DATE-END'] = (str(getDate(max(gti[1])+32)).replace(' ', 'T')[0:19], \
  	                           'End date of observations')
  gtihdu.header['TSTOP'] = (max(gti[1]), '[s] Observation stop time')
  #gtihdu.header['TELAPSE'] = (max(gti[1])-min(gti[0]), '[s] Stop - Start')
  #gtihdu.header['ONTIME'] = (ontime, '[s] On source time')
  #gtihdu.header['EXPOSURE'] = (ontime, '[s] Exposure')
  gtihdu.header['PROCVER'] = (gtihdu.header['PROCVER']+'_hscl_20200220', \
                              'Processing script version number')
  gtihdu.header['CREATOR'] = ('db_hsuf_hscl', 'Software creator of the file')
  gtihdu.header['DATE'] = (str(datetime.now()).replace(' ', 'T')[0:19], \
  	                       'File creation date')

  # define the columns for screening extension
  screen1 = ['HK_SDD'+str(dpu)]
  screen2 = [filt] # expression used for filtering
  col1 = pyfits.Column(name='EXTENSION', format='20A', array=screen1)
  col2 = pyfits.Column(name='EXPRESSION', format='600A', array=screen2)
  cols = pyfits.ColDefs([col1, col2]) # create a ColDefs (column-definitions) object for all columns:
  # create a binary table HDU for the screening extension
  scrhdu = pyfits.BinTableHDU.from_columns(cols)
  # modify keywords that need to change
  scrhdu.header['TTYPE1'] = ('EXTENSION', 'Name of extension to apply screening')
  scrhdu.header['TFORM1'] = ('20A', 'data format of field')
  scrhdu.header['TTYPE2'] = ('EXPRESSION', 'Expression')
  scrhdu.header['TFORM2'] = ('600A', 'data format of field')
  scrhdu.header['EXTNAME'] = ('SCREENING', 'Binary table extension name')
  scrhdu.header['COMMENT'] = 'Screening is done with 64 second bins, factor 8 rebinning.'
  scrhdu.header['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
  scrhdu.header['INSTRUME'] = ('SDD'+str(dpu), 'Instrument name')
  scrhdu.header['OBSERVER'] = ('PHILIP KAARET', 'Principal Investigator')
  scrhdu.header['OBS_ID'] = (uf_evthdr['OBS_ID'], uf_evthdr.comments['OBS_ID']) # 'Observation ID'
  scrhdu.header['OBJECT'] = (uf_evthdr['OBJECT'], uf_evthdr.comments['OBJECT']) # 'Object/target name'
  scrhdu.header['OBJTYPE'] =(uf_evthdr['OBJTYPE'], uf_evthdr.comments['OBJTYPE']) # 'Object type observed'
  scrhdu.header['ORIGIN'] = ('UNIVERSITY OF IOWA', 'Origin of FITS file')
  scrhdu.header['PROCVER'] = (uf_evthdr['PROCVER'], uf_evthdr.comments['PROCVER'])
  # ('dohs_20200124', 'Processing script version number')
  scrhdu.header['SOFTVER'] = (uf_evthdr['SOFTVER'], uf_evthdr.comments['SOFTVER'])
  # ('searchdb_20200124', 'Software version')
  scrhdu.header['CALDBVER'] = (uf_evthdr['CALDBVER'], uf_evthdr.comments['CALDBVER'])
  # ('hsYYYYMMDD', 'CALDB index version used')
  scrhdu.header['TLM2FITS'] = (uf_evthdr['TLM2FITS'], uf_evthdr.comments['TLM2FITS'])
  # ('hkdb_20200124_scidb_20200124', 'Telemetry converter version number')
  scrhdu.header['CREATOR'] = ('db_hsuf_hscl', 'Software creator of the file')
  scrhdu.header['DATE'] = (str(datetime.now()).replace(' ', 'T')[0:19], 'File creation date')

  # join primary with extensions
  hdulist = pyfits.HDUList([prihdu, evthdu, gtihdu, scrhdu])
  # write event file
  print('Writing ', outroot+'_s'+str(dpu)+'_cl.evt')
  hdulist.writeto(outroot+'_s'+str(dpu)+'_cl.evt', overwrite=True)


  ### write spectrum (.pi) file
  # set range and step of PI values, 
  # PI starts at energy = 0, increases in steps of ebin, so energy = pi*ebin
  # this must match program that wrote FITS file
  # these definitions must match with the arf and rmf
  emin = 0.1 # [keV] minimum energy value
  emax = 9.2 # [keV] maximum energy value
  ebin = 0.02 # [keV] energy bin size for PI
  nbin = int((emax-emin)/ebin) # number of bins in histogram
  #print 'Number of bins in spectrum = ', nbin

  # write the output PI file
  # make histogram from PI values
  (s, binedge) = np.histogram(evt_pi, bins=nbin, range=(0, nbin))
  # se = np.sqrt(s+1.0)+1.5 # statistical error via N. Gehrels 1986, ApJ 303, 336

  # definitions for bad data
  emindub = 0.28 # bins below this energy are marked dubious, qual = 2 #JKB 12/14 edit
  imindub = int(np.ceil((emindub - emin)/ebin)) # corresponding index
  emaxdub = 7.18 # bins above this energy are marked dubious, qual = 2 #JKB 12/14 edit 
  imaxdub = int(np.floor((emaxdub - emin)/ebin)) # corresponding index
  # print emindub, imindub, emaxdub, imaxdub
  # quality array 0 = good, 1 = bad, 2 = dubious
  qual = 0*s # make quality array, same size as spectrum, default is good
  # mark energies that are too low or too high as bad
  if imindub > 0 : qual[0:imindub] = 2
  if imaxdub < nbin-1 : qual[imaxdub:nbin] = 2

  # make the FITS header
  # make the primary header and copy header keywords from unfiltered file
  prihdu = pyfits.PrimaryHDU(header=uf_prihdr)
  # modify keywords that need to change
  prihdu.header['DATE-OBS'] = (str(getDate(min(gti[0])+32)).replace(' ', 'T')[0:19], 'Start date of observations') #QQQ
  prihdu.header['DATE-END'] = (str(getDate(max(gti[0])+32)).replace(' ', 'T')[0:19], 'End date of observations') #QQQ
  prihdu.header['DATE'] = (str(datetime.now()).replace(' ', 'T')[0:19], 'File creation date')

  # find channel numbers
  channum = np.arange(nbin)+1
  # define the columns for the FITS spectrum file
  col1 = pyfits.Column(name='CHANNEL', format='I', array=channum)
  col2 = pyfits.Column(name='COUNTS', format='J', array=s) # unit='count'
  col3 = pyfits.Column(name='QUALITY', format='I', array=qual)
  if group == False : # no grouping of bins
    cols = pyfits.ColDefs([col1, col2, col3]) # create a ColDefs (column-definitions) object
  else : # group the data to have at least group counts in each bin
    grpg = 0*s # array with grouping, +1 = start of new bin, -1 = continuing bin
    ig = imindub # start grouping at the lowest energy bin that is good
    grpg[ig] = 1 # mark start of new bin
    cg = s[ig] # number of counts in this bin
    # print ig, s[ig], cg, grpg[ig]
    while ig <= imaxdub : # move through bins until we reach the highest good bin
      ig += 1 # advance to next bin
      if cg >= group : # have enough counts, start a new
        grpg[ig] = 1 # start a new bin
        cg = s[ig] # number of counts in this bin
      else : # need more counts in this bin
        grpg[ig] = -1 # cotinue the bin
        cg += s[ig] # add to counts in the bin
      # print ig, s[ig], cg, grpg[ig]
    col4 = pyfits.Column(name='GROUPING', format='I', array=grpg)
    cols = pyfits.ColDefs([col1, col2, col3, col4]) # create a column-definitions object
  # print grpg
  # create a new binary table HDU object for the spectrum 
  tbhdu = pyfits.BinTableHDU.from_columns(cols)
  # add keywords to spectrum header for columns
  tbhdu.header['TTYPE1'] = ('CHANNEL', 'Channel number')
  tbhdu.header['TFORM1'] = ('I', 'data format of field')
  #tbhdu.header.remove('TUNIT1') # remove keyword in wrong place
  tbhdu.header.insert(10, ('TUNIT1', 'chan', 'physical unit of field'))
  #   tbhdu.header['TUNIT1'] = ('chan', 'physical unit of field')
  tbhdu.header.insert(11, ('TLMIN1', 1, 'minimum legal value'))
  tbhdu.header.insert(12, ('TLMAX1', nbin, 'maximum legal value'))
  tbhdu.header['TTYPE2'] = ('COUNTS', 'Total count in channel')
  tbhdu.header['TFORM2'] = ('J', 'data format of field')
  tbhdu.header.insert(15, ('TUNIT2', 'count', 'physical unit of field'))
  tbhdu.header['TTYPE3'] = ('QUALITY', 'Data quality')
  tbhdu.header['TFORM3'] = ('I', 'data format of field')
  # add keywords to spectrum header
  tbhdu.header['EXTNAME'] = ('SPECTRUM', 'Binary table extension name')
  tbhdu.header['HDUCLASS'] = ('OGIP', 'Format conforms to OGIP standards')
  tbhdu.header['HDUCLAS1'] = ('SPECTRUM', 'First class level')
  tbhdu.header['HDUCLAS2'] = ('TOTAL', 'Second class level')
  tbhdu.header['HDUCLAS3'] = ('COUNT', 'Third class level')
  tbhdu.header['HDUVERS'] = ('1.2.0', 'Version of format (OGIP memo OGIP-92-007')
  tbhdu.header['HDUVERS1'] = ('1.2.0', 'Obsolete included for back compatibility')
  tbhdu.header['AREASCAL'] = (1.0, 'Area scaling factor')
  tbhdu.header['BACKFILE'] = ('none', 'Associated background file')
  tbhdu.header['BACKSCAL'] = (1.0, 'Background file scale factor')
  tbhdu.header['CORRFILE'] = ('none', 'Associated scaling file')
  tbhdu.header['CORRSCAL'] = (1.0, 'Correction file scaling factor')
  tbhdu.header['RESPFILE'] = ('none', 'Associated redist matrix filename')
  #tbhdu.header['ANCRFILE'] = ('none', 'associated ancillary response filename')
  tbhdu.header['ANCRFILE'] = ('none', 'Associated ancillary response filename')
  tbhdu.header['PHAVERS'] = ('1992a', 'Obsolete')
  tbhdu.header['DETCHANS'] = (nbin, 'Total number possible channels')
  tbhdu.header['CHANTYPE'] = ('PI', 'Channel type (PHA, PI, etc.)')
  tbhdu.header['POISSERR'] = (True, 'Poissonian errors to be assumed')
  #tbhdu.header['POISSERR'] = ('T', 'Are Poisson Distribution errors assumed.')
  tbhdu.header['STAT_ERR'] = (0, 'No statistical error specified')
  tbhdu.header['SYS_ERR'] = (0, 'No systematic error specified')
  if group == False : # no grouping of bins
    tbhdu.header['GROUPING'] = (0, 'No grouping data has been specified')
  #tbhdu.header['QUALITY'] = (0, 'No data quality information specified')
  tbhdu.header['TELESCOP'] = ('HALOSAT', 'Telescope (mission) name')
  tbhdu.header['INSTRUME'] = ('SDD'+str(dpu), 'Instrument name')
  tbhdu.header['DATAMODE'] = ('PHOTON', 'Instrument datamode')
  tbhdu.header['FILTER'] = ('NONE', 'No filter was in use')
  tbhdu.header['OBSERVER'] = ('PHILIP KAARET', 'Principal Investigator')
  tbhdu.header['OBS_ID'] = (uf_evthdr['OBS_ID'], uf_evthdr.comments['OBS_ID']) # 'Observation ID'
  tbhdu.header['OBJECT'] = (uf_evthdr['OBJECT'], uf_evthdr.comments['OBJECT']) # 'Object/target name'
  tbhdu.header['OBJTYPE'] =(uf_evthdr['OBJTYPE'], uf_evthdr.comments['OBJTYPE']) # 'Object type observed'
  tbhdu.header['EQUINOX'] = (2000, '[yr] Equinox of celestial coord system')
  tbhdu.header['RADECSYS'] = ('FK5', 'Celestial coord system')
  tbhdu.header['RA_NOM'] = (uf_evthdr['RA_NOM'], uf_evthdr.comments['RA_NOM'])
  # (0.0, '[deg] R.A. of nominal aspect point [J2000]') #QQQ
  tbhdu.header['DEC_NOM'] = (uf_evthdr['DEC_NOM'], uf_evthdr.comments['DEC_NOM'])
  # (0.0, '[deg] Dec. of nominal aspect point [J2000]') #QQQ
  tbhdu.header['RA_OBJ'] = (uf_evthdr['RA_OBJ'], uf_evthdr.comments['RA_OBJ'])
  # (0.0, '[deg] Object Right Ascension [J2000]') #QQQ
  tbhdu.header['DEC_OBJ'] = (uf_evthdr['DEC_OBJ'], uf_evthdr.comments['DEC_OBJ'])
  # (0.0, '[deg] Object Declination [J2000]') #QQQ
  tbhdu.header['DATE-OBS'] = (str(getDate(min(gti[0])+32)).replace(' ', 'T')[0:19], \
  	                           'Start date of observations')
  tbhdu.header['DATE-END'] = (str(getDate(max(gti[1])+32)).replace(' ', 'T')[0:19], \
  	                           'End date of observations')
  tbhdu.header['TELAPSE'] = (max(gti[1])-min(gti[0]), '[s] Stop - Start')
  tbhdu.header['ONTIME'] = (ontime, '[s] Observation time on target')
  tbhdu.header['EXPOSURE'] = (ontime, '[s] Exposure')
  tbhdu.header['CLOCKAPP'] = (True, 'Clock correction applied? (F/T)')
  tbhdu.header['DEADAPP'] = (False, 'Has deadtime been applied to data? (F/T)')
  tbhdu.header['ORIGIN'] = ('UNIVERSITY OF IOWA', 'Origin of FITS file')
  tbhdu.header['TLM2FITS'] = (uf_evthdr['TLM2FITS'], uf_evthdr.comments['TLM2FITS'])
  # ('hkdb_20200124_scidb_20200124', 'Telemetry converter version number')
  tbhdu.header['PROCVER'] = (uf_evthdr['PROCVER'], uf_evthdr.comments['PROCVER'])
  # ('dohs_20200124', 'Processing script version number')
  tbhdu.header['SOFTVER'] = (uf_evthdr['SOFTVER'], uf_evthdr.comments['SOFTVER'])
  # ('searchdb_20200124', 'Software version')
  tbhdu.header['CREATOR'] = ('db_hsuf_hscl', 'Software creator of the file')
  tbhdu.header['CALDBVER'] = (uf_evthdr['CALDBVER'], uf_evthdr.comments['CALDBVER'])
  # ('hsYYYYMMDD', 'CALDB index version used')
  tbhdu.header['DATE'] = (str(datetime.now()).replace(' ', 'T')[0:19], \
  	                       'File creation date')
  # add CHECKSUM and DATASUM using ftools

  # if user requests, then write a keyword with the average rate in a PI band
  if ratekey :
    elow, ehigh = 3.0, 7.0 # energy band for rate
    pilow = max([1, 1+(elow-emin)/ebin])
    pihigh = min([1+nbin, 1+(ehigh-emin)/ebin])
    eband = str(elow)+'-'+str(ehigh)+' keV'
    # select events in band
    q = (pilow <= evt_pi) & (evt_pi <= pihigh)
    # average rate is rate = sum(q) / ontime
    # uncertainty is np.sqrt(sum(q))/ontime
    tbhdu.header['HARDRATE'] = (np.sum(q)/ontime, '[c/s] Average rate in '+eband)
    gerr = (np.sqrt(np.sum(q)+1.0)+1.5) # Gehrels error on number of counts
    tbhdu.header['HARDRERR'] = (gerr/ontime, '[c/s] Error on average rate in '+eband)
    # !!! QQQ do not do flux calculation until checked
    # avgflux = (np.sum((evt_pi[q]-1)*ebin_emin)/ontime
    # avgflux_err = avgflux*gerr/np.sum(q) # error is only from counts statistics, can we use energy?
    # tbhdu.header['HARDFLUX'] = (avgflux, '[keV/s] Average flux in '+eband)
    # tbhdu.header['HARDF_ERR'] = (avgflux/np.sqrt(1+np.sum(q)), 
    #                                '[keV/s] Error on average flux in '+eband)
  
  # join primary with extension
  thdulist = pyfits.HDUList([prihdu, tbhdu])
  print('Writing ', outroot+'_s'+str(dpu)+'.pi') 
  thdulist.writeto(outroot+'_s'+str(dpu)+'.pi', overwrite=True)

  # return the spectrum
  return s, ontime
