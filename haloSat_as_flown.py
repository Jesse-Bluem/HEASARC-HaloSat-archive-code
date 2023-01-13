#Update JKB 10/10/2022 - changed to Py3. Changed file write type 'wb' to 'w', as Py3 requires wb to use byte objects not strings.

import astropy.io.fits as pyfits
import sys
import os
import datetime
import pandas as pd

def getDate(time):
  epoch_date = pd.to_datetime('2000/01/01 00:00:00')
  time1 = time - 37.0
  date = epoch_date + datetime.timedelta(seconds=time1)
  return date

def getTime(date):
  epoch_date = pd.to_datetime('2000/01/01 00:00:00')
  date1 = pd.to_datetime(date)
  time = (date1-epoch_date).total_seconds()
  time += 37
  return time

df_cols=['target_name','target_id','type','start_time', 'end_time', 'ra', 'dec','start_date','end_date']
def makeRow(target_name, target_id,target_type,start_time, end_time, ra, dec,start_date,end_date):
  d = {'target_name': [target_name], 'target_id': [target_id],'type':[target_type], 'start_time': [start_time],
      'end_time': [end_time], 'ra': [ra], 'dec': [dec],'start_date':[start_date],'end_date':[end_date]}
  d = pd.DataFrame(data=d, columns=df_cols)
  return(d)


def getFiles(directory='/home/halosat/Blt_n30/'):
  directory='{0}{1}'.format(directory,'v1_local/')
  files=[]
  for x in os.walk(directory):
    sub_directory = x[0]
 
    if 'unfiltered' in sub_directory:
      try:
        source_id=sub_directory[-17:-11]
        print(source_id[-2:])
        assert(source_id[-2:]=='01')
        print(len(sub_directory),source_id,sub_directory)
        files.append('{0}/hs{1}.att'.format(sub_directory,source_id))
      except:
        print('***Error in getFiles ',sub_directory)
  return files
#return ['/home/halosat/Blt_n30/v1_local/000801/unfiltered/hs000801.att',
  
def makeDf(directory):
  filelist = getFiles(directory)
  df = pd.DataFrame(columns=df_cols)
  for afile in filelist:
    print(afile)
    with pyfits.open(afile) as hdulist:
      print(hdulist['GTI'].header['DATE-END'])
      target_name = hdulist['GTI'].header['OBJECT']
      target_id = hdulist['GTI'].header['OBS_ID']
      target_type = hdulist['GTI'].header['OBJTYPE']
      ra =  hdulist['GTI'].header['RA_NOM']
      dec =  hdulist['GTI'].header['DEC_NOM']
      a = hdulist['GTI'].data
      for i in range(len(a)):
        start,end = a[i]['START'],a[i]['STOP']
        duration = a[i]['TOTAL_TIME']
        start_date,end_date = a[i]['START_DATE'],a[i]['END_DATE']
        if duration > 30:
          d=makeRow(target_name, target_id,target_type,start,end, round(ra,3), round(dec,3),start_date,end_date)
          df=pd.concat([df,d],ignore_index=True)
          #print start,end,duration
          #print df
  return df
    

def makeFits(df):
  prihdr = pyfits.Header()
  prihdr['TELESCOP'] = 'HaloSat'
  prihdr['DATE'] = (str(datetime.datetime.now()), 'FITS file creation time')
  prihdr['DATE-START'] = (str(min(df.start_time.values))[0:19], 'Earliest time for as flown data')
  prihdr['DATE-END'] = (str(max(df.end_time.values))[0:19], 'Latest time for as flown data')
  #prihdr['CHECKSUM'] = 
  #prihdr['DATASUM'] = 
  prihdu = pyfits.PrimaryHDU(header=prihdr)

  df.loc[:,'start_time'] = df.start_time.apply(lambda x: str(x)[0:19])
  df.loc[:,'end_time'] = df.start_time.apply(lambda x: str(x)[0:19])

  col1  = pyfits.Column(name='OBS_ID',                format='6A',  array=df.target_id.values)
  col2  = pyfits.Column(name='OBJECT',         format='22A', array=df.target_name.values)
  col3  = pyfits.Column(name='OBJTYPE',         format='11A', array=df.type.values)
  col4  = pyfits.Column(name='DATE-OBS',          format='19A', array=df.start_time.values,unit='UTC')
  col5  = pyfits.Column(name='DATE-END',            format='19A', array=df.end_time.values,unit='UTC')
  col6  = pyfits.Column(name='RA',                  format='E',   array=df.ra.values,unit='degrees')
  col7  = pyfits.Column(name='DEC',                 format='E',   array=df.dec.values,unit='degrees')

  cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6,col7])
  tbhdu = pyfits.BinTableHDU.from_columns(cols)

  tbhdu.header.comments['TTYPE1'] = 'Observation ID'
  tbhdu.header.comments['TTYPE2'] = 'Target Name'
  tbhdu.header.comments['TTYPE3'] = 'Object Type'
  tbhdu.header.comments['TTYPE4'] = 'Start date of observation (UTC)'
  tbhdu.header.comments['TTYPE5'] = 'End date of observation (UTC)'
  tbhdu.header.comments['TTYPE6'] = 'Object Right ascension [J2000]'
  tbhdu.header.comments['TTYPE7'] = 'Object Declination [J2000]'

  #write header and table to FITS file
  thdulist = pyfits.HDUList([prihdu, tbhdu])
  thdulist.writeto('haloSat_asflown.fits', clobber=True)

def makeTxt(df,directory):
  directory = '{0}{1}'.format(directory,'v1/')
  def get_UT(tai_time):
    ut_time = getDate(tai_time)
    return ut_time.strftime("%Y-%m-%dT%H:%M:%S")

  with open('{0}haloSat_as_flown.txt'.format(directory),'w') as file: #JKB, was 'wb', b = bytes requires special format in Py3.
  
    file.write('<HEADER>\n')
    file.write('#               TABLE: heasarc_QQQ\n')
    file.write('#    TOTAL ROWS: QQQ\n')
    file.write('# table_description = "HaloSat As Flown"\n')
    file.write('#\n')
    file.write('# Table Parameters\n')
    file.write('#\n')
    file.write('# type[:fmt][_unit] [[ucd]] [(index)|(key)] // description [// comment]\n')
    file.write('# field[date_obs] = char20  // GTI Date Start \n')
    file.write('# field[date_end] = char20  // GTI Date End \n')
    file.write('# field[ra] = float8:.3f_degree  // Right Ascension \n')
    file.write('# field[dec] = float8:.3f_degree  // Declination \n')
    file.write('# field[source_id] = char6 // Source ID \n')
    file.write('# field[objtype] = char11 // Source Type \n')
    file.write('# field[object] = char24   // Source Name \n')
    file.write('# \n')
    file.write('# Data Format Specification \n')
    file.write('# \n')
    file.write('#line[1] = date_obs date_end ra dec source_id objtype object\n')
    file.write('# \n')
    file.write('<DATA>\n')

    for i in df.index:
      #row = '{0}|{1}|{2:>7}|{3:>7}|{4}|{5:11}|{6:24}\n'.format(get_UT(df.loc[i,'start_time']),get_UT(df.loc[i,'end_time']),'{0:.3f}'.format(df.loc[i,'ra']),'{0:.3f}'.format(df.loc[i,'dec']),df.loc[i,'target_id'],df.loc[i,'type'],df.loc[i,'target_name'])
      row = '{0}|{1}|{2:>7}|{3:>7}|{4}|{5:11}|{6:24}\n'.format(df.loc[i,'start_date'],df.loc[i,'end_date'],'{0:.3f}'.format(df.loc[i,'ra']),'{0:.3f}'.format(df.loc[i,'dec']),df.loc[i,'target_id'],df.loc[i,'type'],df.loc[i,'target_name'])

      print(row)
      file.write(row)
    file.write('#<END>\n')

if __name__=='__main__':
  directory = sys.argv[1]
  print('\nParent Directory:{0}'.format(directory))

  if True:
    df = makeDf(directory)
    df.sort_values(by=['start_time'],inplace=True)
    df.reset_index(inplace=True)
    makeTxt(df,sys.argv[1])
  else:
    print('Quitting')
  

  



