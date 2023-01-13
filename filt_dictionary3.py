# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 15:56:17 2019
filt_dictionary.py
@author: rringuette
Collected definitions of the filt dictionary

Updates Dec 18 2019:
-	Remove extraneous code
-	Update lcpha definition based on wiki page
-	Change spacing of filt definition for better output
-	Update filt definition based on wiki page

#need alternate lcpha for dark earth clean cuts than elsewhere
def lcdef(lcpha_type=''):
lcpha = fd.lcdef(filter_type='cl_de2') #for dark earth calibations
lcpha = fd.lcdef() #for default

#current uf and cl (mapping) versions of cuts on 8 second bins
def filt_def(cut_note):
filt = filt_def(cut_note)  

#JKB 10/13/2022 - this code was crashing due to UTF 8 error cant decode byte 0x91 in posit 177 - looked like corrupt character in comments? removed characters.
"""

#need alternate lcpha for dark earth clean cuts than elsewhere
def lcdef(lcpha_type=''):
    
    #use as default for all so light curves are written to evt file
    lcpha = [['LC_SCI', 900, 2250], # events in ~0.45 to ~3 keV band
         ['LC_FILT', 600, 12900], #events to be written to file, for code testing
         ['LC_GAP1', 600, 750], #first half of gap between electronic noise band and O band
         ['LC_HARD', 5600, 12900], # ~3 keV to ~7 keV
         ['LC_VLE', 12900, 16384], # ~7 keV to max amplitude
         ['LC_ALL', 0, 16384],   #entire energy range
         ['LC_OXYGEN', 900, 1400], # band with oxygen lines 0.45 to 0.73        
         ['LC_RESET', 0, 100],  #band with reset pulses
         ['LC_ALSI', 2250, 3750], # 1.2 keV to 2.0 keV (Al and Si lines +/- 3*90eV)             
         ['LC_UP', 1400, 2250]]  #0.73 to 1.2 keV       

    return lcpha

#current uf and cl versions of cuts on 8 second bins
def filt_def(cut_note):
    
    filt_dict={}
    #version posted on wiki without background and noise removal 
    filt_dict['inst']='(hk["IN_SAA"] == False) \
          & (hk["SDD_TEMP"] >= -31.0) & (hk["SDD_TEMP"] <= -29.0) \
          & (hk["MON_3P3V"] >= 3.20) & (hk["MON_3P3V"] <= 3.30) \
          & (hk["MON_P5V"] >= 4.80) & (hk["MON_P5V"] <= 5.10) \
          & (hk["MON_M5V"] >= -5.10) & (hk["MON_M5V"] <= -4.80) \
          & (hk["SDDHVMON"] >= -138.0) & (hk["SDDHVMON"] <= -133.0) \
          & (hk["SDD0"] >= 0.060) & (hk["SDD0"] <= 0.080)'
    filt_dict['Base']='(hk["IN_SAA"] == False) \
          & (hk["SDD_TEMP"] >= -31.0) & (hk["SDD_TEMP"] <= -29.0) \
          & (hk["MON_3P3V"] >= 3.20) & (hk["MON_3P3V"] <= 3.30) \
          & (hk["MON_P5V"] >= 4.80) & (hk["MON_P5V"] <= 5.10) \
          & (hk["MON_M5V"] >= -5.10) & (hk["MON_M5V"] <= -4.80) \
          & (hk["SDDHVMON"] >= -138.0) & (hk["SDDHVMON"] <= -133.0) \
          & (hk["SDD0"] >= 0.060) & (hk["SDD0"] <= 0.080) \
          & (hk["LC_ALL"] <= 140) & (hk["LC_GAP1"] <= 1.0) & (hk["TIME_DUR"] >=8.0)'
    #LC_ALL>15 added to uf definition 20190924. Was previously LC_ALL>0
    filt_dict['uf']= filt_dict['Base']+' & (hk["OFFSET"] <= 0.25) \
          & (hk["NADIR_ANGLE"] >= 92.0) & (hk["NADIR_ANGLE"] <= 181.0) \
          & (hk["LC_ALL"] > 15)'
    #standard with dark earth nadir angle cutcut_note='cl'
    filt_dict['uf_de']=filt_dict['Base']+' & (hk["OFFSET"] <= 0.25) \
          & (hk["NADIR_ANGLE"] < 65.0) & (hk["LC_ALL"] > 0)'       
    filt_dict['cl'] = filt_dict['uf']  #use uf and add cuts in filter_pars
    #It looks like requiring the nadir angle to be less than 65 and the 
    #LC_AlSi/LC_UP ratio to be greater than 1.5 is doing the best job for dark earth calibration spectra    
    filt_dict['cl_de2']=filt_dict['uf_de']+' & (hk["LC_ALSI"]/hk["LC_UP"] > 1.5)'
    filt_dict['cl_de']=filt_dict['uf_de']
    filt_dict['cl_halo']=filt_dict['uf']+' & (hk["LC_HARD"] <= 0.25) & (hk["LC_VLE"] <= 2.0)'
    filt_dict['cl_sci']=filt_dict['uf']+' & (hk["LC_HARD"] <= 0.75) & (hk["LC_VLE"] <= 2.0)'

    return filt_dict[cut_note] 
        
