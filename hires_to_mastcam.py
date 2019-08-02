#PROGRAM NAME: hires_to_mastcam_wbayer.pro
#WRITTEN BY: Ryan Anderson, Melissa Rice, Kathleen Hoza
#PURPOSE: to read in high resolution spectra and 
#convolve them to Mastcam bandpasses
# 
# Adapted to Mastcam from Melissa Rice's original hires_to_pancam_S.pro program: 10/26/2010 - Ryan Anderson
# 11/1/2010 - Changed effective wavelengths to centroid of the filter. - Ryan Anderson
# 11/22/2010 - Changed to use rd_tfile.pro to read in data. Added delimeter keyword. - Ryan Anderson
# 12/22/2010 - added microns keyword - Ryan Anderson
# 11/15/2011 - adapted to the new Mastcam passbands provided by Jim Bell
# 3/15/2013 - will now output convolved RGB Bayer Filter reflectances
# 6/18/2019 - translated to Python from IDL - Kathleen Hoza
# INPUTS: input_file (the ASCII file containing the spectra to convolve.
#                     The first column must be wavelength in nm, and
#                     the following columns must contain reflectance)
#          output_file (the program will save the convolved 
#                      wavelengths and reflectances to the given
#                      filename (must end with '.txt'))
#          ncol (the number of columns in the input file - can be
#                approximated to the nearest power of 10)
# 		  microns - set this keyword if the input spectra have wavelengths in microns. 
# 					Wavelengths will be multiplied by 1000 to convert to nm
# 
# OUTPUTS: a file is written with mastcam wavelengths in the first column and the convolved
#          mastcam spectra in the following columns.
# 
# 
from scipy import interpolate
from scipy import integrate
import numpy as np
import os


#READ IN THE HIGH-RESOLUTION SPECTRA AND WAVELENGTHS
microns=False
input_file='Z:\Kristiana\Spectra\MCZ_targets\Full_redo_run\spectra\full_run.csv'
os.chdir('C:\Users\krist\OneDrive\Documents\Python\mastcam_convolve\MCZ_probe')
input_data=np.genfromtxt(input_file, dtype=float, delimiter=',',skip_header=5, unpack=True)
spectrum_wav=np.array(input_data[0])

if microns:
     spectrum_wav=spectrum_wav*1000
#spectrum_wav=spectrum_wav(where(spectrum_wav gt 0))

n=len(input_data) #number of spectra we are convolving+1 for wavelengths

#DEFINE AN ARRAY FOR THE OUTPUTS
output_data=[[]]

for m in range(1,n):
    #spectrum_orig=reform(double(input_data[m,*]))
    #spectrum_orig=spectrum_orig(where(spectrum_wav gt 0))
    spectrum_orig=np.array(input_data[m])
    # spectrum_orig=[1,-1,2,-2]
    # indices_less_zero=[(i,r) for (i,r) in zip(spectrum_orig,x) if i >= 0]
    # print(indices_less_than_zero)
    
    #READ IN THE SOLAR SPECTRUM AND PUT INTO INPUT SPECTRUM WAVELENGTHS

    solar_input=np.genfromtxt('sun_input.txt', unpack=True)

    interpol = interpolate.interp1d(solar_input[0],solar_input[1], kind='linear')
    solar_new=interpol(spectrum_wav) #interpolate values of solar function at spectrum wavelengths
    
    #multiply high-res spectrum by the interpolated solar spectrum
    spectrum=spectrum_orig*solar_new
    
    #fix any crazy negative values in the spectrum - written in IDL, not currently implemented in python
    #zero=where(spectrum LT 0.)
    #IF zero[0] GT 0 THEN spectrum(where(spectrum LT 0.))=0.
    
    #define wavelength array (300-1100 nm with 5nm steps) to use with Mastcam filter data
    wvl=np.arange(350,1040,5)#make_array(161,/index)*5.+300.
    
    #restore merged / normalized Mastcam filter data
    [bayer_wvl, bayerb, bayerg, bayerr]=np.genfromtxt('mastcam/filter0_bayer.txt',unpack=True)
    [filter1_525nm_wvl,filter1_525nm]=np.genfromtxt('mastcam/filter1_525nm_34mm.txt',unpack=True)
    [filter2_440nm_wvl,filter2_440nm]=np.genfromtxt('mastcam/filter2_440nm_34mm.txt',unpack=True)
    [filter3_750nm_wvl,filter3_750nm]=np.genfromtxt('mastcam/filter3_750nm_34mm.txt',unpack=True)
    [filter3_800nm_wvl,filter3_800nm]=np.genfromtxt('mastcam/filter3_800nm_100mm.txt',unpack=True)
    [filter4_675nm_wvl,filter4_675nm]=np.genfromtxt('mastcam/filter4_675nm_34mm.txt',unpack=True)
    [filter4_905nm_wvl,filter4_905nm]=np.genfromtxt('mastcam/filter4_905nm_100mm.txt',unpack=True)
    [filter5_865nm_wvl,filter5_865nm]=np.genfromtxt('mastcam/filter5_865nm_34mm.txt',unpack=True)
    [filter5_935nm_wvl,filter5_935nm]=np.genfromtxt('mastcam/filter5_935nm_100mm.txt',unpack=True)
    [filter6_1035nm_wvl,filter6_1035nm]=np.genfromtxt('mastcam/filter6_1035nm_34mm.txt',unpack=True)
    
    
    #interpolate the the narrow-band filters and solar spectrum to the same wavelength scale
    interpol = interpolate.interp1d(bayer_wvl,bayerr, kind='linear')
    filter0_red=interpol(wvl)
    interpol = interpolate.interp1d(bayer_wvl,bayerg, kind='linear')
    filter0_green=interpol(wvl)
    interpol = interpolate.interp1d(bayer_wvl,bayerb, kind='linear')
    filter0_blue=interpol(wvl)
    
    interpol = interpolate.interp1d(filter1_525nm_wvl,filter1_525nm, kind='linear')
    filter1_525nm=interpol(wvl)
    
    interpol = interpolate.interp1d(filter2_440nm_wvl,filter2_440nm, kind='linear')
    filter2_440nm=interpol(wvl)
    
    interpol = interpolate.interp1d(filter3_750nm_wvl,filter3_750nm, kind='linear')
    filter3_750nm=interpol(wvl)
    
    interpol = interpolate.interp1d(filter3_800nm_wvl,filter3_800nm, kind='linear')
    filter3_800nm=interpol(wvl)
    
    interpol = interpolate.interp1d(filter4_675nm_wvl,filter4_675nm, kind='linear')
    filter4_675nm=interpol(wvl)
    interpol = interpolate.interp1d(filter4_905nm_wvl,filter4_905nm, kind='linear')
    filter4_905nm=interpol(wvl)
    
    interpol = interpolate.interp1d(filter5_865nm_wvl,filter5_865nm, kind='linear')
    filter5_865nm=interpol(wvl)
    interpol = interpolate.interp1d(filter5_935nm_wvl,filter5_935nm, kind='linear')
    filter5_935nm=interpol(wvl)
    
    interpol = interpolate.interp1d(filter6_1035nm_wvl,filter6_1035nm, kind='linear')
    filter6_1035nm=interpol(wvl)

    interpol=interpolate.interp1d(spectrum_wav,solar_new,kind='linear')
    solar_interp=interpol(wvl)
    
    #Integrate over each filter
    filter0_red_int=integrate.trapz(filter0_red,wvl)
    filter0_green_int=integrate.trapz(filter0_green,wvl)
    filter0_blue_int=integrate.trapz(filter0_blue,wvl)
    filter1_525nm_int=integrate.trapz(filter1_525nm,wvl)
    filter2_440nm_int=integrate.trapz(filter2_440nm,wvl)
    filter3_750nm_int=integrate.trapz(filter3_750nm,wvl)
    filter3_800nm_int=integrate.trapz(filter3_800nm,wvl)
    filter4_675nm_int=integrate.trapz(filter4_675nm,wvl)
    filter4_905nm_int=integrate.trapz(filter4_905nm,wvl)
    filter5_865nm_int=integrate.trapz(filter5_865nm,wvl)
    filter5_935nm_int=integrate.trapz(filter5_935nm,wvl)
    filter6_1035nm_int=integrate.trapz(filter6_1035nm,wvl)
    
    #normalize each filter so that the integral is 1
    filter0_red=filter0_red/filter0_red_int
    filter0_green=filter0_green/filter0_green_int
    filter0_blue=filter0_blue/filter0_blue_int
    filter1_525nm=filter1_525nm/filter1_525nm_int
    filter2_440nm=filter2_440nm/filter2_440nm_int
    filter3_750nm=filter3_750nm/filter3_750nm_int
    filter3_800nm=filter3_800nm/filter3_800nm_int
    filter4_675nm=filter4_675nm/filter4_675nm_int
    filter4_905nm=filter4_905nm/filter4_905nm_int
    filter5_865nm=filter5_865nm/filter5_865nm_int
    filter5_935nm=filter5_935nm/filter5_935nm_int
    filter6_1035nm=filter6_1035nm/filter6_1035nm_int
    

    #interpolate the high-res spectrum to the same wavelength scale as the narrow filters
    interpol=interpolate.interp1d(spectrum_wav,spectrum,kind='linear')
    new_spec_rad = interpol(wvl)
    
    
    #multiply high-res spectrum by each filter, integrate, and divide out solar spectrum
    filter0_red=integrate.trapz(filter0_red*new_spec_rad, wvl)/integrate.trapz(filter0_red*solar_interp,wvl)

    filter0_green=integrate.trapz(filter0_green*new_spec_rad, wvl)/integrate.trapz(filter0_green*solar_interp,wvl)
           
    filter0_blue=integrate.trapz(filter0_blue*new_spec_rad, wvl)/integrate.trapz(filter0_blue*solar_interp,wvl)

    filter1_525nm=integrate.trapz(filter1_525nm*new_spec_rad, wvl)/integrate.trapz(filter1_525nm*solar_interp,wvl)

    filter2_440nm=integrate.trapz(filter2_440nm*new_spec_rad, wvl)/integrate.trapz(filter2_440nm*solar_interp,wvl)

    filter3_750nm=integrate.trapz(filter3_750nm*new_spec_rad, wvl)/integrate.trapz(filter3_750nm*solar_interp,wvl)

    filter3_800nm=integrate.trapz(filter3_800nm*new_spec_rad, wvl)/integrate.trapz(filter3_800nm*solar_interp,wvl)

    filter4_675nm=integrate.trapz(filter4_675nm*new_spec_rad, wvl)/integrate.trapz(filter4_675nm*solar_interp,wvl)

    filter4_905nm=integrate.trapz(filter4_905nm*new_spec_rad, wvl)/integrate.trapz(filter4_905nm*solar_interp,wvl)

    filter5_865nm=integrate.trapz(filter5_865nm*new_spec_rad, wvl)/integrate.trapz(filter5_865nm*solar_interp,wvl)

    filter5_935nm=integrate.trapz(filter5_935nm*new_spec_rad, wvl)/integrate.trapz(filter5_935nm*solar_interp,wvl)

    filter6_1035nm=integrate.trapz(filter6_1035nm*new_spec_rad, wvl)/integrate.trapz(filter6_1035nm*solar_interp,wvl)

    #define arrays
    mastcam_spect=[filter0_red,filter0_green,filter0_blue,filter1_525nm,filter2_440nm,filter3_750nm,filter3_800nm,filter4_675nm,filter4_905nm,filter5_865nm,filter5_935nm,filter6_1035nm]
    
    #define mastcam effective wavelengths 
    mastcam_wav = [640.,554.,495.,527.,445.,751.,805.,676.,908.,867.,937.,1012.]
    
    #sort the mastcam arrays in ascending order
    mastcam=[x for _, x in sorted(zip(mastcam_wav,mastcam_spect), key=lambda pair: pair[0])]
    
    output_data.append(mastcam)

mastcam_wav=sorted(mastcam_wav)
output_data[0]=mastcam_wav
headers=[]
with open(input_file,'r') as f:
     for i in range(5):
          headers.append(f.readline())
print(os.getcwd())
with open(input_file.strip('.csv')+'_convolved.csv','w+') as f:
     for line in headers:
          f.write(line)
     for i in range(len(output_data[0])):
          for j in range(len(output_data)):
               f.write(str(output_data[j][i])+',')
          f.write('\n')

