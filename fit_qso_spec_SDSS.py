##########################################
from astropy.io import fits
from scipy import interpolate
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
from sklearn import linear_model
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from astropy.stats import sigma_clip
from numpy.polynomial import Polynomial as P
from scipy.optimize import leastsq
from scipy.interpolate import BSpline, splrep, splev
from lmfit import minimize, Parameters, fit_report
from glob import glob
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as tck
from matplotlib.ticker import NullFormatter

AND = np.logical_and
OR = np.logical_or
####

def open_calibrate_fits(filename,path):
    hdu_raw = fits.open(str(path)+str(filename))
    loglam = hdu_raw[1].data['WAVE']   ##### Changed by Sapna
    flux = hdu_raw[1].data['FLUX']
    sig = hdu_raw[1].data['ERROR']
    c_sdss = ['CONTI_SDSS'][0]
    c_spline = ['CONTI_SPLINE']
    hdu_raw.close()
    return loglam, flux, sig,c_sdss,c_spline

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def running_median(datx,daty,daty_err,bin_size=21):
    xvals,yvals_w,yvals_unw = [],[],[]
    for j in range(len(datx)):
        if j+bin_size < len(datx):
            xvals.append( np.mean(datx[j:j+bin_size]) )
            yvals_w.append( np.average(daty[j:j+bin_size], weights=1.0/daty_err[j:j+bin_size]**2) )
            yvals_unw.append( np.median(daty[j:j+bin_size]) )
        elif j+bin_size >= len(datx):
            k = j
            bin_size = len(datx) - k
            xvals.append( np.mean(datx[j:j+bin_size]) )
            yvals_w.append( np.average(daty[j:j+bin_size], weights=1.0/daty_err[j:j+bin_size]**2) )
            yvals_unw.append( np.median(daty[j:j+bin_size]) )
            k = k-1
        # else:
        #     print('Hi3')
        #     xvals.append( datx[j])
        #     yvals_w.append( daty[j] )
        #     yvals_unw.append( daty[j])
        #     break
    return np.array(xvals),np.array(yvals_w),np.array(yvals_unw)

#############################
def residual_spline(params, x, y, y_err):      
    knot_i = params['knot_i'].value
    model = spline_fit(x,y,int(knot_i))
    return (y - model)/y_err


def spline_fit(ts,ys,ys_err,knot):
    n_interior_knots = knot
    qs = np.linspace(0, 1, n_interior_knots+2)[1:-1]
    knots = np.quantile(ts, qs)
    tck = splrep(ts, ys,w=1.0/ys_err**2, k=3, s=len(ts)-np.sqrt(2*len(ts)), task =0, t=knots) 
    ys_smooth = splev(ts, tck)
    return ys_smooth

def sigma_clip_box(datx,nsig_up,nsig_low , box):
    datx_s = np.array_split(np.array(datx),box)
    #print(len(datx))
    id_clip = []
    t_end = 0
    for ii in range(box):
        
        if ii == 0: t_str,t_end =0, len(datx_s[ii])-1
        else: t_str,t_end = t_end + 1, t_end + len(datx_s[ii])

        id_datx = np.arange(t_str,t_end+1,1)
        #print(ii,t_str,t_end, np.shape(datx[id_datx]),len(datx_s[ii]))
        #print('All',id_datx)
        sigclip = sigma_clip(datx[id_datx], sigma_lower = nsig_low,sigma_upper =nsig_up, maxiters=10,masked=False,axis=0)
        id_nan = np.where( np.isnan(sigclip) == True)
        id_nonan = np.where( np.isnan(sigclip) == False)
        #print('Nonan',id_datx[id_nonan])
        #print('Nan',id_datx[id_nan])        
        id_clip =  np.append(id_clip, id_datx[id_nonan] )
    return id_clip.astype(int)
    
################################################
#rchi_sdss,rchi_new = [],[]
#print(ff,file_in[ff],'Reduced chi-sq b/w clipped and new', sum( (flux_clipped - flux_spline)**2/sig[id_mgii]**2) / len(flux_spline) )
#print(ff,'Reduced chi-sq b/w clipped and SDSS', sum( (flux_clipped - sdss_conti[id_mgii])**2/sig[id_mgii]**2) / len(flux_spline) )    
#rchi_sdss.append(sum( (flux_clipped - sdss_conti[id_mgii])**2/sig[id_mgii]**2) / len(flux_spline))
#rchi_new.append(sum( (flux_clipped - flux_spline)**2/sig[id_mgii]**2) / len(flux_spline))



def main_fit_part(runid,file_in,zemi_q,z_cl,path):

    [wav, flux, sig, sdss_conti,spline_conti] = open_calibrate_fits(filename=file_in,path=path)   ### for getting restframe spectra
    
    w_obs_mgii = (1.0 + z_cl)* 2796.35
    vel = 299792.46*(wav/w_obs_mgii - 1)     
    vrange = 10000.0      
    #lambda_low = (1+zemi_q)*1250.00     ###w_obs_mgii*(-vrange/299792.46 + 1.0)   
    #lambda_up =  (1+zemi_q)*2796.35     ###w_obs_mgii*(vrange/299792.46 + 1.0)     


    lambda_low  = (1.0- 30000.0/299792.46)*(1.0 + 0.32320)*2796.35
    lambda_up   = (1.0+ 30000.0/299792.46)*(1.0 + 0.74980)*2803.53

    id_mgii =    np.where(AND(wav > lambda_low, wav < lambda_up))[0]  

    snr_mgii = np.median(flux[id_mgii]/sig[id_mgii])
    
    ###################### To calculate emission line 
    vprox = 3000.0
    w_ion = [1215.6701,1398.266545, 1549.49265,1908.734 ,2799.9429]
    w_ion_s = ['Lya','SiIV','CIV','CIII','MgII']  
    em_arr = []
    for kk in range(len(w_ion)):
        #if kk == 0: vprox = 10000.0        
        C1 = lambda_low < (1.0 + zemi_q)*w_ion[kk] and  lambda_up > (1.0 + zemi_q)*w_ion[kk]
        C2 = lambda_low < (1.0- vprox/299792.46)*(1.0 + zemi_q)*w_ion[kk] and  lambda_up > (1.0 - vprox/299792.46)*(1.0 + zemi_q)*w_ion[kk]
        C3 = lambda_low < (1.0 + vprox/299792.46)*(1.0 + zemi_q)*w_ion[kk] and  lambda_up > (1.0 + vprox/299792.46)*(1.0 + zemi_q)*w_ion[kk]
        if C1 or C2 or C3: em_arr = np.append(em_arr,'Emission-In')
        else: em_arr = np.append(em_arr,'Emission-Out')
        #print(lambda_low,lambda_up,(1.0- vprox/299792.46)*(1.0 + zemi_q)*w_ion[kk],(1.0+ vprox/299792.46)*(1.0 + zemi_q)*w_ion[kk],(1.0 + zemi_q)*w_ion[kk])        
    id_em = np.where(em_arr == 'Emission-In')[0]
    if np.size(id_em) > 0: em_tag = 'Emission-In'
    else: em_tag = 'Emission-Out'

    SaveFile = open('tmp4_em', 'a')
    SaveFile.write(Fmt % ( file_in,  np.size(id_mgii), np.median(flux[id_mgii]/sig[id_mgii]), np.median(flux[id_mgii]/np.std(flux[id_mgii])), lambda_low, lambda_up, em_tag ))
    SaveFile.close()

#####################################  main input part of the code #####################
nrow = 1
path = '/home/sapnamisra/GalaxyQSO/Proj2_SDSS/PCA-conti-mod/PCA-Spline-Conti-full-2sig-51pix/'


#pair_file = 'main_1-5r500_sample-4-NoBAL.txt'
pair_file = 'main_1-5r500_sample-4-NoGap-NoBAL.txt'
ra_q,dec_q,zemi_q,SNR,Mi,imag,ra_cl,dec_cl,z_cl,M500,R500,q_cl_vel,rhy,rhy_R500,RL500,N500sp,N500,R500_wen = np.loadtxt(pair_file, usecols=[0,1,2,4,5,6,8,9,10,11,12,13,14,15,17,18,19,20],unpack=True,skiprows=1)
pairid,gal_flg_spec,gal_flghot = np.loadtxt(pair_file, usecols=[22,23,24],unpack=True,skiprows=1,dtype=int)
name_q,file_in,cat,name_cl = np.loadtxt(pair_file, usecols=[3,7,16,21],unpack=True,skiprows=1,dtype='str')



SaveFile = open('tmp4_em', 'a')
SaveFile.write('File pix snr_er snr_f lam_low lam_up em-info\n')
SaveFile.close()
Fmt = "%s  %i %10.5f %10.5f  %10.5f  %10.5f %s\n"

for ff in range(len(file_in)):
#for ff in range(0,5):
    print(file_in[ff])
    main_fit_part(ff,file_in[ff],zemi_q[ff],z_cl[ff],path)

