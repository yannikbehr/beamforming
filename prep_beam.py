#!/usr/bin/env mypython

"""
calculate spatial correlation matrix for subsequent use
in beamformer
"""
import os
import sys
import glob
import obspy.sac
from obspy.sac import *
import obspy.signal
import numpy as np
import scipy.signal
import scipy.io
from pylab import plot, show

class BeamFormException(Exception): pass

year = 2001
#sacpath = '/home/data/dev/beamforming/geopsy/dataII'
#sacpath = '/data/wanakaII/yannik/cnipse/sacfiles/2001/Apr/2001_4_30_0_0_0/'
#sacpathII = '/data/wanakaII/yannik/start/sacfiles/2001/Apr/2001_4_30_0_0_0/'
sacpath = '/data/wanakaII/yannik/start/sacfiles/2001/Mar/2001_3_3_0_0_0/'
matpath = './Matfiles_start/' 
JD = '62'
component = 'HHZ'
allstations = glob.glob(os.path.join(sacpath,'ft_grid*.HHZ.SAC'))
#allstations += glob.glob(os.path.join(sacpathII,'ft_grid*.*HZ.SAC'))
if len(allstations) < 1:
    raise BeamFormException("list of input files is empty")
STATIONS = []
Ista = 0

### length of time vector
fftpower = 7
dt = 1
Fsamp = 10
Fmax=1/dt
Nhours=3 #add up ffts over a 3 hour period
Nhours=1 #ffting every hour
Ntimes=int(round(24/Nhours))
times=np.arange(0,24+Nhours,Nhours)
N=Nhours*3600*Fmax
step=N
Fs=Fmax/N #no. of samples all up
TEMP_FILTER='y' #'y' means we want to set a threshold for temporal filtering
Threshold_STD=0.5 #set threshold to be half of mean sigma
freq = Fmax/2.*np.arange(0,2**(fftpower-1))/2**(fftpower-1)
freq_int = [0.02,0.4]
                     
### initialise station no.s
ista = 0
stations = []
INFOM = []
infom = []
info = {}
nfiles = len(allstations)

for _f in allstations:
    sacalltrace = np.zeros((1,24*3600*Fsamp))
    print _f
    x = ReadSac(_f)
    info['staname'] = x.GetHvalue('kstnm').rstrip()
    info['kcomp'] = x.GetHvalue('kcmpnm').rstrip()
    info['slat'] = x.GetHvalue('stla')
    info['slon'] = x.GetHvalue('stlo')
    fs = round(1./x.GetHvalue('delta'))
    info['fs'] = fs
    if info['fs'] == Fsamp:
        if 1000 < x.GetHvalue('npts'):
            x.seis -= x.seis.mean()
            istart = int(round(((x.GetHvalue('nzhour')*60.+x.GetHvalue('nzmin'))*60+x.GetHvalue('nzsec'))*Fsamp))
            sactrace = np.zeros((1,24*3600*Fsamp))
            x.seis = x.seis.reshape(1,x.seis.size)
            sactrace[0,istart:(istart+x.GetHvalue('npts'))]=x.seis[0,:]
            sacalltrace = np.append(sacalltrace,sactrace,axis=0)
            del sactrace
    x.seis = sacalltrace[1]
    del sacalltrace
    if info['staname'] not in stations:
        stations.append(info['staname'])
        infom.append(info)
        istacur=ista
        ista += 1
    if info['staname'] not in STATIONS:
        Ista += 1
        STATIONS.append(info['staname'])
        INFOM.append(info)
        
    x.seis -= x.seis.mean()
    istart = 0
    npts = x.seis.size
    seis0 = np.zeros((1,24*3600*fs))
    seis0[istart:(istart+npts)] = x.seis
    Step = step*fs
    seisband = np.ndarray(shape=(step,Ntimes))
    for ii in range(Ntimes):
        Ilow = (ii)*Step
        Iup = (ii+1)*Step
        seis = seis0[0,Ilow:Iup]
        seisband1 = obspy.signal.filter.bandpassZPHSH(seis,freq_int[0],freq_int[1],df=fs)
        seisband[:,ii] = scipy.signal.resample(seisband1,len(seisband1)/fs)
        seisband[:,ii] -= seisband[:,ii].mean()

    sigmas = np.std(seisband,axis=0)
    I = np.where(sigmas==0)
    Iok = np.where(sigmas>0)
    sigmas[Iok].sort()
    sigmasort = sigmas
    sigma = np.sqrt(sum(sigmasort**2)/len(sigmasort))
    iday = x.GetHvalue('nzjday')

    ###temporal filtering (set a threshold)
    seis = np.ndarray(shape=(step,Ntimes))

    for ii in range(Ntimes):
        seisz = seisband[:,ii]
        if TEMP_FILTER == 'y':
            Threshold_Temp = Threshold_STD*sigma #Arbitrary
            II = np.where(abs(seisz)>Threshold_Temp)
            seisz[II] = Threshold_Temp*np.sign(seisz[II])
            seisz -= seisz.mean()
        seis[:,ii] = seisz

    ### find indices of frequencies that lie within our interval
    frq = freq
    I = np.where((frq>freq_int[0]) & (frq<freq_int[1]))

    Ismall = 2**fftpower
    ipick = np.arange(Ismall)
    Nsub = int(np.floor(N/Ismall)) # Number of time pieces -20 mins long each
    seissmall = np.ndarray(shape=(len(ipick),Nsub,Ntimes))

    for ii in range(Ntimes):
        for jj in range(Nsub):
            seissmall[:,jj,ii] = seis[jj*Ismall+ipick,ii]
    fseis = np.fft.fft(seissmall,n=2**fftpower,axis=0)
    fseis = fseis[0:len(frq),:,:]
    
    ##doing it this way eliminates NaNs appearing when we have low
    ##amplitude
    fseis = np.cos(np.angle(fseis[I,:,:])) + 1j*np.sin(np.angle(fseis[I,:,:]))
    if np.isnan(fseis).any():
        print "NaN found"
        sys.exit(1)

    #indices of first and last frequency used
    Imin=I[0].min()
    Imax=I[0].max()

    outfile = os.path.join(matpath,'%s_%s'%(info['staname'],JD))
    if not os.path.isdir(matpath):
        os.mkdir(matpath)
    #scipy.io.savemat(outfile+'_info',info)
    fseis = np.squeeze(fseis[0,:,:,:])
    scipy.io.savemat(outfile,{'fseis':fseis,'Imin':Imin,'Imax':Imax,'frq':frq,'fftpower':fftpower})

outfile = 'BeamformInputData_start'
scipy.io.savemat(outfile,{'I':I,'JD':JD,'freq':freq,'Nsub':Nsub,'Ntimes':Ntimes,'matpath':matpath,'year':year})
    
