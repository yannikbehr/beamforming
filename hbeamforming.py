#!/usr/bin/env mypython
"""
Horizontal beamforming
"""

import os
import sys
import glob
import obspy.sac
from obspy.sac import *
from obspy.core import read
from pylab import *
import obspy.signal
import scipy.io as sio
from matplotlib import cm, rcParams
import ctypes as C
import pickle
import scipy.interpolate as scint
sys.path.append(os.path.join(os.environ['PROC_SRC'],'misc'))
sys.path.append(os.path.join(os.environ['PROC_SRC'],'NA'))
from dinver_run import get_disp
import progressbar as pg
rcParams = {'backend':'Agg'}

DEBUG = False
def prep_beam_h(files,matfile,nhours=1,fmax=10.,fact=10,new=True):
    if new:
        ntimes = int(round(24/nhours))
        step = nhours*3600*fmax/fact

        slons = array([])
        slats = array([])
        nfiles = len(files)
        seisbandn = zeros((nfiles,ntimes,step))
        seisbande = zeros((nfiles,ntimes,step))
        sigmas = []
        for i,_ftup in enumerate(files):
            trn = read(_ftup[0])[0]
            tre = read(_ftup[1])[0]
            slons = append(slons,trn.stats.sac.stlo)
            slats = append(slats,trn.stats.sac.stla)
            trn.downsample(decimation_factor=fact, strict_length=True)
            tre.downsample(decimation_factor=fact, strict_length=True)
            nptsn = trn.stats.npts
            nptse = tre.stats.npts
            df = trn.stats.sampling_rate
            dt = trn.stats.delta
            seis0n = zeros(24*3600*int(df))
            seis0e = zeros(24*3600*int(df))
            seis0n[0:nptsn] = trn.data
            seis0e[0:nptse] = tre.data
            seis0n -= seis0n.mean()
            seis0e -= seis0e.mean()
            for j in xrange(ntimes):
                ilow = j*step
                iup = (j+1)*step
                seisbandn[i,j,:] = seis0n[ilow:iup]
                seisbandn[i,j,:] -= seisbandn[i,j,:].mean()
                seisbande[i,j,:] = seis0e[ilow:iup]
                seisbande[i,j,:] -= seisbande[i,j,:].mean()

        if 1:
            fftpower = 7
            ismall = 2**fftpower
            ipick = arange(ismall)
            n=nhours*3600*df
            nsub = int(np.floor(n/ismall)) # Number of time pieces -20 mins long each
            #seissmall = zeros((len(ipick),ntimes,nsub,nfiles))
            seissmalln = zeros((nfiles,ntimes,nsub,len(ipick)))
            seissmalle = zeros((nfiles,ntimes,nsub,len(ipick)))
            for ii in xrange(nfiles):
                for jj in xrange(ntimes):
                    for kk in xrange(nsub):
                        #seissmall[:,jj,kk,ii] = seisband[kk*ismall+ipick,jj,ii]
                        seissmalln[ii,jj,kk,:] = seisbandn[ii,jj,kk*ismall+ipick]
                        seissmalle[ii,jj,kk,:] = seisbande[ii,jj,kk*ismall+ipick]

        LonLref= 165
        LonUref= 179.9
        LatLref= -48
        LatUref= -34
        stacoord=vstack((slons,slats))
        ##Find the stations which belong to this grid
        idx = where((stacoord[0] >= LonLref) & (stacoord[0] <= LonUref) & \
                    (stacoord[1] >= LatLref) & (stacoord[1] <= LatUref))
        #nb. this simple 'mean' calc only works if we don't cross lat=0 or lon=180
        meanlat = slats.mean()
        meanlon = slons.mean()

        sio.savemat(matfile,{'seissmalln':seissmalln,'seissmalle':seissmalle,'slats':slats,'slons':slons,'dt':dt})
        return seissmalln, seissmalle, meanlat, meanlon, slats, slons, dt
    else:
        a = sio.loadmat(matfile)
        seissmalln = a['seissmalln']
        seissmalle = a['seissmalle']
        slats = a['slats']
        slons = a['slons']
        dt = a['dt'][0][0]
        slats = slats.reshape(slats.shape[0],)
        slons = slons.reshape(slons.shape[0],)
        meanlat = slats.mean()
        meanlon = slons.mean()
        return seissmalln, seissmalle, meanlat, meanlon, slats, slons, dt


def calc_steer(slats,slons):
    #theta= arange(0,362,2)
    theta= arange(0,365,5)
    theta = theta.reshape((theta.size,1))
    sta_origin_dist = array([])
    sta_origin_bearing = array([])
    for lat,lon in zip(slats,slons):
        dist, az, baz = obspy.signal.rotate.gps2DistAzimuth(meanlat,meanlon,lat,lon)
        sta_origin_dist = append(sta_origin_dist,dist)
        sta_origin_bearing = append(sta_origin_bearing,az)
    sta_origin_x = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
    sta_origin_y = sta_origin_dist*sin(sta_origin_bearing*pi/180.)
    #zeta_x = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
    #zeta_y = sta_origin_dist*sin(sta_origin_bearing*pi/180.)
    zeta_x = -cos(theta*pi/180.)
    zeta_y = -sin(theta*pi/180.)
    #dot product betwen zeta and x
    zetax = zeta_x*sta_origin_x+zeta_y*sta_origin_y
    #slowness in s/km
    #slowness = arange(0.03,0.505,0.005)
    slowness = arange(0.125,0.51,0.01)
    slowness = slowness.reshape((1,slowness.size))
    return zetax,theta,slowness,sta_origin_x,sta_origin_y


def mkfilelist(filesN, filesE):
    """
    Order the file lists for the north and east component.
    """
    newlist = []
    for _fN in filesN:
        for _fE in filesE:
            a = os.path.basename(_fN).split('_')
            stN = a[2].split('.')[0]
            a = os.path.basename(_fE).split('_')
            stE = a[2].split('.')[0]
            if stN == stE:
                newlist.append((_fN,_fE))
                break
            
    return newlist

def syntrace(dist,wtype='rayleigh'):
    model = ['4\n','3.  6.0 3.5 2.5 100. 200.\n',
             '2.  3.4 2.0 2.3 100. 200.\n',
             '5.  6.5 3.8 2.5 100. 200.\n',
             '0.  8.0 4.7 3.3 500. 900.\n']
    rayc,lovc = get_disp(model,0.02,1.0,nmode=1)
    rayu,lovu = get_disp(model,0.02,1.0,gv=True)
    indlc = where(lovc[:,0] > 0.)
    indrc = where(rayc[:,0] > 0.)
    indlu = where(lovu[:,0] > 0.)
    indru = where(rayu[:,0] > 0.)
    if 0:
        plot(1./lovc[indlc,0][0],1./lovc[indlc,1][0],label='Love phase')
        plot(1./rayc[indrc,0][0],1./rayc[indrc,1][0],label='Rayleigh phase')
        plot(1./lovu[indlu,0][0],1./lovu[indlu,1][0],label='Love group')
        plot(1./rayu[indru,0][0],1./rayu[indru,1][0],label='Rayleigh group')
        xlabel('Period [s]')
        ylabel('Velocity [km/s]')
        legend(loc='lower right')


    ################# calculate a synthetic seismogram by summing the
    ################# contributions of wavepackages centered around
    ################# frequency intervals of 0.01 Hz from 0 to 1 Hz
    df = 0.01
    fint = arange(0.02,1.01,0.01)
    if wtype == 'rayleigh':
        frc = rayc[indrc,0][0]
        fru = rayu[indru,0][0]
        c = 1./rayc[indrc,1][0]
        u = 1./rayu[indru,1][0]
    elif wtype == 'love':
        frc = lovc[indlc,0][0]
        fru = lovu[indlu,0][0]
        c = 1./lovc[indlc,1][0]
        u = 1./lovu[indlu,1][0]
    else:
        print "incorrect wave type [rayleigh or love]"
    dom = 2*pi*df
    om0 = 0.005+fint[:-1]*2*pi
    om0 = fint*2*pi
    x = dist
    t = linspace(80,208,128)
    repc = scint.splrep(2*pi*frc,c)
    repu = scint.splrep(2*pi*fru,u)
    fsum = 0

    for _w in om0:
        y = dom/2.*(t-(x/scint.splev(_w,repu)))
        f0 = dom/pi*sin(y)/y*cos(_w*t-_w*x/scint.splev(_w,repc))
        fsum += f0

    return t,fsum

def syntest(theta,zetax):
    dtheta = int(unique(diff(theta[:,0])))
    ind = int(round(225/dtheta))
    nsources, nstations = zetax.shape
    rtraces = zeros((nstations,1,1,128))
    ltraces = zeros((nstations,1,1,128))
    for i,ddiff in enumerate(zetax[ind,:]/1000.):
        t,fsum = syntrace(500.+ddiff,wtype='rayleigh')
        rtraces[i,0,0,:] = fsum
        t,fsum = syntrace(500.+ddiff,wtype='love')
        ltraces[i,0,0,:] = fsum
    if 1:
        figure()
        subplot(2,1,1)
        plot(t,rtraces[0,0,0,:],label='Rayleigh')
        legend(loc='upper left')
        subplot(2,1,2)
        plot(t,ltraces[0,0,0,:],label='Love')
        legend(loc='upper left')
        xlabel('Time [s]')
    phi = theta[ind]*pi/180.-pi
    n = rtraces*sin(phi) + ltraces*cos(phi)
    e = rtraces*cos(phi) - ltraces*sin(phi)
    return n,e

def beamforming(seisn,seise,slowness,zetax,theta,dt,new=True,matfile=None,freq_int=(0.1,0.4)):
    if new:
        if not DEBUG:
            widgets = ['horizontal beamforming: ', pg.Percentage(), ' ', pg.Bar('#'),
                       ' ', pg.ETA()]
            pbar = pg.ProgressBar(widgets=widgets, maxval=theta.size).start()
        _p = 6.
        _f = 1./_p
        periods = arange(4.,11.)
        nstat, ntimes, nsub, nfft = seisn.shape
        nsources = theta.size
        freq = fftfreq(nfft,dt)
        I = np.where((freq>freq_int[0]) & (freq<freq_int[1]))
        beamr = zeros((nsources,slowness.size,ntimes,nfft))
        beamt = zeros((nsources,slowness.size,ntimes,nfft))
        df = dt/nfft
        idx = [int(1./(p*df)) for p in periods]
        ind = int(_f/df)
        N = fft(seisn,n=nfft,axis=3)
        E = fft(seise,n=nfft,axis=3)
        for i,az in enumerate(theta):
            if not DEBUG:
                pbar.update(i)
            daz = az*pi/180.-pi
            R = cos(daz)*N + sin(daz)*E
            T = -sin(daz)*N + cos(daz)*E
            dist = zetax[i,:]
            #for ww in [ind]:
            for ww in idx:
                FF = freq[ww]
                omega = 2*pi*FF
                for tt in xrange(ntimes):
                #for tt in [0]:
                    for TT in xrange(nsub):
                    #for TT in [0]:
                        Yt = asmatrix(squeeze(T[:,tt,TT,ww]))
                        Yr = asmatrix(squeeze(R[:,tt,TT,ww]))
                        #Y = squeeze(R[:,tt,TT,ww])
                        YtT = Yt.T.copy()
                        YrT = Yr.T.copy()
                        covr = dot(YrT,conjugate(Yr))
                        covt = dot(YtT,conjugate(Yt))
                        for cc in xrange(slowness.shape[1]):
                            velocity = 1./slowness[0][cc]*1000
                            e = exp(-1j*dist*omega/velocity)
                            eT = e.T.copy()
                            beamr[i,cc,tt,ww] += (abs(asarray(dot(conjugate(eT),dot(covr,e).T)))**2)/nsub
                            beamt[i,cc,tt,ww] += (abs(asarray(dot(conjugate(eT),dot(covt,e).T)))**2)/nsub
                            #beam[i,cc,tt] = abs(dot(Y,conjugate(e)))**2
        if not DEBUG:
            pbar.finish()
        sio.savemat(matfile,{'beamr':beamr,'beamt':beamt})
    else:
        a = sio.loadmat(matfile)
        beamr = a['beamr']
        beamt = a['beamt']
    return beamr, beamt


def polar_plot(beam,theta,slowness):
    theta = theta[:,0]
    slowness = slowness[0,:]
    tre = squeeze(beam[:,:,:,21])
    tre = tre.mean(axis=2)
    tre = tre-tre.max()
    fig = figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1,projection='polar')
    cmap = cm.get_cmap('jet')
    ax.contourf((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                100,cmap=cmap,antialiased=True,
                linewidths=0.1,linstyles='dotted')
    ax.contour((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                100,cmap=cmap)
    ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                      labels=['90','45','0','315','270','225','180','135'])
    ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
    ax.grid(True)
    #ax.set_title(os.path.basename(_d))
    ax.set_rmax(0.5)

def polar_plot_test(beam,theta,slowness):
    theta = theta[:,0]
    slowness = slowness[0,:]
    tre = squeeze(beam[:,:,:,21])
    #tre = tre.mean(axis=2)
    tre = tre-tre.max()
    fig = figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1,projection='polar')
    cmap = cm.get_cmap('jet')
    ax.contourf((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                100,cmap=cmap,antialiased=True,
                linewidths=0.1,linstyles='dotted')
    ax.contour((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                100,cmap=cmap)
    ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                      labels=['90','45','0','315','270','225','180','135'])
    ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
    ax.grid(True)
    #ax.set_title(os.path.basename(_d))
    ax.set_rmax(0.5)


if __name__ == '__main__':
    if 0:
        datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
        #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
        filesN = glob.glob(os.path.join(datdir,'ft_grid*.HHN.SAC'))
        filesE = glob.glob(os.path.join(datdir,'ft_grid*.HHE.SAC'))
        matfile = 'prep_beam_h_2001_2_22.mat'
        nlist = mkfilelist(filesN, filesE)
        if DEBUG:
            print "preparing data"
        seisn, seise, meanlat, meanlon, slats, slons, dt = prep_beam_h(nlist,matfile,
                                                                       nhours=1,fmax=10.,
                                                                       fact=10,new=False)
        zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
        matfile = 'beam_h_2001_2_22.mat'
        if DEBUG:
            print "beamforming"
        beamr,beamt = beamforming(seisn,seise,slowness,zetax,theta,dt,
                                  new=True,freq_int=(0.1,0.4),matfile=matfile)
        if 1:
            if DEBUG:
                print "plotting"
            polar_plot(beamr,theta,slowness)
            polar_plot(beamt,theta,slowness)
            show()
    if 1:
        datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
        #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
        filesN = glob.glob(os.path.join(datdir,'ft_grid*.HHN.SAC'))
        filesE = glob.glob(os.path.join(datdir,'ft_grid*.HHE.SAC'))
        matfile = 'prep_beam_h_2001_2_22.mat'
        nlist = mkfilelist(filesN, filesE)
        seisn, seise, meanlat, meanlon, slats, slons, dt = prep_beam_h(nlist,matfile,
                                                                       nhours=1,fmax=10.,
                                                                       fact=10,new=True)
        zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
        seisn, seise = syntest(theta,zetax)
        beamr,beamt = beamforming(seisn,seise,slowness,zetax,theta,dt,
                                  new=True,freq_int=(0.1,0.4),matfile=matfile)
        polar_plot_test(beamr,theta,slowness)
        polar_plot_test(beamt,theta,slowness)
        show()
