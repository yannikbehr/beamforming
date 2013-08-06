#!/usr/bin/env mypython
"""
Horizontal beamforming
"""

import os
import sys
import glob
from obspy.sac import *
from obspy.core import read
if os.environ.has_key('SGE_TASK_ID'):
    import matplotlib
    matplotlib.use('Agg')
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from pylab import *
import obspy.signal
from obspy.signal.invsim import cosTaper
import scipy.io as sio
from matplotlib import cm, rcParams
import scipy.interpolate as scint
import progressbar as pg
rcParams = {'backend':'Agg'}

DEBUG = False
def prep_beam_h(files,matfile,nhours=1,fmax=10.,fact=10,new=True,onebit=False,fftpower=7,laura=False):
    if new:
        ntimes = int(round(24/nhours))
        step = int(nhours*3600*fmax/fact)
        slons = array([])
        slats = array([])
        nfiles = len(files)
        seisbandn = zeros((nfiles,ntimes,step))
        seisbande = zeros((nfiles,ntimes,step))
        sigmas = []
        if laura:
            freqs = (fmax/fact)/2.*np.arange(1,2**(fftpower-1)+1)/2**(fftpower-1)
        else:
            freqs = fftfreq(2**fftpower,1./(fmax/fact))
        for i,_ftup in enumerate(files):
            trn = read(_ftup[0])[0]
            tre = read(_ftup[1])[0]
            slons = append(slons,trn.stats.sac.stlo)
            slats = append(slats,trn.stats.sac.stla)
            trn.data -= trn.data.mean()
            trn.filter("bandpass",freqmin=0.02,freqmax=0.4,corners=4,zerophase=True)
            tre.data -= tre.data.mean()
            tre.filter("bandpass",freqmin=0.02,freqmax=0.4,corners=4,zerophase=True)
            if (trn.stats.sampling_rate - 1.0) > 0.0001 and (tre.stats.sampling_rate - 1.0) > 0.0001:
                trn.downsample(decimation_factor=fact, strict_length=True,no_filter=True)
                tre.downsample(decimation_factor=fact, strict_length=True,no_filter=True)
            else:
                continue
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
            if onebit:
                seisbandn[i,j,:] = sign(seisbandn[i,j,:])
                seisbande[i,j,:] = sign(seisbande[i,j,:])


        ismall = 2**fftpower
        ipick = arange(ismall)
        n=nhours*3600*df
        taper = cosTaper(len(ipick))
        nsub = int(np.floor(n/ismall)) # Number of time pieces -20 mins long each
        seissmalln = zeros((nfiles,ntimes,nsub,len(ipick)))
        seissmalle = zeros((nfiles,ntimes,nsub,len(ipick)))
        for ii in xrange(nfiles):
            for jj in xrange(ntimes):
                for kk in xrange(nsub):
                    seissmalln[ii,jj,kk,:] = seisbandn[ii,jj,kk*ismall+ipick]*taper
                    seissmalle[ii,jj,kk,:] = seisbande[ii,jj,kk*ismall+ipick]*taper

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

        sio.savemat(matfile,{'seissmalln':seissmalln,'seissmalle':seissmalle,'slats':slats,'slons':slons,'dt':dt,'files': files,'freqs':freqs})
        return seissmalln, seissmalle, freqs, meanlat, meanlon, slats, slons, dt
    else:
        a = sio.loadmat(matfile,struct_as_record=False)
        seissmalln = a['seissmalln']
        seissmalle = a['seissmalle']
        slats = a['slats']
        slons = a['slons']
        dt = a['dt'][0][0]
        freqs = a['freqs']
        slats = slats.reshape(slats.shape[0],)
        slons = slons.reshape(slons.shape[0],)
        meanlat = slats.mean()
        meanlon = slons.mean()
        return seissmalln, seissmalle, freqs, meanlat, meanlon, slats, slons, dt


def calc_steer(slats,slons):
    theta= arange(0,362,2)
    #theta= arange(0,365,5)
    theta = theta.reshape((theta.size,1))
    sta_origin_dist = array([])
    sta_origin_bearing = array([])
    meanlat = slats.mean()
    meanlon = slons.mean()
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
    slowness = arange(0.03,0.605,0.005)
    #slowness = arange(0.125,0.51,0.01)
    slowness = slowness.reshape((1,slowness.size))
    return zetax,theta,slowness,sta_origin_x,sta_origin_y


def mkfilelist(filesN, filesE):
    """
    Order the file lists for the north and east component.
    """
    newlist = []
    for _fN in filesN:
        for _fE in filesE:
            trN = SacIO(_fN,headonly=True)
            stN = trN.kstnm.rstrip()
            dtN = trN.delta
            trE = SacIO(_fE,headonly=True)
            stE = trE.kstnm.rstrip()
            dtE = trE.delta
            if stN == stE:
                if dtE != dtN:
                    print "sampling intervals are not identical for %s and %s"%(_fN,_fE)
                    break
                else:
                    newlist.append((_fN,_fE))
                    break
            
    return newlist,1./dtN

def beamforming(seisn,seise,freqs,slowness,zetax,theta,dt,periods,new=True,matfile=None,freq_int=(0.02,0.4)):
    if new:
        if not DEBUG:
            widgets = ['horizontal beamforming: ', pg.Percentage(), ' ', pg.Bar('#'),
                       ' ', pg.ETA()]
            pbar = pg.ProgressBar(widgets=widgets, maxval=theta.size).start()
        nstat, ntimes, nsub, nfft = seisn.shape
        nsources = theta.size
        beamr = zeros((nsources,slowness.size,ntimes,nfft))
        beamt = zeros((nsources,slowness.size,ntimes,nfft))
        N = fft(seisn,n=nfft,axis=3)
        E = fft(seise,n=nfft,axis=3)
        ind = np.where((freqs>freq_int[0])&(freqs<freq_int[1]))[0]
        N = N[:,:,:,ind]
        E = E[:,:,:,ind]
        freqs = freqs[ind]
        indices = [argmin(abs(freqs - 1./p)) for p in periods]

        for i,az in enumerate(theta):
            if not DEBUG:
                pbar.update(i)
            daz = az*pi/180.-pi
            R = cos(daz)*N + sin(daz)*E
            T = -sin(daz)*N + cos(daz)*E
            ### ignore amplitude and keep only phase
            R = exp(angle(R)*1j)
            T = exp(angle(T)*1j)
            dist = zetax[i,:]
            #for ww in [ind]:
            for ww in indices:
                FF = freqs[ww]
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
                            beamr[i,cc,tt,ww] += (1./(nstat*nstat)*abs(asarray(dot(conjugate(eT),dot(covr,e).T)))**2)/nsub
                            beamt[i,cc,tt,ww] += (1./(nstat*nstat)*abs(asarray(dot(conjugate(eT),dot(covt,e).T)))**2)/nsub
                            #beam[i,cc,tt] = abs(dot(Y,conjugate(e)))**2
        if not DEBUG:
            pbar.finish()
        sio.savemat(matfile,{'beamr':beamr,'beamt':beamt,'theta':theta,'slowness':slowness,'freqs':freqs})
    else:
        a = sio.loadmat(matfile,struct_as_record=False)
        beamr = a['beamr']
        beamt = a['beamt']
        freqs = a['freqs']
    return beamr, beamt, freqs


def polar_plot(beam,theta,freqs,slowness,dt,nfft,wtype,fout=None):
    df = dt/nfft
    periods = [6.]
    idx = [int(1./(p*df)) for p in periods]
    idx = [argmin(abs(freqs - 1./p)) for p in periods]
    theta = theta[:,0]
    slowness = slowness[0,:]
    for ind in idx:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre.mean(axis=2)
        #tre = tre-tre.max()
        fig = figure(figsize=(6,6))
        #ax = fig.add_subplot(1,1,1,projection='polar')
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax  = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
        cmap = cm.get_cmap('jet')
        ax.contourf((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                    100,cmap=cmap,antialiased=True,
                    linstyles='dotted')
        ax.contour((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                   100,cmap=cmap)
        ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                          labels=['90','45','0','315','270','225','180','135'])
        ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
        ax.grid(True)
        ax.set_title("%s %ds period"%(wtype,1./(ind*df)))
        ax.set_rmax(0.5)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=tre.min(), vmax=tre.max()))
        if fout is not None:
            savefig(fout)

def polar_plot_resp(beam,theta,slowness,dt,nfft):
    df = dt/nfft
    periods = [6.]
    idx = [int(1./(p*df)) for p in periods]
    theta = theta[:,0]
    slowness = slowness[0,:]
    for ind in idx:
        tre = squeeze(beam[:,:,ind])
        tre = tre-tre.max()
        fig = figure(figsize=(6,6))
        ax = fig.add_subplot(1,1,1,projection='polar')
        cmap = cm.get_cmap('jet')
        ax.contourf((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                    100,cmap=cmap,antialiased=True,
                    linstyles='dotted')
        ax.contour((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                   100,cmap=cmap)
        ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                          labels=['90','45','0','315','270','225','180','135'])
        ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
        ax.grid(True)
        ax.set_title("%s %ds period"%('Array response',1./(ind*df)))
        ax.set_rmax(0.5)

def arr_resp(nfft,dt,indices,slowness,zetax,theta,sta_origin_x,sta_origin_y,
             new=True,matfile=None,src=False):
    """
    calculate array response
    """
    if new:
        theta1 = 90
        zeta_x = -cos(theta1*pi/180.)
        zeta_y = -sin(theta1*pi/180.)
        zeta_src = zeta_x*sta_origin_x + zeta_y*sta_origin_y
        c1 = 3000
        beam = zeros((zetax.shape[0],slowness.size,nfft))
        freqs = fftfreq(nfft,dt)
        for ww in indices:
            FF = freqs[ww]
            for cc in xrange(slowness.shape[1]):
                omega = 2*pi*FF
                velocity = 1./slowness[0][cc]*1000
                e_steer=exp(-1j*zetax*omega/velocity).T
                if src:
                    e_src = exp(-1j*zeta_src*omega/c1).T
                    Y = multiply(ones((zetax.shape[1],1),'complex128'),atleast_2d(e_src).T)
                    Y = ones((zetax.shape[1],1),'complex128')
                    R = dot(Y,conjugate(Y).T)
                else:
                    R = ones((zetax.shape[1],zetax.shape[1]))
                beam[:,cc,ww] = atleast_2d(sum(abs(asarray(dot(conjugate(e_steer.T),dot(R,e_steer))))**2,axis=1))
        sio.savemat(matfile,{'beam':beam})
        return beam
    else:
        beam = sio.loadmat(matfile)['beam']
        return beam

def response(datdir,nprep=False,nbeam=False,doplot=True):
    """
    Calculate the theoretical array response
    """
    filesN = glob.glob(os.path.join(datdir,'ft_*.*HN.SAC'))
    filesE = glob.glob(os.path.join(datdir,'ft_*.*HE.SAC'))
    if datdir.endswith('/'):
        datdir = datdir[0:-1]
    temp = os.path.basename(datdir)
    matfile1 = "%s_%s.mat"%('prep_beam_h',temp)
    nlist,sample_f = mkfilelist(filesN, filesE)
    sample_f = int(round(sample_f))
    newprep = not os.path.isfile(matfile1)
    if nprep:
        newprep = True
    seisn, seise, meanlat, meanlon, slats, slons, dt = prep_beam_h(nlist,matfile1,
                                                                   nhours=1,fmax=sample_f,
                                                                   fact=sample_f,new=newprep)
    print meanlat, meanlon
    zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
    matfile2 = "%s_%s.mat"%('array_response',temp)
    nstat, ntimes, nsub, nfft = seisn.shape
    #slowness = arange(0.03,0.505,0.005)
    slowness = slowness.reshape((1,slowness.size))
    indices = arange(nfft)
    newbeam = not os.path.isfile(matfile2)
    if nbeam:
        newbeam = True
    beam = arr_resp(nfft,dt,indices,slowness,zetax,theta,sta_origin_x,sta_origin_y,
                    new=newbeam,matfile=matfile2,src=False)
    polar_plot_resp(beam,theta,slowness,dt,nfft)
    show()
    
def main(datdir,nprep=False,nbeam=False,doplot=True,save=False,nostat=20):
    filesN = glob.glob(os.path.join(datdir,'[!^ft]*.*HN.SAC'))
    filesE = glob.glob(os.path.join(datdir,'[!^ft]*.*HE.SAC'))
    if len(filesN) < nostat or len(filesE) < nostat:
        print "not enough files in ",datdir
        return
    if datdir.endswith('/'):
        datdir = datdir[0:-1]
    temp = os.path.basename(datdir)
    matfile1 = "%s_%s.mat"%('prep_beam_h',temp)
    nlist,sample_f = mkfilelist(filesN, filesE)
    sample_f = int(round(sample_f))
    newprep = not os.path.isfile(matfile1)
    if nprep:
        newprep = True
    seisn, seise, freqs, meanlat, meanlon, slats, slons, dt = prep_beam_h(nlist,matfile1,
                                                                   nhours=1,fmax=sample_f,
                                                                   fact=sample_f,new=newprep,
                                                                   laura=False)
    zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
    matfile2 = "%s_%s.mat"%('beam_h',temp)
    newbeam = not os.path.isfile(matfile2)
    if nbeam:
        newbeam = True
    nsources,ntimes,nsub,nfft = seisn.shape
    #periods = [6.]
    periods = [4.,5.,6.,7.,8.,9.,10.,12.,15.,18.]
    beamr,beamt,freqs = beamforming(seisn,seise,freqs,slowness,zetax,theta,dt,periods,
                              new=newbeam,matfile=matfile2)
    if doplot:
        fout = None
        if save:
            fout = matfile2.replace('.mat','_radial.png')
        polar_plot(beamr,theta,freqs,slowness,dt,nfft,'rayleigh',fout=fout)
        if save:
            fout = matfile2.replace('.mat','_transverse.png')
        polar_plot(beamt,theta,freqs,slowness,dt,nfft,'love',fout=fout)
        if not save:
            show()

def proc_main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-t","--test",dest="test",action="store_true",
                      help="Run a synthetic test using the given network layout.",
                      default=False)
    parser.add_option("-b","--beam",dest="beam",action="store_true",
                      help="Recalculate beam.",
                      default=False)
    parser.add_option("-d","--data",dest="data",action="store_true",
                      help="Prepare data.",
                      default=False)
    parser.add_option("-s","--save",dest="save",action="store_true",
                      help="Save output to file.",
                      default=False)
    parser.add_option("--noplot",dest="plot",action="store_false",
                      help="Don't plot anything.",
                      default=True)
    parser.add_option("--resp",dest="aresp",action="store_true",
                      help="Calculate array response.",
                      default=False)
    parser.add_option("--nstat",dest="nstat",
                      help="Minimum number of stations.",
                      default=20)
    
    (opts,args) = parser.parse_args()
    if opts.test:
        datdir = args[0]
        test(datdir,nbeam=opts.beam,nprep=opts.data,doplot=opts.plot)
    elif opts.aresp:
        datdir = args[0]
        response(datdir,nbeam=opts.beam,nprep=opts.data,doplot=opts.plot)
    else:
        datdir = args[0]
        main(datdir,nbeam=opts.beam,nprep=opts.data,doplot=opts.plot,
             save=opts.save,nostat=int(opts.nstat))
    

if __name__ == '__main__':
    try:
        __IPYTHON__
    except:
        proc_main()
    else:
        print 'running in ipython'
        #datdir = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/sacfiles/5Hz/2002/Jul/2002_7_12_0_0_0'
        #main(datdir,nbeam=False,nprep=False,doplot=True,save=False,nostat=20)
