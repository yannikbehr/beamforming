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

sys.path.append(os.path.join(os.environ['PROC_SRC'],'misc'))
sys.path.append(os.path.join(os.environ['PROC_SRC'],'NA'))
from herman_surface_wave import herman_syn
from dinver_run import get_disp
import progressbar as pg
rcParams = {'backend':'Agg'}

DEBUG = False
def prep_beam_h(files,matfile,nhours=1,fmax=10.,fact=10,new=True,onebit=False):
    if new:
        ntimes = int(round(24/nhours))
        step = int(nhours*3600*fmax/fact)
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
            if onebit:
                seisbandn[i,j,:] = sign(seisbandn[i,j,:])
                seisbande[i,j,:] = sign(seisbande[i,j,:])


        if 1:
            fftpower = 7
            ismall = 2**fftpower
            ipick = arange(ismall)
            n=nhours*3600*df
            taper = cosTaper(len(ipick))
            nsub = int(np.floor(n/ismall)) # Number of time pieces -20 mins long each
            #seissmall = zeros((len(ipick),ntimes,nsub,nfiles))
            seissmalln = zeros((nfiles,ntimes,nsub,len(ipick)))
            seissmalle = zeros((nfiles,ntimes,nsub,len(ipick)))
            for ii in xrange(nfiles):
                for jj in xrange(ntimes):
                    for kk in xrange(nsub):
                        #seissmall[:,jj,kk,ii] = seisband[kk*ismall+ipick,jj,ii]
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

        sio.savemat(matfile,{'seissmalln':seissmalln,'seissmalle':seissmalle,'slats':slats,'slons':slons,'dt':dt})
        return seissmalln, seissmalle, meanlat, meanlon, slats, slons, dt
    else:
        a = sio.loadmat(matfile,struct_as_record=False)
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
            trN = SacIO(_fN,headonly=True)
            stN = trN.kstnm.rstrip()
            dtN = trN.delta
            trE = SacIO(_fE,headonly=True)
            stE = trE.kstnm.rstrip()
            dtE = trE.delta
            if dtE != dtN:
                print "sampling intervals are not identical for %s and %s"%(_fN,_fE)
                return
            if stN == stE:
                newlist.append((_fN,_fE))
                break
            
    return newlist,1./dtN

def syntrace(dist,wtype='rayleigh',wmodel='Herrmann'):
    if wmodel == 'Aki':
        model = ['4\n','3.  6.0 3.5 2.5 100. 200.\n',
                 '2.  3.4 2.0 2.3 100. 200.\n',
                 '5.  6.5 3.8 2.5 100. 200.\n',
                 '0.  8.0 4.7 3.3 500. 900.\n']
        model = ['4\n','3.  6.0 3.5 2.5 100. 200.\n',
                 '2.  6.2 3.6 2.5 100. 200.\n',
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
        t = linspace(0,256,256)
        repc = scint.splrep(2*pi*frc,c)
        repu = scint.splrep(2*pi*fru,u)
        fsum = 0

        for _w in om0:
            y = dom/2.*(t-(x/scint.splev(_w,repu)))
            f0 = dom/pi*sin(y)/y*cos(_w*t-_w*x/scint.splev(_w,repc))
            fsum += f0
        return t,fsum

    elif wmodel == 'Herrmann':
        model = ['3. 6.0 3.5 2.5 100.0 200.0 0.0 0.0 1.0 1.0\n',
                 '2. 3.4 2.0 2.3 100. 200. 0.0 0.0 1.0 1.0\n',
                 '5. 6.5 3.8 2.5 100. 200. 0.0 0.0 1.0 1.0\n',
                 '0. 8.0 4.7 3.3 500.0 900.0 0.0 0.0 1.0 1.0\n']
        model = ['3. 6.0 3.5 2.5 100.0 200.0 0.0 0.0 1.0 1.0\n',
                 '2. 6.2 3.6 2.5 100. 200. 0.0 0.0 1.0 1.0\n',
                 '5. 6.5 3.8 2.5 100. 200. 0.0 0.0 1.0 1.0\n',
                 '0. 8.0 4.7 3.3 500.0 900.0 0.0 0.0 1.0 1.0\n']
        z,r,t = herman_syn(model,dist,npts=256)
        time = linspace(0,256,256)
        if wtype == 'rayleigh':
            return time,z
        if wtype == 'love':
            return time,t
        else:
            print "incorrect wave type [rayleigh or love]"
    else:
        print "wmodel has to be either Herrmann or Aki"

def syntest(theta,zetax):
    dtheta = int(unique(diff(theta[:,0])))
    ind = int(round(225/dtheta))
    nsources, nstations = zetax.shape
    npts = 256
    rtraces = zeros((nstations,1,1,npts))
    ltraces = zeros((nstations,1,1,npts))
    taper = cosTaper(npts)
    for i,ddiff in enumerate(zetax[ind,:]/1000.):
        t,fsum = syntrace(500.+ddiff,wtype='rayleigh',wmodel='Herrmann')
        rtraces[i,0,0,:] = fsum*taper
        t,fsum = syntrace(500.+ddiff,wtype='love',wmodel='Herrmann')
        ltraces[i,0,0,:] = fsum*taper
    if 1:
        figure()
        subplot(2,1,1)
        plot(t,rtraces[0,0,0,:],label='Rayleigh')
        rtr = SacIO()
        rtr.fromarray(rtraces[0,0,0,:],distkm=500.+zetax[ind,0]/1000.)
        rtr.WriteSacBinary('rayleigh_synthetic.sac')
        legend(loc='upper left')
        subplot(2,1,2)
        plot(t,ltraces[0,0,0,:],label='Love')
        ltr = SacIO()
        ltr.fromarray(ltraces[0,0,0,:],distkm=500.+zetax[ind,0]/1000.)
        ltr.WriteSacBinary('love_synthetic.sac')
        legend(loc='upper left')
        xlabel('Time [s]')
    phi = theta[ind]*pi/180.-pi
    n = rtraces*sin(phi) + ltraces*cos(phi)
    e = rtraces*cos(phi) - ltraces*sin(phi)
    return n,e

def beamforming(seisn,seise,slowness,zetax,theta,dt,indices,new=True,matfile=None,freq_int=(0.1,0.4)):
    if new:
        if not DEBUG:
            widgets = ['horizontal beamforming: ', pg.Percentage(), ' ', pg.Bar('#'),
                       ' ', pg.ETA()]
            pbar = pg.ProgressBar(widgets=widgets, maxval=theta.size).start()
        nstat, ntimes, nsub, nfft = seisn.shape
        nsources = theta.size
        freq = fftfreq(nfft,dt)
        I = np.where((freq>freq_int[0]) & (freq<freq_int[1]))
        beamr = zeros((nsources,slowness.size,ntimes,nfft))
        beamt = zeros((nsources,slowness.size,ntimes,nfft))
        N = fft(seisn,n=nfft,axis=3)
        E = fft(seise,n=nfft,axis=3)
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
                            beamr[i,cc,tt,ww] += (1./(nstat*nstat)*abs(asarray(dot(conjugate(eT),dot(covr,e).T)))**2)/nsub
                            beamt[i,cc,tt,ww] += (1./(nstat*nstat)*abs(asarray(dot(conjugate(eT),dot(covt,e).T)))**2)/nsub
                            #beam[i,cc,tt] = abs(dot(Y,conjugate(e)))**2
        if not DEBUG:
            pbar.finish()
        sio.savemat(matfile,{'beamr':beamr,'beamt':beamt})
    else:
        a = sio.loadmat(matfile,struct_as_record=False)
        beamr = a['beamr']
        beamt = a['beamt']
    return beamr, beamt


def polar_plot(beam,theta,slowness,dt,nfft,wtype,fout=None):
    df = dt/nfft
    periods = [6.]
    idx = [int(1./(p*df)) for p in periods]
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

def polar_plot_test(beam,theta,slowness,indices,wtype,dt,nfft,resp=False): 
    df = dt/nfft
    theta = theta[:,0]
    slowness = slowness[0,:]
    p = []
    cmat = zeros((nfft,slowness.size))
    for ind in indices:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre/tre.max()
        p.append(1./(ind*df))
        cmat[ind,:] = tre[45,:]
                                
    figure()
    contourf(p,1./slowness,cmat[1::].T,100)

    if 1:
        xmin, xmax = xlim()
        model = ['4\n','3.  6.0 3.5 2.5 100. 200.\n',
                 '2.  3.4 2.0 2.3 100. 200.\n',
                 '5.  6.5 3.8 2.5 100. 200.\n',
                 '0.  8.0 4.7 3.3 500. 900.\n']
        model = ['4\n','3.  6.0 3.5 2.5 100. 200.\n',
                 '2.  6.2 3.6 2.5 100. 200.\n',
                 '5.  6.5 3.8 2.5 100. 200.\n',
                 '0.  8.0 4.7 3.3 500. 900.\n']
        rayc,lovc = get_disp(model,0.02,1.0,nmode=1)
        rayu,lovu = get_disp(model,0.02,1.0,gv=True)
        indlc = where(lovc[:,0] > 0.)
        indrc = where(rayc[:,0] > 0.)
        indlu = where(lovu[:,0] > 0.)
        indru = where(rayu[:,0] > 0.)
        if wtype == 'rayleigh':
            plot(1./rayc[indrc,0][0],1./rayc[indrc,1][0],label='Theoretical Rayleigh phase',color='black')
        elif wtype == 'love':
            plot(1./lovc[indlc,0][0],1./lovc[indlc,1][0],label='Theoretical Love phase',color='black')
        else:
            print 'wtype either has to be rayleigh or love'
        xlabel('Period [s]')
        ylabel('Velocity [km/s]')
        legend(loc='upper right')
        xlim(xmin,xmax)
    colorbar()
    ax = gca()
    ax.autoscale_view(tight=True)
    xlim(1,20)
    periods = [6.,20.]
    idx = [int(1./(p*df)) for p in periods]
    for ind in idx:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre/tre.max()
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
        ax.set_title("%s %ds period"%(wtype,1./(ind*df)))
        ax.grid(True)
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
    slowness = arange(0.03,0.505,0.005)
    slowness = slowness.reshape((1,slowness.size))
    indices = arange(nfft)
    newbeam = not os.path.isfile(matfile2)
    if nbeam:
        newbeam = True
    beam = arr_resp(nfft,dt,indices,slowness,zetax,theta,sta_origin_x,sta_origin_y,
                    new=newbeam,matfile=matfile2,src=False)
    polar_plot_resp(beam,theta,slowness,dt,nfft)
    show()
    
def test(datdir,nprep=False,nbeam=False,doplot=True):
    """
    Run synthetic test.
    """
    filesN = glob.glob(os.path.join(datdir,'ft_*.*HN.SAC'))
    filesE = glob.glob(os.path.join(datdir,'ft_*.*HE.SAC'))
    if len(filesN) < 2 or len(filesE) < 2:
        print "not enough files in ",datdir
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
    zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
    seisn, seise = syntest(theta,zetax)
    matfile2 = "%s_%s.mat"%('test_beam_h',temp)
    nstat, ntimes, nsub, nfft = seisn.shape
    indices = arange(nfft)
    newbeam = not os.path.isfile(matfile2)
    if nbeam:
        newbeam = True
    beamr,beamt = beamforming(seisn,seise,slowness,zetax,theta,dt,indices,
                              new=newbeam,freq_int=(0.1,0.4),matfile=matfile2)
    if doplot:
        polar_plot_test(beamr,theta,slowness,indices[1::],'rayleigh',dt,nfft)
        polar_plot_test(beamt,theta,slowness,indices[1::],'love',dt,nfft)
        show()

def main(datdir,nprep=False,nbeam=False,doplot=True,save=False,nostat=20):
    filesN = glob.glob(os.path.join(datdir,'ft_*.*HN.SAC'))
    filesE = glob.glob(os.path.join(datdir,'ft_*.*HE.SAC'))
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
    seisn, seise, meanlat, meanlon, slats, slons, dt = prep_beam_h(nlist,matfile1,
                                                                   nhours=1,fmax=sample_f,
                                                                   fact=sample_f,new=newprep)
    zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
    matfile2 = "%s_%s.mat"%('beam_h',temp)
    newbeam = not os.path.isfile(matfile2)
    
    if nbeam:
        newbeam = True
    nsources,ntimes,nsub,nfft = seisn.shape
    df = dt/nfft
    periods = [6.]
    #periods = [4.,5.,6.,7.,8.,9.,10.,12.,15.,18.]
    indices = [int(1./(p*df)) for p in periods]
    beamr,beamt = beamforming(seisn,seise,slowness,zetax,theta,dt,indices,
                              new=newbeam,freq_int=(0.1,0.4),matfile=matfile2)
    if doplot:
        fout = None
        if save:
            fout = matfile2.replace('.mat','_radial.png')
        polar_plot(beamr,theta,slowness,dt,nfft,'rayleigh',fout=fout)
        if save:
            fout = matfile2.replace('.mat','_transverse.png')
        polar_plot(beamt,theta,slowness,dt,nfft,'love',fout=fout)
        if not save:
            show()

if __name__ == '__main__':
    from optparse import OptionParser
    #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
    #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
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
