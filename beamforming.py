#!/usr/bin/env mypython

"""
next try to rewrite laura's beamformer
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
from obspy.core.util import gps2DistAzimuth
from obspy.signal.invsim import cosTaper
from obspy.signal.filter import bandpass
import scipy.io as sio
from matplotlib import cm, rcParams
import ctypes as C
import pickle
import scipy.interpolate as scint
import progressbar as pg
from smooth import smooth
rcParams = {'backend':'Agg'}

DEBUG = False
def prep_beam(files, matfile, nhours=1, fmax=10., threshold_std=0.5, onebit=False,
              tempfilter=False, specwhite=True, timenorm=False, fact=10, new=True,
              fftpower=7, freq_int=(0.02, 0.4), laura=False):
    if new:
        ntimes = int(round(24 / nhours))
        step = nhours * 3600 * fmax / fact
        stations = []
        slons = array([])
        slats = array([])
        nfiles = len(files)
        seisband = zeros((nfiles, ntimes, step))
        if laura:
            freqs = (fmax / fact) / 2.*np.arange(1, 2 ** (fftpower - 1) + 1) / 2 ** (fftpower - 1)
        else:
            freqs = fftfreq(2 ** fftpower, 1. / (fmax / fact))
        for i, _f in enumerate(files):
            tr = read(_f)[0]
            staname = tr.stats.station
            kcomp = tr.stats.channel
            slat = tr.stats.sac.stla
            slon = tr.stats.sac.stlo
            slons = append(slons, tr.stats.sac.stlo)
            slats = append(slats, tr.stats.sac.stla)
            if staname not in stations:
                stations.append(staname)
            tr.data -= tr.data.mean()
            tr.filter("bandpass", freqmin=0.02, freqmax=0.4, corners=4, zerophase=True)
            if (tr.stats.sampling_rate - 1.0) > 0.0001:
                tr.decimate(fact, strict_length=True, no_filter=True)
            else:
                continue
            npts = tr.stats.npts
            df = tr.stats.sampling_rate
            dt = tr.stats.delta
            seis0 = zeros(24 * 3600 * int(df))
            taper = cosTaper(tr.stats.npts)
            #detrend(tr.data)
            istart = int(round(((tr.stats.starttime.hour * 60 + tr.stats.starttime.minute) * 60\
                          + tr.stats.starttime.second) * df))
            if timenorm:
                smoothdata = smooth(abs(tr.data), window_len=257, window='flat')
                tr.data /= smoothdata
            #seis0[istart:(istart+npts)] = tr.data*taper
            try:
                seis0[istart:(istart + npts)] = tr.data
            except ValueError, e:
                print _f
                raise ValueError

            seis0 -= seis0.mean()
            for j in xrange(ntimes):
                ilow = j * step
                iup = (j + 1) * step
                seisband[i, j, :] = seis0[ilow:iup]
                seisband[i, j, :] -= seisband[i, j, :].mean()
                if onebit:
                    seisband[i, j, :] = sign(seisband[i, j, :])
            if tempfilter:
                sigmas = seisband[i, :, :].std(axis=1, ddof=1)
                sgm = ma.masked_equal(array(sigmas), 0.).compressed()
                sigma = sqrt(sum(sgm ** 2) / sgm.size)
                threshold = threshold_std * sigma
                seisband[i] = where(abs(seisband[i]) > threshold, threshold * sign(seisband[i]), seisband[i])
                seisband[i] = apply_along_axis(lambda e: e - e.mean(), 1, seisband[i])


        ismall = 2 ** fftpower
        ipick = arange(ismall)
        taper = cosTaper(len(ipick))
        n = nhours * 3600 * df
        nsub = int(np.floor(n / ismall)) # Number of time pieces -20 mins long each
        seissmall = zeros((nfiles, ntimes, nsub, len(ipick)))
        for ii in xrange(nfiles):
            for jj in xrange(ntimes):
                for kk in xrange(nsub):
                    seissmall[ii, jj, kk, :] = seisband[ii, jj, kk * ismall + ipick] * taper

        fseis = fft(seissmall, n=2 ** fftpower, axis=3)
        if np.isnan(fseis).any():
            print "NaN found"
            return
        ind = np.where((freqs > freq_int[0]) & (freqs < freq_int[1]))[0]
        fseis = fseis[:, :, :, ind]
        if specwhite:
            fseis = exp(angle(fseis) * 1j)

        LonLref = 165
        LonUref = 179.9
        LatLref = -48
        LatUref = -34
        stacoord = vstack((slons, slats))
        ##Find the stations which belong to this grid
        idx = where((stacoord[0] >= LonLref) & (stacoord[0] <= LonUref) & \
                    (stacoord[1] >= LatLref) & (stacoord[1] <= LatUref))
        #nb. this simple 'mean' calc only works if we don't cross lat=0 or lon=180
        meanlat = slats.mean()
        meanlon = slons.mean()

        sio.savemat(matfile, {'fseis':fseis, 'seissmall':seissmall, 'slats':slats, 'slons':slons, 'dt':dt, 'files':files, 'freqs':freqs[ind]})
        return fseis, freqs[ind], meanlat, meanlon, slats, slons, dt, seissmall
    else:
        a = sio.loadmat(matfile)
        fseis = a['fseis']
        seissmall = a['seissmall']
        slats = a['slats']
        slons = a['slons']
        dt = a['dt'][0][0]
        freqs = a['freqs']
        slats = slats.reshape(slats.shape[0],)
        slons = slons.reshape(slons.shape[0],)
        meanlat = slats.mean()
        meanlon = slons.mean()
        return fseis, freqs, meanlat, meanlon, slats, slons, dt, seissmall


def calc_steer(slats, slons):
    theta = arange(0, 362, 2)
    #theta= arange(0,365,5)
    theta = theta.reshape((theta.size, 1))
    sta_origin_dist = array([])
    sta_origin_bearing = array([])
    meanlat = slats.mean()
    meanlon = slons.mean()
    for lat, lon in zip(slats, slons):
        dist, az, baz = gps2DistAzimuth(meanlat, meanlon, lat, lon)
        sta_origin_dist = append(sta_origin_dist, dist)
        sta_origin_bearing = append(sta_origin_bearing, az)
    sta_origin_x = sta_origin_dist * cos(sta_origin_bearing * pi / 180.)
    sta_origin_y = sta_origin_dist * sin(sta_origin_bearing * pi / 180.)
    #zeta_x = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
    #zeta_y = sta_origin_dist*sin(sta_origin_bearing*pi/180.)
    zeta_x = -cos(theta * pi / 180.)
    zeta_y = -sin(theta * pi / 180.)
    #dot product betwen zeta and x
    zetax = zeta_x * sta_origin_x + zeta_y * sta_origin_y
    #slowness in s/km
    slowness = arange(0.03, 0.505, 0.005)
    #slowness = arange(0.125,0.51,0.01)
    slowness = slowness.reshape((1, slowness.size))
    return zetax, theta, slowness, sta_origin_x, sta_origin_y

def beamforming(seis, freqs, slowness, theta, zetax, nsources, dt, indices, new=True, matfile=None, laura=False):
    if new:
        nstat, ntimes, nsub, nfft = seis.shape
        beam = zeros((nsources, slowness.size, ntimes, nfft))
        if not DEBUG:
            widgets = ['vertical beamforming: ', pg.Percentage(), ' ', pg.Bar('#'),
                       ' ', pg.ETA()]
            pbar = pg.ProgressBar(widgets=widgets, maxval=len(indices) * ntimes).start()
            count = 0
        for ww in indices:
            FF = freqs[ww]
            for tt in xrange(ntimes):
                if not DEBUG:
                    count += 1
                    pbar.update(count)
            #for tt in [0]:
                for TT in xrange(nsub):
                #for TT in [0]:
                    Y = asmatrix(squeeze(seis[:, tt, TT, ww]))
                    YT = Y.T.copy()
                    R = dot(YT, conjugate(Y))
                    for cc in xrange(slowness.shape[1]):
                        omega = 2 * pi * FF
                        velocity = 1. / slowness[0][cc] * 1000
                        e_steer = exp(-1j * zetax * omega / velocity)
                        e_steerT = e_steer.T.copy()
                        if laura:
                            beam[:, cc, tt, ww] += sum(abs(asarray(dot(conjugate(e_steer), dot(R, e_steerT)))) ** 2, axis=1) / nsub
                        else:
                            beam[:, cc, tt, ww] += 1. / (nstat * nstat) * diag(abs(asarray(dot(conjugate(e_steer), dot(R, e_steerT)))) ** 2) / nsub
        if not DEBUG:
            pbar.finish()
        if matfile is not None:
            sio.savemat(matfile, {'beam':beam, 'theta':theta, 'slowness':slowness, 'freqs':freqs})
    else:
        beam = sio.loadmat(matfile)['beam']
    return beam

def arr_resp(nfft, dt, nstat, indices, slowness, zetax, theta, sta_origin_x, sta_origin_y,
             new=True, matfile=None, src=False, src_param=(90, 3000)):
    """
    calculate array response
    """
    if new:
        if src_param is not None:
            theta1, c1 = src_param
            zeta_x = -cos(theta1 * pi / 180.)
            zeta_y = -sin(theta1 * pi / 180.)
            zeta_src = zeta_x * sta_origin_x + zeta_y * sta_origin_y
        freqs = fftfreq(nfft, dt)
        beam = zeros((zetax.shape[0], slowness.size, nfft))
        for ww in indices:
            FF = freqs[ww]
            for cc in xrange(slowness.shape[1]):
                omega = 2 * pi * FF
                velocity = 1. / slowness[0][cc] * 1000
                e_steer = exp(-1j * zetax * omega / velocity).T
                if src:
                    e_src = exp(1j * zeta_src * omega / c1).T
                    Y = multiply(ones((zetax.shape[1], 1), 'complex128'), atleast_2d(e_src).T)
                    #Y = ones((zetax.shape[1],1),'complex128')
                    R = dot(Y, conjugate(Y).T)
                else:
                    R = ones((zetax.shape[1], zetax.shape[1]))
                #beam[:,cc,ww] = atleast_2d(sum(abs(asarray(dot(conjugate(e_steer.T),dot(R,e_steer))))**2,axis=1))

                beam[:, cc, ww] += 1. / (nstat * nstat) ** 2 * diag(abs(asarray(dot(conjugate(e_steer.T), dot(R, e_steer)))) ** 2)
                #print abs(R.max()), abs(R.min()), abs(beam[:, cc, ww].max()), abs(beam[:, cc, ww].min())
        if matfile is not None:
            sio.savemat(matfile, {'beam':beam})
        return beam
    else:
        beam = sio.loadmat(matfile)['beam']
        return beam


def polar_plot(beam, theta, freqs, slowness, dt, nfft, wtype, fout=None):
    df = dt / nfft
    periods = [6.]
    idx = [argmin(abs(freqs - 1. / p)) for p in periods]
    theta = theta[:, 0]
    slowness = slowness[0, :]
    for ind in idx:
        tre = squeeze(beam[:, :, :, ind])
        tre = tre[:, :, :].mean(axis=2)
        #tre = tre-tre.max()
        fig = figure(figsize=(6, 6))
        #ax = fig.add_subplot(1,1,1,projection='polar')
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
        cmap = cm.get_cmap('jet')
        ax.contourf((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
                    100, cmap=cmap, antialiased=True,
                    linstyles='dotted')
        ax.contour((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
                   100, cmap=cmap)
        ax.set_thetagrids([0, 45., 90., 135., 180., 225., 270., 315.],
                          labels=['90', '45', '0', '315', '270', '225', '180', '135'])
        ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5], labels=['0.1', '0.2', '0.3', '0.4', '0.5'], color='r')
        ax.grid(True)
        ax.set_title("%s %ds period" % (wtype, 1. / (freqs[ind])))
        ax.set_rmax(0.5)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=tre.min(), vmax=tre.max()))
        if fout is not None:
            savefig(fout)

def polar_plot_resp(beam, theta, slowness, dt, nfft, periods=[6.], polar=True):
    df = dt / nfft
    idx = [int(1. / (p * df)) for p in periods]
    theta = theta[:, 0]
    slowness = slowness[0, :]
    for ind in idx:
        tre = squeeze(beam[:, :, ind])
        #tre = tre-tre.max()
        fig = figure(figsize=(6, 6))
        if polar:
            ax = fig.add_subplot(1, 1, 1, projection='polar')
        else:
            ax = fig.add_subplot(1, 1, 1)
        cmap = cm.get_cmap('jet')
        if polar:
            ax.contourf((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
                        100, cmap=cmap, antialiased=True,
                        linstyles='dotted')
            ax.contour((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
                       100, cmap=cmap)
        else:
            ax.contourf(theta, slowness, tre.T, 100, cmap=cmap, antialiased=True,
                        linstyles='dotted')
            ax.contour(theta, slowness, tre.T, 100, cmap=cmap)

        if polar:
            ax.set_thetagrids([0, 45., 90., 135., 180., 225., 270., 315.],
                              labels=['90', '45', '0', '315', '270', '225', '180', '135'])
            ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5], labels=['0.1', '0.2', '0.3', '0.4', '0.5'], color='r')
            ax.set_rmax(0.5)
        else:
            ax.set_xlabel('Azimuth [degrees]')
            ax.set_ylabel('Slowness [s/km]')
        ax.grid(True)
        ax.set_title("%s %ds period" % ('Array response', 1. / (ind * df)))

def response(datdir, nprep=False, nbeam=False, doplot=True):
    """
    Calculate the theoretical array response
    """
    files = glob.glob(os.path.join(datdir, 'ft_*.*HZ.SAC'))
    if len(files) < 2:
        print "not enough files in ", datdir
        return
    if datdir.endswith('/'):
        datdir = datdir[0:-1]
    temp = os.path.basename(datdir)
    matfile1 = "%s_%s.mat" % ('prep_beam', temp)
    trZ = SacIO(files[0], headonly=True)
    sample_f = int(round(1. / trZ.delta))
    newprep = not os.path.isfile(matfile1)
    if nprep:
        newprep = True
    fseis, meanlat, meanlon, slats, slons, dt, seissmall = prep_beam(files, matfile1, onebit=False,
                                                                     nhours=1, fmax=sample_f,
                                                                     fact=sample_f,
                                                                     tempfilter=True, new=newprep)
    zetax, theta, slowness, sta_origin_x, sta_origin_y = calc_steer(slats, slons)
    matfile2 = "%s_%s.mat" % ('array_response', temp)
    nstat, ntimes, nsub, nfft = fseis.shape
    slowness = arange(0.03, 0.505, 0.005)
    slowness = slowness.reshape((1, slowness.size))
    indices = arange(nfft)
    newbeam = not os.path.isfile(matfile2)
    if nbeam:
        newbeam = True
    beam = arr_resp(nfft, dt, nstat, indices, slowness, zetax, theta, sta_origin_x, sta_origin_y,
                    new=newbeam, matfile=matfile2, src=False)
    polar_plot_resp(beam, theta, slowness, dt, nfft)
    show()

def main(datdir, files, nprep=False, nbeam=False, doplot=True, save=False, nostat=20):
    #files = glob.glob(os.path.join(datdir,'ft_*.*HZ.SAC'))
    #files = glob.glob(os.path.join(datdir,'[!^ft]*.*HZ.SAC'))
    if len(files) < nostat:
        print "not enough files in ", datdir
        print len(files), nostat
        return
    if datdir.endswith('/'):
        datdir = datdir[0:-1]
    temp = os.path.basename(datdir)
    matfile1 = "%s_%s.mat" % ('prep_beam', temp)
    trZ = SacIO(files[0], headonly=True)
    sample_f = int(round(1. / trZ.delta))
    newprep = not os.path.isfile(matfile1)
    if nprep:
        newprep = True
    fseis, freqs, meanlat, meanlon, slats, slons, dt, seissmall = prep_beam(files, matfile1,
                                                                           nhours=1, fmax=sample_f,
                                                                           new=newprep, laura=False,
                                                                           fact=sample_f,
                                                                           onebit=False,
                                                                           tempfilter=True,
                                                                           specwhite=False,
                                                                           timenorm=False)
    zetax, theta, slowness, sta_origin_x, sta_origin_y = calc_steer(slats, slons)
    matfile2 = "%s_%s.mat" % ('beam', temp)
    newbeam = not os.path.isfile(matfile2)
    if nbeam:
        newbeam = True
    nsources, ntimes, nsub, nfft = fseis.shape
    df = dt / nfft
    periods = [6.]
    periods = [6., 1. / 0.148]
    #periods = [4.,5.,6.,7.,8.,9.,10.,12.,15.,18.]
    indices = [argmin(abs(freqs - 1. / p)) for p in periods]
    beam = beamforming(fseis, freqs, slowness, theta, zetax, theta.size, dt, indices,
                              new=newbeam, matfile=matfile2, laura=False)
    if doplot:
        fout = None
        if save:
            fout = matfile2.replace('.mat', '_vertical.png')
        polar_plot(beam, theta, freqs, slowness, dt, nfft, 'rayleigh', fout=fout)
        if not save:
            show()


def proc_main():
    from optparse import OptionParser
    #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
    #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
    parser = OptionParser()
    parser.add_option("-t", "--test", dest="test", action="store_true",
                      help="Run a synthetic test using the given network layout.",
                      default=False)
    parser.add_option("-b", "--beam", dest="beam", action="store_true",
                      help="Recalculate beam.",
                      default=False)
    parser.add_option("-d", "--data", dest="data", action="store_true",
                      help="Prepare data.",
                      default=False)
    parser.add_option("-s", "--save", dest="save", action="store_true",
                      help="Save output to file.",
                      default=False)
    parser.add_option("--noplot", dest="plot", action="store_false",
                      help="Don't plot anything.",
                      default=True)
    parser.add_option("--resp", dest="aresp", action="store_true",
                      help="Calculate array response.",
                      default=False)
    parser.add_option("--nstat", dest="nstat",
                      help="Minimum number of stations.",
                      default=20)

    (opts, args) = parser.parse_args()

    flist = sys.stdin.read().split('\n')
    flist.pop()
    if opts.test:
        datdir = args[0]
        test(datdir, nbeam=opts.beam, nprep=opts.data, doplot=opts.plot)
    elif opts.aresp:
        datdir = args[0]
        response(datdir, nbeam=opts.beam, nprep=opts.data, doplot=opts.plot)
    else:
        datdir = args[0]
        main(datdir, flist, nbeam=opts.beam, nprep=opts.data, doplot=opts.plot,
             save=opts.save, nostat=int(opts.nstat))


if __name__ == '__main__':
    try:
        __IPYTHON__
    except:
        proc_main()
    else:
        print 'running in ipython'
        datdir = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/sacfiles/5Hz/2002/Jul/2002_7_12_0_0_0'
        main(datdir, nbeam=True, nprep=True, doplot=False, save=False, nostat=20)

