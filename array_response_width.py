#!/usr/bin/env mypython
"""
Calculate an azimuth vs. array-response width plot.
"""
import os
import sys
sys.path.append('/Volumes/GeoPhysics_05/users-data/yannik78/proc-scripts_git/beamforming')
from beamforming import arr_resp, calc_steer, polar_plot_resp
import scipy.io as sio
from obspy.core.utcdatetime import UTCDateTime
import glob
from pylab import *


def read_matfile(matfile):
    a = sio.loadmat(matfile)
    fseis = a['fseis']
    seissmall = a['seissmall']
    slats = a['slats']
    slons = a['slons']
    dt = a['dt'][0][0]
    slats = slats.reshape(slats.shape[0],)
    slons = slons.reshape(slons.shape[0],)
    meanlat = slats.mean()
    meanlon = slons.mean()
    return fseis, meanlat, meanlon, slats, slons, dt, seissmall

def array_av_resp(matfile, fout='array_response.mat', new=True):
    if new:
        maxind = 0
        nsmax = 0
        nstats = []
        fseis, meanlat, meanlon, slats, slons, dt, seissmall = read_matfile(matfile)
        zetax, theta, slowness, sta_origin_x, sta_origin_y = calc_steer(slats, slons)
        nstat, ntimes, nsub, nfft = fseis.shape
        indices = arange(nfft)
        df = dt / nfft
        periods = [4., 5., 6., 7., 8., 9., 10., 12., 15., 18.]
        periods = [6.]
        velocities = array([2000, 3000, 4000, 5000, 6000])
        azimuths = theta.T[0]
        idx = [int(1. / (p * df)) for p in periods]
        errs = zeros((azimuths.size, velocities.size, len(idx)))
        for a, _az in enumerate(azimuths)[0:1]:
            for c, _c in enumerate(velocities[1:2]):
                print _az, _c
                beam = arr_resp(nfft, dt, nstat, indices, slowness, zetax, theta,
                                sta_origin_x, sta_origin_y,
                                new=True, matfile=None, src=True, src_param=(_az, _c))
                polar_plot_resp(beam, theta, slowness, dt, nfft, polar=False, periods=[7.])
                for i, ind in enumerate(idx):
                    tre = squeeze(beam[:, :, ind])
                    _tmp, _sl = unravel_index(tre.argmax(), tre.shape)
                    hcrv = tre[:, _sl]
                    figure()
                    plot(theta, hcrv)
                    idx2 = where(hcrv > hcrv.max() / 2.)
                    minerr = theta[idx2[0][0], 0]
                    maxerr = theta[idx2[0][-1], 0]
                    err = maxerr - minerr
                    if err > 180.:
                        err = 360 - err
                    errs[a, c, i] = err
        sio.savemat(fout, {'az':azimuths, 'slowness':slowness,
                          'periods':array(periods), 'error':errs})
    else:
        a = sio.loadmat(fout)
        azimuths = a['az'].T[0]
        slowness = a['slowness'].T[0]
        periods = a['periods'].T[0]
        errs = a['error']
    return azimuths, slowness, periods, errs

if __name__ == '__main__':
    #dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams'
    #dirn = '/Volumes/Wanaka_01/yannik/start/beamforming'
    taranaki_max = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/prep_beam_2002_8_7_0_0_0.mat'
    taranaki_min = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/prep_beam_2002_3_23_0_0_0.mat'
    start_max = '/Volumes/Wanaka_01/yannik/start/beamforming/prep_beam_2001_2_24_0_0_0.mat'
    start_min = '/Volumes/Wanaka_01/yannik/start/beamforming/prep_beam_2001_1_19_0_0_0.mat'
    matfile = taranaki_max
    #matfile = start_max
    azimuths, slowness, periods, errs = array_av_resp(matfile,
                                                   fout='array_response_taranaki_max.mat',
                                                   new=True)
    plot(azimuths, errs[:, 0, 0])
    show()

    #polar_plot_resp(beam,theta,slowness,dt,nfft,polar=False,periods=[7.])
    if 0:
        for matfile in matfiles[0:1]:
            fseis, meanlat, meanlon, slats, slons, dt, seissmall = read_matfile(matfile)
            zetax, theta, slowness, sta_origin_x, sta_origin_y = calc_steer(slats, slons)
            nstat, ntimes, nsub, nfft = fseis.shape
            indices = arange(nfft)
            slowness = arange(0.03, 0.505, 0.01)
            slowness = arange(0.00, 0.505, 0.01)
            slowness = slowness.reshape((1, slowness.size))
            a = os.path.basename(matfile).split('_')
            year = a[-6]
            month = a[-5]
            day = a[-4]
            matfile2 = os.path.join(dirout, 'array_resp_%s_%s_%s_0_0_0.mat' % (year, month, day))
            #beam = arr_resp(nfft,dt,nstat,indices,slowness,zetax,theta,sta_origin_x,sta_origin_y,
            #                new=False,matfile=matfile2,src=True,src_param=(210,3000))
            #polar_plot_resp(beam,theta,slowness,dt,nfft,polar=False,periods=[7.])
            if 1:
                df = dt / nfft
                periods = [4., 5., 6., 7., 8., 9., 10., 12., 15., 18.]
                #periods = [7.]
                #velocities = [2000,3000,4000,5000,6000]
                slowness = arange(0.125, 0.51, 0.01)
                velocities = 1. / slowness * 1000.
                slowness = slowness.reshape((1, slowness.size))
                azimuths = theta.T[0]
                idx = [int(1. / (p * df)) for p in periods]
                errs = empty((azimuths.size, velocities.size, len(idx)))
                for a, _az in enumerate(azimuths):
                    for c, _c in enumerate(velocities):
                        print _az, _c
                        beam = arr_resp(nfft, dt, nstat, indices, slowness, zetax, theta,
                                        sta_origin_x, sta_origin_y,
                                        new=True, matfile=None, src=True, src_param=(_az, _c))
                        for i, ind in enumerate(idx):
                            tre = squeeze(beam[:, :, ind])
                            _tmp, _sl = unravel_index(tre.argmax(), tre.shape)
                            hcrv = tre[:, _sl]
                            idx2 = where(hcrv > hcrv.max() / 2.)
                            minerr = theta[idx2[0][0], 0]
                            maxerr = theta[idx2[0][-1], 0]
                            errs[a, c, i] = maxerr - minerr
                matfile = 'array_response_taranaki.mat'
                sio.savemat(matfile, {'az':azimuths, 'slowness':slowness, 'periods':array(periods),
                                     'error':errs})
            #polar_plot_resp(beam,theta,slowness,dt,nfft,polar=False,periods=[7.])
            #plot(theta[2:-2],errs[2:-2],label=str(periods[i])+' s')
            #xlabel('Azimuth [degrees]')
            #ylim(0,30)

            if 0:
                df = dt / nfft
                periods = [7., 9., 12., 15., 18.]
                periods = [7.]
                idx = [int(1. / (p * df)) for p in periods]
                for i, ind in enumerate(idx):
                    tre = squeeze(beam[:, :, ind])
                    _az, _sl = unravel_index(tre.argmax(), tre.shape)
                    hcrv = tre[:, _sl]
                    idx = where(hcrv > hcrv.max() / 2.)
                    minerr = theta[idx[0][0], 0]
                    maxerr = theta[idx[0][-1], 0]
                    print maxerr - minerr
                    plot(theta, hcrv, label=str(periods[i]) + ' s')
                    hlines([hcrv.max() / 2], minerr, maxerr)
                legend(loc='upper left')
                xlabel('Azimuth [degrees]')

            if 0:
                df = dt / nfft
                periods = [7.]
                idx = [int(1. / (p * df)) for p in periods]
                for i, ind in enumerate(idx):
                    for c in [2000, 3000, 4000, 5000, 6000]:
                        print c
                        beam = arr_resp(nfft, dt, nstat, indices, slowness, zetax,
                                        theta, sta_origin_x, sta_origin_y,
                                        new=True, matfile=matfile2, src=True, src_param=(210, c))
                        tre = squeeze(beam[:, :, ind])
                        _az, _sl = unravel_index(tre.argmax(), tre.shape)
                        plot(theta, tre[:, _sl], label=str(c) + ' m/s')
                legend(loc='upper left')
                xlabel('Azimuth [degrees]')


            if 0:
                df = dt / nfft
                periods = [4., 5., 6., 7., 8., 9., 10., 12., 15., 18.]
                #periods = [7.]
                idx = [int(1. / (p * df)) for p in periods]
                for ind in idx:
                    errs = []
                    for az in xrange(36):
                        tre = squeeze(beam[:, :, ind])
                        crv = r_[tre[az + 36, :][::-1], tre[az, :]]
                        hcrv = (tre[az + 36, :] + tre[az, :]) / 2
                        idx = where(hcrv < hcrv.max() / 2)
                        minerr = slowness[0, idx[0][0]]
                        #plot(r_[-slowness[0,:][::-1],slowness[0,:]],crv)
                        #hlines([hcrv.max()/2],-minerr,minerr)
                        wdth = 2 * minerr
                        errs.append(wdth)
                    plot(theta[0:36], errs)
