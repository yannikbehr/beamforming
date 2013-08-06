#!/usr/bin/env mypython
"""
Calculate an azimuth vs. array-response width plot.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
from matplotlib.colors import Normalize
import numpy as np
import os
import sys
sys.path.append('/Volumes/GeoPhysics_05/users-data/yannik78/proc-scripts_git/beamforming')
from beamforming import arr_resp, calc_steer
import scipy.io as sio
from obspy.core.utcdatetime import UTCDateTime
import glob

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

def polar_plot_resp(ax, beam, theta, slowness, dt, nfft, nstat, periods=[6.]):
    df = dt / nfft
    idx = [int(1. / (p * df)) for p in periods]
    theta = theta[:, 0]
    slowness = slowness[0, :]
    for ind in idx:
        tre = np.squeeze(beam[:, :, ind])
        cmap = cm.get_cmap('jet')
        ax.contourf((theta[::-1] + 90.) * np.pi / 180., slowness, tre.T,
                    100, cmap=cmap, antialiased=True,
                    linstyles='dotted', norm=Normalize(vmin=0.0, vmax=0.14))
        ax.contour((theta[::-1] + 90.) * np.pi / 180., slowness, tre.T,
                   100, cmap=cmap)
        ax.set_thetagrids([0, 45., 90., 135., 180., 225., 270., 315.],
                          labels=['90', '45', '0', '315', '270', '225', '180', '135'])
        ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5], labels=['0.1', '0.2', '0.3', '0.4', '0.5'], color='r')
        ax.set_rmax(0.5)
        ax.set_title('%d s, %d stations' % (periods[0], nstat), va='bottom')
        ax.grid(True)

def array_av_resp(matfile, period):
    maxind = 0
    nsmax = 0
    fseis, meanlat, meanlon, slats, slons, dt, seissmall = read_matfile(matfile)
    zetax, theta, slowness, sta_origin_x, sta_origin_y = calc_steer(slats, slons)
    nstat, ntimes, nsub, nfft = fseis.shape
    indices = np.arange(nfft)
    df = dt / nfft
    if period > 6.:
        velocities = 1000. / np.array([0.2, 0.4])
    else:
        velocities = 1000. / np.array([0.2, 0.3, 0.4])
    azimuths = theta.T[0]
    azimuths = np.array([0, 90, 180, 270])
    df = dt / nfft
    indices = [int(1. / (period * df))]
    beamall = None
    for a, _az in enumerate(azimuths):
        for c, _c in enumerate(velocities):
            print _az, _c
            beam = arr_resp(nfft, dt, nstat, indices, slowness, zetax, theta,
                            sta_origin_x, sta_origin_y,
                            new=True, matfile=None, src=True, src_param=(_az, _c))
            if beamall is None:
                beamall = beam
            else:
                beamall += beam
    beamall /= azimuths.size * velocities.size
    print '--> ', beamall.max(), beamall.min()
    print '-->', azimuths.size * velocities.size
    return beamall, theta, slowness, dt, nfft, nstat

if __name__ == '__main__':
    datadir = '/home/behry/uni/data/beamformer_movie'
    savedir = '/home/behry/uni/noise_dir_paper/pix'
    taranaki_max = os.path.join(datadir, 'prep_beam_2002_8_7_0_0_0.mat')
    taranaki_max = os.path.join(datadir, 'prep_beam_2002_7_29_0_0_0.mat')
    taranaki_min = os.path.join(datadir, 'prep_beam_2002_3_23_0_0_0.mat')
    taranaki_min = os.path.join(datadir, 'prep_beam_2002_3_17_0_0_0.mat')
    rcParams['figure.subplot.left'] = 0.05
    rcParams['figure.subplot.right'] = 0.95
    rcParams['figure.subplot.top'] = 0.92
    rcParams['figure.subplot.bottom'] = 0.05
    rcParams['figure.subplot.wspace'] = 0.25
    rcParams['figure.subplot.hspace'] = 0.3

    matfile = taranaki_max
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(2, 2, 1, projection='polar')
    period = 6.
    beamall, theta, slowness, dt, nfft, nstat = array_av_resp(matfile, period)
    polar_plot_resp(ax, beamall, theta, slowness, dt, nfft, nstat,
                    periods=[period])
    if True:
        ax = fig.add_subplot(2, 2, 2, projection='polar')
        matfile = taranaki_min
        beamall, theta, slowness, dt, nfft, nstat = array_av_resp(matfile, period)
        polar_plot_resp(ax, beamall, theta, slowness, dt, nfft, nstat,
                        periods=[period])

        ax = fig.add_subplot(2, 2, 3, projection='polar')
        matfile = taranaki_max
        period = 12.
        beamall, theta, slowness, dt, nfft, nstat = array_av_resp(matfile, period)
        polar_plot_resp(ax, beamall, theta, slowness, dt, nfft, nstat,
                        periods=[period])

        ax = fig.add_subplot(2, 2, 4, projection='polar')
        matfile = taranaki_min
        beamall, theta, slowness, dt, nfft, nstat = array_av_resp(matfile, period)
        polar_plot_resp(ax, beamall, theta, slowness, dt, nfft, nstat,
                        periods=[period])

    fig.text(0.05, 0.95, 'A', weight='bold')
    fig.text(0.55, 0.95, 'B', weight='bold')
    fig.text(0.55, 0.46, 'D', weight='bold')
    fig.text(0.05, 0.46, 'C', weight='bold')
    fout = os.path.join(savedir, 'array_response_panels.pdf')
    plt.savefig(fout)
    #plt.show()
