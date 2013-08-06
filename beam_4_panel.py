#!/usr/bin/env python
"""
Plot output from beamforming averaged over the
whole deployment time.

Created on Sep 12, 2012

@author: behry
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.io as sio
import numpy as np
import os
from matplotlib.colorbar import ColorbarBase
from matplotlib import rcParams
from matplotlib.colors import Normalize
from mpl_toolkits.basemap import cm as cmb



def polar_plot(ax, beam, theta, slowness, freqs, period, wtype, cmap):
    ind = np.argmin(np.abs(freqs - 1. / period))
    tre = np.squeeze(beam[:, :, :, ind])
    tre = tre.mean(axis=2)
    tre = 10 * np.log10(np.abs(tre))
    max = tre.max()
    min = tre.min()
    #tre = (tre - min) / (max - min)
    tre -= tre.max()
    #print tre.min(), tre.max()
    inds = np.where(slowness > 0.2)
    ax.contourf((theta[::-1] + 90.) * np.pi / 180., slowness[inds], tre[:, inds[0]].T,
                100, cmap=cmap, antialiased=True,
                linstyles='dotted', norm=Normalize(vmin= -4, vmax=0))
    ax.contour((theta[::-1] + 90.) * np.pi / 180., slowness[inds], tre[:, inds[0]].T,
                100, cmap=cmap, norm=Normalize(vmin= -4, vmax=0))
    ax.set_thetagrids([0, 45., 90., 135., 180., 225., 270., 315.],
                      labels=['90', '45', '0', '315', '270', '225', '180', '135'],
                      frac=1.11)
    s_cutoff = 1. / (4.5 - (25. / 3.) * 1. / period)
    ax.plot((theta[::-1] + 90.) * np.pi / 180., np.repeat(s_cutoff, theta.size), 'k--')
    print s_cutoff, period

    if wtype == 'z':
        if period < 10.:
            ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5],
                          labels=['0.1', '0.2', '0.3', '0.4', '0.5'], color='r')
        else:
            ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5],
                          labels=['0.1', '0.2', '0.3', '0.4', '0.5'], color='k')
        ax.set_rmax(0.5)
        ax.set_title('%d s, vertical' % (period), va='bottom')
    elif wtype == 'h':
        ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                      labels=['0.1', '0.2', '0.3', '0.4', '0.5', '0.6'], color='k')
        ax.set_rmax(0.6)
        ax.set_title('%d s, transverse' % (period), va='bottom')
    ax.grid(True)


if __name__ == '__main__':
    rcParams['figure.subplot.left'] = 0.05
    rcParams['figure.subplot.right'] = 0.95
    rcParams['figure.subplot.top'] = 0.92
    rcParams['figure.subplot.bottom'] = 0.15
    rcParams['figure.subplot.wspace'] = 0.25
    rcParams['figure.subplot.hspace'] = 0.32

    datadir = '/home/behry/uni/data/beamformer_movie'
    savedir = '/home/behry/uni/noise_dir_paper/pix'
    fout = os.path.join(savedir, 'beamformer_4_panel.pdf')
    beamfnh = 'average_beam_horizontal_2nd_run.mat'
    beamfnz = 'average_beam_vertical.mat'
    a = sio.loadmat(os.path.join(datadir, beamfnz))
    avbeam = a['avbeam']
    theta = a['theta'].T[0]
    slowness = a['slowness'].T[0]
    freqs = a['freqs']

    b = sio.loadmat(os.path.join(datadir, beamfnh))
    avbeamh = b['avbeam']
    thetah = b['theta'].T[0]
    slownessh = b['slowness'].T[0]
    freqsh = b['freqs']

    fig = plt.figure(figsize=(8, 8))
    cmap = cm.get_cmap('gist_rainbow_r')
    cmap = cm.get_cmap('jet')
    period = 6.
    ax = fig.add_subplot(2, 2, 1, projection='polar')
    polar_plot(ax, avbeam, theta, slowness, freqs, period, 'z', cmap)

    period = 12.
    ax = fig.add_subplot(2, 2, 3, projection='polar')
    polar_plot(ax, avbeam, theta, slowness, freqs, period, 'z', cmap)

    period = 6.
    ax = fig.add_subplot(2, 2, 2, projection='polar')
    polar_plot(ax, avbeamh, thetah, slownessh, freqsh, period, 'h', cmap)

    period = 12.
    ax = fig.add_subplot(2, 2, 4, projection='polar')
    polar_plot(ax, avbeamh, thetah, slownessh, freqsh, period, 'h', cmap)

    fig.text(0.05, 0.95, 'A', weight='bold')
    fig.text(0.55, 0.95, 'B', weight='bold')
    fig.text(0.55, 0.51, 'D', weight='bold')
    fig.text(0.05, 0.51, 'C', weight='bold')
    # colorbar
    if True:
        cax = fig.add_axes([0.25, 0.06, 0.5, 0.03])
        #cb = ColorbarBase(cax, cmap=cmap, norm=Normalize(vmin=0, vmax=1),
        #             orientation='horizontal', ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
        cb = ColorbarBase(cax, cmap=cmap, norm=Normalize(vmin= -4, vmax=0),
                     orientation='horizontal', ticks=[-4, -3, -2, -1, 0.])
        cb.set_label('dB')
    plt.savefig(fout)
    #plt.show()
