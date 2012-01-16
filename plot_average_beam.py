#!/usr/bin/env mypython
"""
Plot output from beamforming averaged either over a month or over the
whole deployment time.
"""
import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io as sio
import sys
import glob
import os
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib import rcParams
from obspy.core.utcdatetime import UTCDateTime
from optparse import OptionParser

def mycmp(beama, beamb):
    tmp = os.path.basename(beama).split('_')
    year = int(tmp[-6])
    month = int(tmp[-5])
    day = int(tmp[-4])
    datea = UTCDateTime(year,month,day,0,0,0,0).getTimeStamp()
    tmp = os.path.basename(beamb).split('_')
    year = int(tmp[-6])
    month = int(tmp[-5])
    day = int(tmp[-4])
    dateb = UTCDateTime(year,month,day,0,0,0,0).getTimeStamp()
    if datea > dateb: return 1
    if datea < dateb: return -1
    if datea == dateb: return 0

def getmonth(beam):
    tmp = os.path.basename(beam).split('_')
    year = int(tmp[-6])
    month = int(tmp[-5])
    day = int(tmp[-4])
    return month

def polar_plot(beam,theta,slowness,freqs,wtype,fout=None):
    #periods = [4.,5.,6.,7.,8.,9.,10.]
    period = 6.
    ind = argmin(abs(freqs - 1./period))
    #tre = squeeze(beam[:,:,:,ind])
    #tre = tre.mean(axis=2)
    tre = beam
    tre = 10*log10(abs(tre))
    tre = tre-tre.max()
    tre = tre.reshape(1,tre.shape[0])
    tre = np.repeat(tre,6,axis=0)
    fig = figure(figsize=(6,6))
    cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
    ax  = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
    cmap = cm.get_cmap('jet')
    ax.contourf((theta[::-1]+90.)*pi/180.,np.arange(0.55,0.61,0.01),tre,
                100,cmap=cmap,antialiased=True,
                linstyles='dotted')
    ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                      labels=['90','45','0','315','270','225','180','135'])
    ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
    ax.grid(True)
    #ax.set_title("%s %ds period"%(wtype,1./(ind*df)))
    ax.set_rmax(0.6)
    ax.set_rmin(0.55)
    ColorbarBase(cax, cmap=cmap,
                 norm=Normalize(vmin=tre.min(), vmax=tre.max()))
    if fout is not None:
        savefig(fout)
    show()

def polar_plot_panel(beam,theta,slowness,freqs,wtype,fout=None):
    rcParams['figure.subplot.left'] = 0.05
    rcParams['figure.subplot.right'] = 0.95
    rcParams['figure.subplot.top'] = 0.95
    rcParams['figure.subplot.bottom'] = 0.05
    rcParams['figure.subplot.wspace'] = 0.29
    periods = [6.]
    periods = [5.,6.,7.,8.,9.,10.,12.,15.,18.]
    idx = [argmin(abs(freqs - 1./p)) for p in periods]
    cnt = 1
    cmap = cm.get_cmap('jet')
    fig = figure(figsize=(8,8))
    for ind in idx:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre.mean(axis=2)
        tre = 10*log10(tre)
        tre = tre-tre.max()
        ax = fig.add_subplot(3,3,cnt,projection='polar')
        inds = where(slowness > 0.2)
        ax.contourf((theta[::-1]+90.)*pi/180.,slowness[inds],tre[:,inds[0]].T,
                    100,cmap=cmap,antialiased=True,
                    linstyles='dotted')
        ax.contour((theta[::-1]+90.)*pi/180.,slowness[inds],tre[:,inds[0]].T,
                   100,cmap=cmap)
        ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                          labels=['90','45','0','315','270','225','180','135'])
        ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
        v_cutoff = 4.5 - (25./3.) * freqs[ind]
        ax.plot(theta*pi/180.,ones(theta.size)*1./v_cutoff,'k')
        ax.grid(True)
        ax.set_title("%s %ds period"%(wtype,periods[cnt-1]))
        ax.set_rmax(0.6)
        cnt += 1
    if fout is not None:
        savefig(fout)

def monthly_average(fl,months,comp,new=True):
    fl.sort(cmp=mycmp)
    if new:
        for month in months.keys():
            print month
            avbeam = None
            cnt = 0
            for _f in fl:
                if getmonth(_f) == months[month]:
                    print _f
                    if comp == 'transverse':
                        beam = sio.loadmat(_f)['beamt']
                    elif comp == 'vertical':
                        beam = sio.loadmat(_f)['beam']
                    else:
                        print 'comp: ',comp
                        print 'comp has to be transverse or vertical'
                        return
                    
                    if avbeam is None:
                        avbeam = beam
                    else:
                        avbeam += beam
                    cnt += 1
            if avbeam is None: continue
            avbeam /= cnt
            sio.savemat(os.path.join(dirn,'average_beam_%s_%s.mat'%(comp,month)),{'avbeam':avbeam})
            dt = 1.0
            theta= arange(0,365,5)
            nfft = 128
            slowness = arange(0.125,0.51,0.01)
            fout = os.path.join(dirn,'average_beam_%s_%s.png'%(comp,month))
            polar_plot_panel(avbeam,theta,slowness,dt,nfft,'%s %s'%(comp,month),fout=fout)
    else:
        for month in months.keys():
            matfile = os.path.join(dirn,'average_beam_%s_%s.mat'%(comp,month))
            print matfile
            avbeam = sio.loadmat(matfile)['avbeam']
            dt = 1.0
            theta= arange(0,365,5)
            nfft = 128
            slowness = arange(0.125,0.51,0.01)
            fout = os.path.join(dirn,'average_beam_%s_%s.png'%(comp,month))
            polar_plot_panel(avbeam,theta,slowness,dt,nfft,'%s %s'%(comp,month),fout=fout)

def average(fl,comp,new=True,single_vel=True,fout='./average_beam.mat'):
    if not os.path.isfile(fout):
        new = True
    if new:
        avbeam = None
        for _f in fl:
            print _f
            a = sio.loadmat(_f)
            if comp == 'transverse':
                beam = a['beamt']
            elif comp == 'vertical':
                beam = a['beam']
            else:
                print 'comp: ',comp
                print 'comp has to be transverse or vertical'
                return
            theta= a['theta'].T[0]
            slowness = a['slowness'][0]
            freqs = a['freqs'].T[0]
            if single_vel:
                vel = 1./slowness
                vel1 = 2.47
                p = 6.
                idv = argmin(abs(vel - vel1))
                idf = argmin(abs(freqs - 1./p))
                beam1 = squeeze(beam[:,idv,:,idf])
                b1meana = beam1.mean(axis=1)
                if avbeam is None:
                    avbeam = b1meana
                else:
                    avbeam += b1meana
            else:
                if avbeam is None:
                    avbeam = beam
                else:
                    avbeam += beam
        avbeam /= len(fl)
        sio.savemat(fout,{'avbeam':avbeam,'theta':theta,'slowness':slowness,'freqs':freqs})
    else:
        a = sio.loadmat(fout)
        avbeam = a['avbeam']
        theta= a['theta'].T[0]
        slowness = a['slowness'].T[0]
        freqs = a['freqs']
    if single_vel:
        polar_plot(avbeam,theta,slowness,freqs,'rayleigh')
    else:
        polar_plot_panel(avbeam,theta,slowness,freqs,'%s '%(comp),fout=fout.replace('.mat','.png'))


def main():
    usage = "usage: files | %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--dep",dest="dep",
                      help="Chose between 'Taranaki' and 'START'",
                      default='Taranaki')
    parser.add_option("-o","--fout",dest="fout",
                      help="Output filename",
                      default='./average_beam.mat')
    parser.add_option("-a","--average",dest="average",
                      help="Choose between 'monthly' or 'all'",
                      default='all')
    parser.add_option("-c","--comp",dest="comp",
                      help="Choose between 'vertical' or 'transverse'",
                      default='vertical')
    parser.add_option("-n","--new",dest="new",action='store_true',
                      help="Recalculate average beam",
                      default=False)
    parser.add_option("-l","--laura",dest="laura",action='store_true',
                      help="Calculate average beam as done by Laura Brooks",
                      default=False)
    (opts,args) = parser.parse_args()
    if len(args) != 0:
        parser.error("incorrect number of arguments")
    if opts.dep == 'Taranaki':
        months = {'march':3,'april':4,'may':5,'june':6,'july':7,'august':8,'september':9}
    if opts.dep == 'START':
        months = {'january':1,'february':2,'march':3,'april':4,'may':5,'june':6}

    fl = sys.stdin.read().split('\n')
    fl.pop()

    if opts.average =='monthly':
        monthly_average(fl,months,opts.comp,new=opts.new)
    if opts.average =='all':
        average(fl,opts.comp,fout=opts.fout,new=opts.new,single_vel=opts.laura)


if __name__ == '__main__':
    main()
### Bash lines to print the output of monthly_average
# for i in average*transverse*png;do convert $i -resize 75% `echo $i|cut -d. -f1`_small.png;done
# for i in average*transverse*small.png; do lpr -PCO505 $i;done
