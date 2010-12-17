from obspy.core import UTCDateTime, Stream, read
from obspy.signal.array_analysis import sonic
from obspy.signal import cornFreq2Paz
import glob
from obspy.core.util import AttribDict

if 0:
    #
    # Load data
    fl = glob.glob('data/*.SAC')
    sts = Stream()
    for _f in fl:
        s = read(_f)
        tr = s[0]
        ts = tr.stats.starttime
        nstats =  AttribDict()
        nstats.latitude = tr.stats.sac.stla
        nstats.longitude = tr.stats.sac.stlo
        nstats.elevation = tr.stats.sac.stel/1000.
        tr.stats['coordinates'] = nstats
        tr.trim(ts,ts+3600)
        sts += s


if 0:
    #
    # Execute sonic
    kwargs = dict(
            # slowness grid: X min, X max, Y min, Y max, Slow Step
            sll_x=-.5, slm_x=.5, sll_y=-.5, slm_y=.5, sl_s=0.01,
            # sliding window propertieds
            win_len=360.0, win_frac=1.,
            # frequency properties
            frqlow=.1, frqhigh=.2, prewhiten=0,
            # restrict output
            semb_thres=-1e9, vel_thres=-1e9, verbose=True, timestamp='mlabhour',
            stime=ts, etime=ts+3600)
    outs = sonic(sts, **kwargs)


if 0:
    import matplotlib.pyplot as plt
    labels = 'rel.power abs.power baz slow'.split()

    fig = plt.figure()
    for i, lab in enumerate(labels):
        ax = fig.add_subplot(4,1,i+1)
        ax.scatter(outs[:,0], outs[:,i+1], c=outs[:,1], alpha=0.6, edgecolors='none')
        ax.set_ylabel(lab)
        ax.xaxis_date()

    fig.autofmt_xdate()
    fig.subplots_adjust(top=0.95, right=0.95, bottom=0.2, hspace=0)



if 1:
    import numpy as np
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize

    cmap = cm.hot_r
    pi = np.pi

    #
    # make output human readable, adjust backazimuth to values between 0 and 360
    t, rel_power, abs_power, baz, slow = outs.T
    baz[baz < 0.0] += 360

    # choose number of fractions in plot (desirably N is an integer!)
    N = 30
    abins = np.arange(N+1)*360./N
    sbins = np.linspace(0, .5, N+1)

    # sum rel power in bins given by abins and sbins
    hist, baz_edges, sl_edges = np.histogram2d(baz, slow,
                                                       bins=[abins, sbins], weights=rel_power)

    # transfrom to gradient
    baz_edges = baz_edges/180*np.pi

    # add polar and colorbar axes
    fig = plt.figure(figsize=(8,8))
    cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
    ax  = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)

    dh = abs(sl_edges[1] - sl_edges[0])
    dw = abs(baz_edges[1] - baz_edges[0])
    # circle through backazimut
    for i, row in enumerate(hist):
        bars = ax.bar(left=(pi/2 - (i+1)*dw)*np.ones(N),
                      height=dh*np.ones(N),
                      width=dw, bottom=dh*np.arange(N),
                      color=cmap(row/hist.max()))
        ax.set_xticks([pi/2, 0, 3./2*pi, pi])
        ax.set_xticklabels(['N', 'E' , 'S', 'W'])
        # set slowness limits
        ax.set_ylim(0, .5)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=hist.min(), vmax=hist.max()))
        
    plt.show()
