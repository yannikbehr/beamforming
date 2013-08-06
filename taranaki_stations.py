#!/usr/bin/env python
"""
Plot the number of stations available per day for the
Taranaki dataset.

Created on Sep 10, 2012

@author: behry
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime
from matplotlib.dates import date2num, num2date

fin = '/home/behry/uni/data/taranaki_dataset/number_of_stations_per_day.txt'
dates, stats = np.loadtxt(fin, converters={0:lambda x: date2num(UTCDateTime(x))}, unpack=True)
print num2date(dates[stats.argmin()])
idx = np.argsort(dates)
print stats.max()
plt.plot_date(dates[idx], stats[idx], 'k-')
plt.ylabel('Number of stations')
#plt.savefig('number_of_stations_taranaki.png')
plt.show()
