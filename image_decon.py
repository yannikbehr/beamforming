#!/usr/bin/env mypython

"""
test different ways of image convolution and deconvolution
"""

from pylab import *

def conv(mat,kern):
    ### zero_padding
    row,col = mat.shape
    kr, kc = kern.shape
    m = zeros((row+2,col+2))
    m[1:-1,1:-1] = mat
    nrow, ncol = m.shape
    m[where(m == 0)] = 127
    #imshow(m)
    res = array([])
    for _r in xrange(nrow-kr+1):
        for _c in xrange(ncol-kc+1):
            tmat = m[_r:_r+kr,_c:_c+kc]
            res = append(res,sum(tmat*kern))
    return res.reshape(row,col)

def conv_fft(mat,kern):
    ### zero padding
    row,col = mat.shape
    nfftx = 2**(int(log(col)/log(2))+1)
    nffty = 2**(int(log(row)/log(2))+1)
    padx = nfftx/2
    pady = nffty/2
    print row,col,nfftx,nffty, padx, pady
    m = zeros((2*nfftx,2*nffty))
    m[pady:pady+row,padx:padx+col] = mat
    print m[:,8], mat[:,0]

a = random_integers(0,255,size=64).reshape((8,8))
#kern = ones((3,3))*1/9
kern = array([[0,0,0],[0,1,0],[0,0,-1]])
kern = rot90(rot90(kern))
kern = array([[2,0,0],[0,-1,0],[0,0,-1]])
#a = imread('stinkbug.png')
#b = (255*a).mean(axis=2)
conv_fft(a,kern)
#res = conv(b,kern)
#figure()
#imshow(res,cmap='gist_gray')
#imshow(a)


