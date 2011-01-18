#!/usr/local/python2/bin/python

import ctypes as C
from pylab import *


lib = C.cdll.LoadLibrary('./test_ndim_array.so')

if 0:
    ### test 2dimension arrays
    a = arange(6.).reshape((2,3))
    dim1,dim2 = a.shape
    ptr = C.POINTER(C.c_double)
    data = (ptr*dim1)(*[row.ctypes.data_as(ptr) for row in a])
    errcode = lib.test_2dim_array(C.byref(data))
    print a[0][0]

if 0:
    ### test 4 dimensionl arrays
    a = random(192.).reshape((2,4,4,6))
    dim1,dim2,dim3,dim4 = a.shape
    ptr = C.c_void_p
    voids1 = []
    print a[1][3][2][3]
    for i in xrange(dim1):
        voids2 = []
        for j in xrange(dim2):
            row = a[i][j]
            p = (ptr * dim3)(*[col.ctypes.data_as(ptr) for col in row])
            voids2.append(C.cast(p, C.c_void_p))
        p2 = (ptr*dim2)(*voids2)
        voids1.append(C.cast(p2,C.c_void_p))
    data = (ptr*dim1)(*voids1)
    errcode = lib.test_4dim_array(C.byref(data))
    print a[1][1][2][3]

if 1:
    ### test complex arrays
    a = ones((3,3,4,2),'complex')
    a[2,1,2,0] = complex(2,12)
    fnc = lib.test_complex_array
    fnc.argtypes = [\
        ctypeslib.ndpointer(complex,ndim=4,flags='aligned, contiguous'),
        C.POINTER(C.c_long),
        C.POINTER(C.c_long)]
    errcode = fnc(a,a.ctypes.strides,a.ctypes.shape)
