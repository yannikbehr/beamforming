
FFT = $(shell /usr/local/python2/bin/python -c "import numpy, os; print os.path.dirname(numpy.fft.__file__)")
CC = gcc
ICC = icc

beam: beam_c.c
	$(CC) -funroll-all-loops -Wall -shared -fPIC -o beam_c.so $(FFT)/fftpack_lite.so beam_c.c 

beamicc:
	$(ICC) -O2 -fpic -shared -obeam_c.so $(FFT)/fftpack_lite.so beam_c.c 