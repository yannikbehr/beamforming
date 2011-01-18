
FFT = $(shell /usr/local/python2/bin/python -c "import numpy, os; print os.path.dirname(numpy.fft.__file__)")
CC = gcc

beam: beam_c.c
	$(CC) -shared -fPIC -o beam_c.so $(FFT)/fftpack_lite.so beam_c.c 