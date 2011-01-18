#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/************************************************************************/
/* numpy fftpack_lite definitions, must be linked with fftpack_lite.so  */
/* e.g. link with -L/path/to/fftpack_lite -l:fftpack_lite.so            */       
/* numpy has no floating point transfrom compiled, therefore we need to */
/* cast the result to double                                            */    
/************************************************************************/
extern void rfftf(int N, double* data, double* wrk);
extern void rffti(int N, double* wrk);

int beam(double *traces, int stride1, int stride2, int stride3, int stride4, 
	 int nfft, int nstat, int ntimes, int nsub, int digfreq, double flow, 
	 double fhigh){

  double *fftpack_work = 0;
  int wlow,whigh,i,j,k,l;
  double df;
  double **window;
  fftpack_work = (double *)calloc((2*nfft+15), sizeof(double));
  rffti(nfft, fftpack_work);
  df = digfreq/(float)nfft;
  wlow = (int)(flow/df+0.5);
  if (wlow < 1) {
    /*******************************************************/
    /* never use spectral value at 0 -> this is the offset */
    /*******************************************************/
    wlow = 1;
  }
  whigh = (int)(fhigh/df+0.5);
  if (whigh>(nfft/2-1)) {
    /***************************************************/
    /* we avoid using values next to nyquist frequency */
    /***************************************************/
    whigh = nfft/2-1;
  }

  /****************************************************************************/
  /* first we need the fft'ed window of traces for the stations of this group */
  /****************************************************************************/
  for(k=0;l<ntimes;k++){
    for(l=0;l<nsub;l++){
      window = (double **)calloc(nstat, sizeof(double *));
      for(j=0;j<nstat;j++) {
	window[j] = (double *)calloc(nfft+1, sizeof(double));
	for(i=0;i<nfft;i++){
	  window[j][i] = traces[j*stride1/sizeof(double)+k*stride2/sizeof(double)+l*stride3/sizeof(double)+i*stride4/sizeof(double)];
	}
	rfftf(nfft, window[j]+1, fftpack_work);
	window[j][0] = window[j][1];
	window[j][1] = 0.0;
      }
    for (j=0;j<nstat;j++) {
        free((void *)window[j]);
    }
    free((void *)window);
    }
  }

  return 1;
}
