#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
typedef struct {double real; double imag;} cdouble;

/************************************************************************/
/* numpy fftpack_lite definitions, must be linked with fftpack_lite.so  */
/* e.g. link with -L/path/to/fftpack_lite -l:fftpack_lite.so            */       
/* numpy has no floating point transfrom compiled, therefore we need to */
/* cast the result to double                                            */    
/************************************************************************/
extern void rfftf(int N, double* data, double* wrk);
extern void rffti(int N, double* wrk);

int beam_fft(double ****traces, double ***beam, double **zetax,
	     double *slowness, int nslow, int nsources, int nfft, 
	     int nstat, int ntimes, int nsub, int digfreq, double flow, double fhigh){

  double *fftpack_work = 0;
  int wlow,whigh,i,j,m,n,k,l;
  double df;
  double **window;
  double re;
  double im;
  double sumre;
  double sumim;
  int w;
  double vel;
  double wtau;

  printf("%g\n",traces[3][4][6][100]);
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
  
  window = (double **)calloc(nstat, sizeof(double *));
  for(j=0;j<nstat;j++) {
    window[j] = (double *)calloc(nfft+1, sizeof(double));
  }
  for(m=0;m<ntimes;m++){
    for(n=0;n<nsub;n++){
      printf("ntimes:%d nsub:%d\n",m,n); 
      /****************************************************************************/
      /* first we need the fft'ed window of traces for the stations of this group */
      /****************************************************************************/
	//memcpy((window[j]+1),traces[j][k][l][i],nfft*sizeof(double));
      for(j=0;j<nstat;j++){
	for(i=0;i<nfft;i++){
	  window[j][i] = traces[j][m][n][i];
	}
	rfftf(nfft, window[j]+1, fftpack_work);
	window[j][0] = window[j][1];
	window[j][1] = 0.0;
      }
    

      /*************************************************************/
      /* we start with loop over angular frequency
      /*************************************************************/
      for (w=wlow;w<=whigh;w++) {
	/***********************************/
	/* now we loop over slowness index */
	/***********************************/
	for (k=0;k<nslow;k++) {
	  vel = 1./slowness[k]*1000.;
	  /************************************/
	  /* now we loop over azimuth index */
	  /************************************/
	  for (l=0;l<nsources;l++) {
	    /********************************************/
	    /* this is the loop over the stations group */
	    /********************************************/
	    sumre = sumim = 0.;
	    for (j=0;j<nstat;j++) {
	      wtau = (float) (-2.*M_PI*df*(float)w*zetax[l][j]/vel);
	      re = window[j][2*w];
	      im = window[j][2*w+1];
	      sumre += (float) (re*cos(wtau)-im*sin(wtau));
	      sumim += (float) (im*cos(wtau)+re*sin(wtau));
	    }
	    beam[w][l][k] += (sumre*sumre+sumim*sumim)/nsub;
	  }
	}
      }
    }
    //beam[w][l][k] /= ntimes;
  }
  for (j=0;j<nstat;j++) {
    free((void *)window[j]);
  }
  free((void *)window);
  return 1;
}


int beam(cdouble *traces, long *stride, long *shape, double ***beam, double **zetax,
	 double *slowness, int nslow, int nsources, int nfft, 
	 int nstat, int ntimes, int nsub, int digfreq){

  double *fftpack_work = 0;
  int wlow,whigh,i,j,m,n,k,l;
  double df;
  double **window;
  double re;
  double im;
  double sumre;
  double sumim;
  int w;
  double vel;
  double wtau;
  double flow,fhigh;
  int M,N,K,L,S0,S1,S2,S3;

  //flow = freq_int[0];
  //fhigh = freq_int[1];
  flow = 0.02;
  fhigh = 0.4;
  M = shape[0]; N = shape[1]; K = shape[2]; L = shape[3];
  S0 = stride[0]/sizeof(cdouble);
  S1 = stride[1]/sizeof(cdouble);
  S2 = stride[2]/sizeof(cdouble);
  S3 = stride[3]/sizeof(cdouble);

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
  printf("%d %d\n",nstat,nfft);
  printf("%d %d %d %d\n",M,N,K,L);
  printf("%d %d %d %d\n",S0,S1,S2,S3);
  printf("%g\t",traces[8*S0+7*S1+12*S2+77*S3].real);
  printf("%g\n",traces[8*S0+7*S1+12*S2+77*S3].imag);
  for(m=0;m<ntimes;m++){
    for(n=0;n<nsub;n++){
      printf("ntimes:%d nsub:%d\n",m,n); 
      /****************************************************************************/
      /* first we need the fft'ed window of traces for the stations of this group */
      /****************************************************************************/
	//memcpy((window[j]+1),traces[j][k][l][i],nfft*sizeof(double));

      /*************************************************************/
      /* we start with loop over angular frequency
      /*************************************************************/
      for (w=21;w<22;w++) {
	/***********************************/
	/* now we loop over slowness index */
	/***********************************/
	for (k=0;k<nslow;k++) {
	  vel = 1./slowness[k]*1000.;
	  /************************************/
	  /* now we loop over azimuth index */
	  /************************************/
	  for (l=0;l<nsources;l++) {
	    /********************************************/
	    /* this is the loop over the stations group */
	    /********************************************/
	    sumre = sumim = 0.;
	    for(int q=0;q<nsources;q++){
	      for (j=0;j<M;j++) {
		for(int r=0;r<M;r++){
		  wtau1 = (double) (2.*M_PI*df*(double)w*zetax[l][j]/vel);
		  wtau2 = (double) (2.*M_PI*df*(double)w*zetax[l][j]/vel);
	      re = traces[j*S0+m*S1+n*S2+w*S3].real;
	      im = traces[j*S0+m*S1+n*S2+w*S3].imag;
	      sumre += (double) (re*cos(wtau)-im*sin(wtau));
	      sumim += (double) (im*cos(wtau)+re*sin(wtau));
	    }
	    beam[w][l][k] += (sumre*sumre+sumim*sumim)/nsub;
	  }
	}
      }
    }
    //beam[w][l][k] /= ntimes;
  }
  return 1;
}
