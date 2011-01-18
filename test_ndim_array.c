#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
typedef struct {double real; double imag;} cdouble;

int test_2dim_array(double **arr){
  int i,j,k;
  for(i=0;i<2;i++){
    for(j=0;j<3;j++){
      printf("%g\n",arr[i][j]);
    }
  }
  arr[0][0] = 111.11;
  return 1;
}


int test_4dim_array(double ****arr){
  int i,j,k,l;
  printf("%g\n",arr[1][3][2][3]);
  for(i=0;i<2;i++){
    for(j=0;j<4;j++){
      for(k=0;k<4;k++){
	for(l=0;l<6;l++){
	  printf("%g\n",arr[i][j][k][l]);
	}
      }
    }
  }
  arr[1][1][2][3] = 222.22;
  return 1;
}


int test_complex_array(cdouble *arr, long *strides, long *shape){
  
  int M,N,K,L,S0,S1,S2,S3,i,j,k,l;
  M = shape[0]; N = shape[1]; K = shape[2]; L = shape[3];
  S0 = strides[0]/sizeof(cdouble);
  S1 = strides[1]/sizeof(cdouble);
  S2 = strides[2]/sizeof(cdouble);
  S3 = strides[3]/sizeof(cdouble);
  /*for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      for(k=0;k<K;k++){
	for(l=0;l<L;l++){
	  printf("%g\t",arr[i*S0+j*S1+k*S2+l*S3].real);
	  printf("%g\n",arr[i*S0+j*S1+k*S2+l*S3].imag);
	}
      }
    }
    }*/
  printf("%g\t",arr[2*S0+1*S1+2*S2+0*S3].real);
  printf("%g\n",arr[2*S0+1*S1+2*S2+0*S3].imag);
  return 1;
}


