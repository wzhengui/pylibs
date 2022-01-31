#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define MAX 1000000


void gaussj(double **a, int n, double **b, int m)
{
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv, temp;

  indxc = ivector(1,n);
  indxr = ivector(1,n);
  ipiv = ivector(1,n);
  for (j=1;j<=n;j++) ipiv[j] = 0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j]!=1)
	for (k=1;k<=n;k++){
	  if (ipiv[k]==0){
	    if (fabs(a[j][k]) >= big){
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1"); 
	}     
    ++(ipiv[icol]);
    if (irow != icol){
      for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]);
      for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0) nrerror("gaussj:Singular Matrix-2");
    pivinv = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol){
	dum = a[ll][icol];
	a[ll][icol] = 0.0;
	for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}

char *number2string(int N,char *char_number)
{
  int decimal;

  if (N<= 9){
    char_number[0] = '0'+N;
    char_number[1] = '\0';
  }
  else {
    decimal = N/10;
    char_number[0] = '0'+decimal;
    char_number[1] = '0'+N-decimal*10;
    char_number[2] = '\0';
  }
  return (char_number);
}

double time[MAX], eta[MAX];
main(int argc, char *argv[]){

  struct tide {
      char  label[3];
      double freq;
      double amp;
      double phase;
  };

  FILE  *fid_input, *fid_output, *fid_const;
  char  fname[200];
  double **M, **F;
  struct tide *values;
  double sum1, sum2, sum3, sum4;
  int   N = 0, i = 0, j = 0, l = 0, N_tides;

  void gaussj(double **, int, double **, int);

  if (argc!=4){
    fprintf(stderr,"Usage:%s <input time series> <tidal constituent file> <output file>\n",argv[0]);
    exit(1);
  }

  sprintf(fname,"%s",argv[2]);
  if ((fid_const = fopen(fname,"r")) == NULL){
     fprintf(stderr,"Could not open file %s\n",fname);
     exit(1);
  }
  fscanf(fid_const,"%d",&N_tides);

  values = (struct tide *)calloc(N_tides,sizeof(struct tide));
  M = dmatrix(1,2*N_tides+1,1,2*N_tides+1);
  F = dmatrix(1,2*N_tides+1,1,1);

  for (i=0;i<N_tides;i++){
     fscanf(fid_const,"%s %lf",values[i].label,&(values[i].freq));
  }

  sprintf(fname,"%s",argv[1]);
  if ((fid_input = fopen(fname,"r")) == NULL){
     fprintf(stderr,"Could not open file %s\n",fname);
     exit(1);
  }

  sprintf(fname,"%s",argv[3]);
  if ((fid_output = fopen(fname,"w")) == NULL){
     fprintf(stderr,"Could not open file %s\n",fname);
     exit(1);
  }

  i = 0;
  while(!feof(fid_input)){
     fscanf(fid_input,"%lf%lf\n",&time[i],&eta[i]);
     i++;      
  }
  N = i;

  sum1 = 0.0;
  for (j=0;j<N;j++){
    sum1 += eta[j];
  }
  sum1 /= (double) N;
  fprintf(fid_output,"%s %lf %lf\n","Z0",sum1,0.0);
  for (j=0;j<N;j++){
    eta[j] -= sum1;
  }

  for (j=0;j<N_tides;j++){
     for (l=0;l<N_tides;l++){
       sum1 = 0.0;sum2 = 0.0;sum3 = 0.0;sum4 = 0.0;
       for (i=0;i<N;i++){
	  sum1 += cos(values[l].freq*time[i])*cos(values[j].freq*time[i]);
	  sum2 += sin(values[l].freq*time[i])*cos(values[j].freq*time[i]);
	  sum3 += cos(values[l].freq*time[i])*sin(values[j].freq*time[i]);
	  sum4 += sin(values[l].freq*time[i])*sin(values[j].freq*time[i]);
	}
	M[j+1][l+1] = sum1;
	M[j+1][l+1+N_tides] = sum2;
	M[j+1+N_tides][l+1] = sum3;
	M[j+1+N_tides][l+1+N_tides] = sum4;
     }
     sum1 = 0.0;
     sum2 = 0.0;

     for (i=0;i<N;i++){
       sum1 += eta[i]*cos(values[j].freq*time[i]);
       sum2 += eta[i]*sin(values[j].freq*time[i]);
     }
     F[j+1][1] = sum1;
     F[j+1+N_tides][1] = sum2;
  }

  gaussj(M,2*N_tides,F,1);

  for (j=0;j<N_tides;j++){
     values[j].amp = sqrt(pow(F[j+1][1],2.0)+pow(F[j+1+N_tides][1],2.0));
     values[j].phase = atan2(F[j+1+N_tides][1],F[j+1][1]);
     fprintf(fid_output,"%s %lf %lf\n",values[j].label,values[j].amp,values[j].phase);
  }
  fclose(fid_output);
  
}

