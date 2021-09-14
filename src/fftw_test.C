#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <sndfile.h>
#include <fftw3.h>

int main(void)
{
double array[] = {0.1, 0.6, 0.1, 0.4, 0.5, 0, 0.8, 0.7, 0.8, 0.6, 0.1,0};
//double array2[] = {1, 6, 1, 4, 5, 0, 8, 7, 8, 6, 1,0};
double *out;
double *err;
int i,size = 12;


fftw_complex *out_cpx;

fftw_plan fft;
fftw_plan ifft;
out_cpx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*size*10);
out = (double *) malloc(size*10*sizeof(double));
//err = (double *) malloc(size*sizeof(double));

fft = fftw_plan_dft_r2c_1d(size, array, out_cpx, FFTW_ESTIMATE);  //Setup fftw plan for fft
ifft = fftw_plan_dft_c2r_1d(size*10, out_cpx, out, FFTW_ESTIMATE);   //Setup fftw plan for ifft

fftw_execute(fft);
fftw_execute(ifft);

printf("Input:    \tOutput:    \n");
printf("%i\n",sizeof(out));
for(i=0;i<size*10;i++)
{
//err[i] = abs(array[i] - out[i]);    
printf("%f\t%f\n",array[i],out[i]/size);
}

fftw_destroy_plan(fft);
fftw_destroy_plan(ifft);
fftw_free(out_cpx);
free(err);
free(out);
return 0;
}
