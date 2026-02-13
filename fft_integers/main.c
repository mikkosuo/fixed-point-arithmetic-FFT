#include <stdio.h>
#include "twiddle.h"
#include "implementations.h"
#include "math.h"
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include "test_data.h"


int main()
{
  int32_t N = 1024;
  printf("Started fft integers test.\n");

  FILE* fp_real = fopen("output_main/output_real.txt","w");
  FILE* fp_imag = fopen("output_main/output_imag.txt","w");
  dump_complex_signal_16(signal, N, fp_real, fp_imag);
  printf("Mallocing \n");
  int32_t (*fft_output)[2] = calloc((int)(N), sizeof(*fft_output));
  printf("Computing fft \n");
  //fft_fixed_q15(signal, N, fft_output);
  int32_t (*signal32)[2] = malloc(N*sizeof(*signal32));
  for(int i = 0; i < N; i++)
  {
    signal32[i][0] = (int32_t)signal[i][0];
    signal32[i][1] = (int32_t)signal[i][1];
  }
  fft_fixed_q15_inplace(signal32, N);
  //ifft_fixed_q15_inplace(signal32, N);
  printf("Printing sinusoidal output\n");

  FILE* fp_fft_real = fopen("output_main/output_fft_real.txt","w");
  FILE* fp_fft_imag = fopen("output_main/output_fft_imag.txt","w");
  dump_complex_signal_32(signal32, N, fp_fft_real, fp_fft_imag);
  free(fft_output);
  free(signal32);
  return 0;
}
