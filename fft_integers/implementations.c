#include "implementations.h"
#include "twiddle.h"
#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>


void bit_reverse_copy(int16_t (*Input)[2], int N, int32_t (*X)[2])
{
    int nBits = log2_int(N);

    for (int i = 0; i < N; i++)
    {
        int j = 0;
        int x = i;
        for (int b = 0; b < nBits; b++)
        {
            j = (j << 1) | (x & 1);
            x >>= 1;
        }

        X[j][0] = Input[i][0];
        X[j][1] = Input[i][1];
    }
}

void fft_fixed_q15(int16_t (*Input)[2], int N, int32_t (*X)[2])
{
    int32_t Q = 15;

    bit_reverse_copy(Input, N, X);

    int32_t stage = 1;
    while(stage < N)
    {
        int32_t halfStage = stage;
        stage *= 2;
        
        int32_t step = N / stage;
        for(int32_t start = 0; start < N; start+=stage)
        {
            for(int32_t k = 0; k < halfStage; k++)
            {
                int16_t wr = twiddle[k*step][0];
                int16_t wi = twiddle[k*step][1];

                int32_t a_real = X[start+k][0];
                int32_t a_imag = X[start+k][1];

                int32_t b_real = X[start+k+halfStage][0];
                int32_t b_imag = X[start+k+halfStage][1];

                // t = W * b
                int32_t t_real = (b_real * (int32_t)wr - b_imag * (int32_t)wi) >> Q;
                int32_t t_imag = (b_real * (int32_t)wi + b_imag * (int32_t)wr) >> Q;

                // Butterfly
                X[start + k][0] = a_real + t_real;
                X[start + k][1] = a_imag + t_imag;

                X[start + k + halfStage][0] = a_real - t_real;
                X[start + k + halfStage][1] = a_imag - t_imag;
            }
        }

        // Scale to prevent overflow
        for(int32_t i = 0; i < N; i++)
        {
            X[i][0] >>= 1;
            X[i][1] >>= 1;
        }
    }
}

void fft_fixed_q15_inplace(int32_t (*X)[2], int N)
{
    int32_t Q = 15;

    // 1. In-place Bit-Reversal
    // This shuffles X in-place so we don't need a secondary Input buffer.
    for (int i = 0; i < N; i++) 
    {
        int j = reverse_bits(i, N);
        if (i < j)
        {
            // Swap X[i] and X[j]
            int32_t temp_r = X[i][0];
            int32_t temp_i = X[i][1];
            X[i][0] = X[j][0];
            X[i][1] = X[j][1];
            X[j][0] = temp_r;
            X[j][1] = temp_i;
        }
    }

    // 2. FFT Stages
    int32_t stage = 1;
    while(stage < N)
    {
        int32_t halfStage = stage;
        stage *= 2;
        int32_t step = N / stage;

        for(int32_t start = 0; start < N; start += stage)
        {
            for(int32_t k = 0; k < halfStage; k++)
            {
                int16_t wr = twiddle[k*step][0];
                int16_t wi = twiddle[k*step][1];

                int32_t a_real = X[start+k][0];
                int32_t a_imag = X[start+k][1];

                int32_t b_real = X[start+k+halfStage][0];
                int32_t b_imag = X[start+k+halfStage][1];

                // t = W * b (using your verified 32-bit cast logic)
                int32_t t_real = (b_real * (int32_t)wr - b_imag * (int32_t)wi) >> Q;
                int32_t t_imag = (b_real * (int32_t)wi + b_imag * (int32_t)wr) >> Q;

                // Butterfly
                X[start + k][0] = a_real + t_real;
                X[start + k][1] = a_imag + t_imag;

                X[start + k + halfStage][0] = a_real - t_real;
                X[start + k + halfStage][1] = a_imag - t_imag;
            }
        }

        // Scale to prevent overflow
        for(int32_t i = 0; i < N; i++)
        {
            X[i][0] >>= 1;
            X[i][1] >>= 1;
        }
    }
}

void fft_fixed_q15_inplace_no_scaling(int32_t (*X)[2], int N)
{
    int32_t Q = 15;

    // In-place Bit-Reversal
    for (int i = 0; i < N; i++) {
        int j = reverse_bits(i, N);
        if (i < j) {
            // Swap X[i] and X[j]
            int32_t temp_r = X[i][0];
            int32_t temp_i = X[i][1];
            X[i][0] = X[j][0];
            X[i][1] = X[j][1];
            X[j][0] = temp_r;
            X[j][1] = temp_i;
        }
    }

    int32_t stage = 1;
    while(stage < N)
    {
        int32_t halfStage = stage;
        stage *= 2;
        int32_t step = N / stage;

        for(int32_t start = 0; start < N; start += stage)
        {
            for(int32_t k = 0; k < halfStage; k++)
            {
                int16_t wr = twiddle[k*step][0];
                int16_t wi = twiddle[k*step][1];

                int32_t a_real = X[start+k][0];
                int32_t a_imag = X[start+k][1];

                int32_t b_real = X[start+k+halfStage][0];
                int32_t b_imag = X[start+k+halfStage][1];

                // t = W * b
                int32_t t_real = (b_real * (int32_t)wr - b_imag * (int32_t)wi) >> Q;
                int32_t t_imag = (b_real * (int32_t)wi + b_imag * (int32_t)wr) >> Q;

                // Butterfly
                X[start + k][0] = a_real + t_real;
                X[start + k][1] = a_imag + t_imag;

                X[start + k + halfStage][0] = a_real - t_real;
                X[start + k + halfStage][1] = a_imag - t_imag;
            }
        }
    }
}

void ifft_fixed_q15_inplace(int32_t (*X)[2], int N) 
{
    // 1. Conjugate
    for(int i = 0; i < N; i++) X[i][1] = -X[i][1];

    // 2. In-place FFT
    fft_fixed_q15_inplace_no_scaling(X, N);

    // 3. Conjugate again
    for(int i = 0; i < N; i++) X[i][1] = -X[i][1];
}

void dump_complex_signal_16(int16_t (*A)[2], int Size, FILE* fp_real, FILE* fp_imag)
{
  fprintf(fp_real, "[");
  fprintf(fp_imag, "[");
  for (int i = 0; i < Size; i++) 
  {
    fprintf(fp_real, "%d,", A[i][0]);
    fprintf(fp_imag, "%d,", A[i][1]);
  }
  fprintf(fp_real, "]");
  fprintf(fp_imag, "]");
}


void dump_complex_signal_32(int32_t (*A)[2], int Size, FILE* fp_real, FILE* fp_imag)
{
  fprintf(fp_real, "[");
  fprintf(fp_imag, "[");
  for (int i = 0; i < Size; i++) 
  {
    fprintf(fp_real, "%d,", A[i][0]);
    fprintf(fp_imag, "%d,", A[i][1]);
  }
  fprintf(fp_real, "]");
  fprintf(fp_imag, "]");
}

void dump_long_int_array(long int* A, int Size, FILE* fp)
{
  fprintf(fp, "[");
  for(int i = 0; i < Size; i++)
    fprintf(fp, "%ld,", A[i]);
  fprintf(fp, "]");
}
void dump_int_array(int* A, int Size, FILE* fp)
{
  fprintf(fp, "[");
  for(int i = 0; i < Size; i++)
    fprintf(fp, "%d,", A[i]);
  fprintf(fp, "]");
}
