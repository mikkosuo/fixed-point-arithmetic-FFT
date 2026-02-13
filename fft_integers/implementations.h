#pragma once
#include <inttypes.h>
#include <stdio.h>

static inline int32_t log2_int(int32_t N)
{
    int32_t n = 0;
    while( N >>= 1)
        n++;
    return n;
}

static inline int32_t reverse_bits(int32_t i, int32_t N)
{
    int32_t reversed = 0;
    int32_t n = i;
    
    // We only need to reverse the bits within the range of N
    // If N=8, we reverse 3 bits (log2(8))
    for (int32_t bit = 1; bit < N; bit <<= 1)
    {
        reversed <<= 1;        // Shift result left
        reversed |= (n & 1);   // Add the LSB of n to the result
        n >>= 1;               // Shift n right to get the next bit
    }
    return reversed;
}

// A faster version that avoids the loop if you know your bit-depth
static inline int32_t reverse_bits_fast(uint32_t x, int log2N) {
    x = ((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1);
    x = ((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2);
    x = ((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4);
    x = ((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8);
    x = (x >> 16) | (x << 16);
    return x >> (32 - log2N);
}

void bit_reverse_copy(int16_t (*Input)[2], int32_t N, int32_t (*X)[2]);
void fft_fixed_q15(int16_t (*A)[2], int N_size, int32_t (*C)[2]);
void fft_fixed_q15_inplace(int32_t (*X)[2], int N);
void fft_fixed_q15_inplace_no_scaling(int32_t (*X)[2], int N);
void ifft_fixed_q15_inplace(int32_t (*X)[2], int N);

void dump_complex_signal_16(int16_t (*A)[2], int Size, FILE* fp_real, FILE* fp_imag);
void dump_complex_signal_32(int32_t (*A)[2], int Size, FILE* fp_real, FILE* fp_imag);
void dump_long_int_array(long int* A, int Size, FILE* Fp);
void dump_int_array(int* A, int Size, FILE* Fp);