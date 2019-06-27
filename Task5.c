#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#define BASE (1LL << 31)
#define FIXED_MAX ((1L << 31)-1)
#define FIXED_MIN (-1L << 31)

#define M_PI (3.14159265358979323846)

#define FILTER_LEN (1 << 7)
#define BUFF_LEN (1 << 7)
#define FILTER_MASK (FILTER_LEN - 1)

double coeffscalc(float filterfreq, int32_t samplerate, float Q, float gain)
{
	double a0, a1, a2, b1, b2, norm;

	double V = pow(10, abs(gain) / 20);
	double K = tan(M_PI * filterfreq / samplerate);
	double K2 = K * K;

	norm = 1 / (1 + K / Q + K2);
	a0 = K2 * norm;
	a1 = 2 * a0;
	a2 = a0;
	b1 = 2 * (K2 - 1) * norm;
	b2 = (1 - K / Q + K2) * norm;
}