#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#define BASE (1LL << 30)
#define FIXED_MAX ((1L << 31)-1)
#define FIXED_MIN (-1L << 31)

#define M_PI (3.14159265358979323846)

#define FILTER_LEN (3)
#define BUFF_LEN (1 << 7)
#define FILTER_MASK (FILTER_LEN - 1)

int32_t FloatToFixed(double x)
{
	if (x >= 2)
		return FIXED_MAX;

	else if (x < -2)
		return FIXED_MIN;

	double base = (double)BASE;
	int32_t res = (int32_t)(x * base);

	return res;
}

void coeffscalc(int32_t *coeffs, float filterfreq, int32_t samplerate, float Q)
{
	double a0, a1, a2, b1, b2, norm;

	double K = tan(M_PI * filterfreq / samplerate);
	double K2 = K * K;

	norm = 1.0 / (1.0 + K / Q + K2);
	a0 = K2 * norm;
	a1 = 2 * a0;
	a2 = a0;
	b1 = 2 * (K2 - 1) * norm;
	b2 = (1 - K / Q + K2) * norm;

	//coeffs[0] = a0;
	//coeffs[1] = a1;
	//coeffs[2] = a2;
	//coeffs[3] = b1;
	//coeffs[4] = b2;

	coeffs[0] = FloatToFixed(a0);
	coeffs[1] = FloatToFixed(a1);
	coeffs[2] = FloatToFixed(a2);
	coeffs[3] = FloatToFixed(b1);
	coeffs[4] = FloatToFixed(b2);
}

int32_t IIR(int32_t* coeffs, int32_t *buffer, int16_t sample, int64_t err)
{
	// 0000 0000 0000 0000 . 000,0 0000 0000 0000 . 0000 0000 0000 0000 . 0000 0000 0000 0000
	int64_t acc = err;
	buffer[0] = buffer[1];
	buffer[1] = buffer[2];
	buffer[2] = (int32_t)sample;
	buffer[3] = buffer[4];
	buffer[4] = buffer[5];


	err = (int64_t)buffer[2] * coeffs[0] + (int64_t)buffer[1] * coeffs[1] + (int64_t)buffer[0] * coeffs[2] - (int64_t)buffer[4] * coeffs[3] - (int64_t)buffer[3] * coeffs[4];
	buffer[5] = (int32_t)(err >> 30);
	err &= (int64_t)0x000000003fffffff;
	return buffer[5];
}

int16_t sweep_signal(int32_t samplerate, float start_freq, float end_freq, float amplitude, int32_t lenght, int t)
{
	int32_t signal;

	float w1 = 2 * M_PI * start_freq;
	float w2 = 2 * M_PI * end_freq;
	float tmp1 = log(w2 / w1);

	float n;
	float tmp2;

	n = (float)t / lenght;
	tmp2 = exp(n * tmp1) - 1.0;
	signal = FloatToFixed(amplitude * sin(w1 * lenght * tmp2 / (samplerate * tmp1)));
	return (int16_t)(signal >> 15);
}

//HEADER
typedef struct
{
	uint8_t         chunkID[4];
	uint32_t        ChunkSize;
	uint8_t         Format[4];

	uint8_t         Subchunk1ID[4];
	uint32_t        Subchunk1Size;
	uint16_t        AudioFormat;
	uint16_t        NumChannels;
	uint32_t        SampleRate;
	uint32_t        ByteRate;
	uint16_t        blockAlign;
	uint16_t        bitsPerSample;

	uint8_t         Subchunk2ID[4];
	uint32_t        Subchunk2Size;
} WavHeader;

void InitHeader(WavHeader *header, int32_t lenght)
{
	header->chunkID[0] = 'R'; header->chunkID[1] = 'I'; header->chunkID[2] = 'F'; header->chunkID[3] = 'F';
	header->Format[0] = 'W'; header->Format[1] = 'A'; header->Format[2] = 'V'; header->Format[3] = 'E';
	header->Subchunk1ID[0] = 'f'; header->Subchunk1ID[1] = 'm'; header->Subchunk1ID[2] = 't'; header->Subchunk1ID[3] = ' ';
	header->Subchunk2ID[0] = 'd'; header->Subchunk2ID[1] = 'a'; header->Subchunk2ID[2] = 't'; header->Subchunk2ID[3] = 'a';
	header->Subchunk1Size = 16;
	header->AudioFormat = 1;
	header->NumChannels = 2;
	header->SampleRate = 48000;
	header->ByteRate = (header->SampleRate * header->SampleRate * header->NumChannels) / 8;
	header->bitsPerSample = 16;
	header->blockAlign = header->NumChannels * header->bitsPerSample / 8;
	header->Subchunk2Size = lenght * header->NumChannels * header->bitsPerSample / 8;
	header->ChunkSize = 4 + (8 + header->Subchunk1Size) + (8 + header->Subchunk2Size);
}

int main(int argc, char* argv[])
{
	float time = 1;
	float freq = 20;
	float end_freq = 20000;
	float amplitude = 0.5;
	int32_t samplerate = 48000;
	float Fc = 200;
	float Q = 0.707;
	int32_t lenght = 0;
	int16_t signal;
	int8_t type;

	printf("start\n");

	lenght = (int32_t)(time * samplerate);

	WavHeader header;

	FILE *file_out;

	InitHeader(&header, lenght);

	file_out = fopen("test_signal.wav", "wb");

	fwrite(&header, sizeof(header), 1, file_out);

	int32_t sample_buffer[6] = { 0, 0, 0, 0, 0, 0 };

	int32_t coeffs[5] = {0, 0, 0, 0, 0};

	coeffscalc(coeffs, Fc, samplerate, Q);

	int64_t acc = 0;

	int16_t buffer[BUFF_LEN * 2];

	int16_t out;

	int n;
	for (int t = 0; t < lenght;)
	{
		for (int j = 0; j < BUFF_LEN; j++, t++)
		{
			
			signal = sweep_signal(samplerate, freq, end_freq, amplitude, lenght, t);

			out = IIR(coeffs, sample_buffer, signal, acc);

			buffer[j * 2] = (int16_t)out;
			buffer[j * 2 + 1] = signal;

			n = j;
		}
		fwrite(buffer, 2 * n * sizeof(int16_t), 1, file_out);
	}

	fclose(file_out);

	puts("\n!Finish!\n");
	system("pause");
	return 0;
}