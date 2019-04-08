#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <endian.h>


using namespace std;


struct header_file
{
	char chunk_id[4];
	int chunk_size;
	char format[4];
	char subchunk1_id[4];
	int subchunk1_size;
	short int audio_format;
	short int num_channels;
	int sample_rate;			// sample_rate denotes the sampling rate.
	int byte_rate;
	short int block_align;
	short int bits_per_sample;
	char subchunk2_id[4];
	long int subchunk2_size;			// subchunk2_size denotes the number of samples.
};

//Chunks
struct chunk_t
{
	char ID[4]; //"data" = 0x61746164
	unsigned long size;  //Chunk data bytes
};

float *input_signal;
unsigned long input_size = 0;
float *ir_signal;
unsigned long ir_size = 0;
struct header_file headers[3];
unsigned long size[3];
float ** signals;
int nextInd = 0;

int readFile(char *name)
{
	FILE *file;
	short int *buffer;
	unsigned long fileLen;
	//	Open file
	file = fopen(name, "rb");
	if (!file)
	{
		fprintf(stderr, "Unable to open file %s", name);
		exit(-1);
	}
	//	Read header
	fread(&headers[nextInd], sizeof(struct header_file), 1, file);
	printf(" Size of Header file : %d  bytes\n", sizeof(struct header_file));
	printf(" Sampling rate of the input wave file : %d Hz \n", headers[nextInd].sample_rate);
	printf(" Bits per sample in wave file : %d \n", headers[nextInd].bits_per_sample);

	//	Get file length
	fseek(file, 0, SEEK_END);
	fileLen = ftell(file);
	printf("Total bytes for samples: %d\n", fileLen - sizeof(struct header_file));
	int count = (fileLen - sizeof(struct header_file)) / 2;
	printf("Number of samples (16bits or 2bytes per sample): %d\n", count);

	fseek(file, sizeof(struct header_file), SEEK_SET);

	short int *value = new short int[count];
	//memset(value, 0, sizeof(short int) * samples_count);
	printf("value initialized\n");

	signals[nextInd] = new float[count];
	//memset(value, 0, sizeof(float) * samples_count);

	printf("signal initialized\n");

	//Reading data
	for (int i = 0; i < count; i++)
	{
		fread(&value[i], 2, 1, file);
		signals[nextInd][i] = value[i] / 32768.0;
	}
	printf("signals[200] = %f\n", signals[nextInd][200]);

	size[nextInd] = count;

	nextInd++;

	return count;
}

/*****************************************************************************
*
*    Function:     convolve
*
*    Description:  Convolves two signals, producing an output signal.
*                  The convolution is done in the time domain using the
*                  "Input Side Algorithm" (see Smith, p. 112-115).
*
*    Parameters:   x[] is the signal to be convolved
*                  N is the number of samples in the vector x[]
*                  h[] is the impulse response, which is convolved with x[]
*                  M is the number of samples in the vector h[]
*                  y[] is the output signal, the result of the convolution
*                  P is the number of samples in the vector y[].  P must
*                       equal N + M - 1
*
*****************************************************************************/

void convolve(float x[], int N, float h[], int M, float y[], int P)
{
	int n, m;

	/*  Make sure the output buffer is the right size: P = N + M - 1  */
	if (P != (N + M - 1)) {
		printf("Output signal vector is the wrong size\n");
		printf("It is %-d, but should be %-d\n", P, (N + M - 1));
		printf("Aborting convolution\n");
		return;
	}

	/*  Clear the output buffer y[] to all zero values  */
	for (n = 0; n < P; n++)
		y[n] = 0.0;

	int counter = 0;
	/*  Do the convolution  */
	/*  Outer loop:  process each input value x[n] in turn  */
	for (n = 0; n < N; n++) {
		if ((counter = (counter + 1) % 1000) == 0)
			printf("working on loop n = %d out of %d\n", n, N);
		/*  Inner loop:  process x[n] with each sample of h[]  */
		for (m = 0; m < M; m++)
			y[n + m] += x[n] * h[m];
	}
}

inline void endian_swap(unsigned int& x)
{
	x = (x >> 24) |
		((x << 8) & 0x00FF0000) |
		((x >> 8) & 0x0000FF00) |
		(x << 24);
}

int main(int argc, char ** argv)
{
	if (argc != 4) {
		printf("Usage: ./Assignment4 [input file] [IR file] [output file]\n");
		return 0;
	}
	signals = new float*[3];

	input_size = readFile(argv[1]);
	//printf("input_size is %lu\n", size[0]);
	bool useImpulse = true;
	if (!useImpulse) {
		ir_size = readFile(argv[2]);
	} else {
		signals[1] = new float[1];
		signals[1][0] = 1.0;
		size[1] = 1;//*/
	}
	//printf("ir_size is %lu\n", size[1]);
	
	//printf("input_size is %lu\n", input_size);
	//printf("input_signal[200] is %f\n", signals[0][200]);
	//printf("ir_signal[200] is %f\n", signals[1][200]);
	size[2] = size[0] + size[1] - 1;
	signals[2] = new float[size[2]];

	convolve(signals[0], size[0], signals[1], size[1], signals[2], size[2]);
	printf("finished the convolution\n");
	short int *buffer = new short int[sizeof(struct header_file) / 2 + size[2]];
	printf("initialized the ouput buffer\n");

	//file = fopen(argv[1], "rb");
	//fread(&buffer, sizeof(struct header_file), 1, file);
	//fclose(file);
	printf("changing the size of the header\n");
	headers[0].subchunk2_size = htole32(size[2]);
	memcpy(buffer, &headers[0], sizeof(struct header_file));

	printf("converting the signal back into short ints\n");
	for (unsigned long i = 0; i < size[2]; i++) {
		buffer[sizeof(struct header_file) + i] = (short int)(signals[2][i] * 32768);
	}
	//memcpy(buffer + 40, &size[2], 4);
	
	FILE *file;
	file = fopen(argv[3], "w");
	fwrite(&buffer, sizeof(struct header_file) + size[2], 1, file);
	fclose(file);

	printf("finished\n");
	return 0;
}
