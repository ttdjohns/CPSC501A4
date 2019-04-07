#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

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
		exit(-1);
	}

	/*  Clear the output buffer y[] to all zero values  */
	for (n = 0; n < P; n++)
		y[n] = 0.0;

	/*  Do the convolution  */
	/*  Outer loop:  process each input value x[n] in turn  */
	for (n = 0; n < N; n++) {
		/*  Inner loop:  process x[n] with each sample of h[]  */
		for (m = 0; m < M; m++)
			y[n + m] += x[n] * h[m];
	}
}

int main(int argc, char ** argv)
{
	/// parse command line arguments
	if (argc != 4 && argc != 3) {
		printf("Usage: ./Assignment4 [input file] [IR file] [output file]\n");
		return 0;
	}

	FILE* inputFile = fopen(argv[1], "rb");
	unsigned long inputLen;
	struct header_file inputHeader;
	if (inputFile == NULL) {
		printf("Could not open input file\n");
		exit(-1);
	}

	//	Read header
	fread(&inputHeader, sizeof(), 1, inputFile);
	//printf(" Size of inputHeader file : %d  bytes\n", sizeof(struct header_file));
	//printf(" Sampling rate of the input wave file : %d Hz \n", inputHeader.sample_rate);
	//printf(" Bits per sample in wave file : %d \n", inputHeader.bits_per_sample);

	//	Get file length
	fseek(inputFile, 0, SEEK_END);
	inputLen = ftell(inputFile);
	//printf("Total bytes for samples: %d\n", inputLen);
	int countInput = inputLen / 2;
	//printf("Number of samples (16bits or 2bytes per sample): %d\n", countInput);

	fseek(inputFile, 0, SEEK_SET);

	short int *inputBuffer;
	//	Allocate memory
	inputBuffer = (short int *)malloc(inputLen + 1);
	if (!inputBuffer)
	{
		fprintf(stderr, "Memory error!");
		fclose(inputFile);
		return;
	}

	//	Read file contents into buffer
	fread(inputBuffer, inputLen, 1, inputFile);
	fclose(inputFile);
	/*short int input_signal[inputCount];
	for (int i = 0; i <= countInput; i++)
	{
		input_signal[i] = ((short int*)inputBuffer)[i];
		printf("inputBuffer[%d] = %d\n", i, input_signal[i]);
	}*/


	FILE* irFile = fopen(argv[2], "rb");
	unsigned long irLen;
	struct header_file irHeader;
	if (irFile == NULL) {
		printf("Could not open IR file\n");
		exit(-1);
	}

	//	Read header
	fread(&irHeader, sizeof(), 1, irFile);
	//printf(" Size of inputHeader file : %d  bytes\n", sizeof(struct header_file));
	//printf(" Sampling rate of the input wave file : %d Hz \n", irHeader.sample_rate);
	//printf(" Bits per sample in wave file : %d \n", irHeader.bits_per_sample);

	//	Get file length
	fseek(irFile, 0, SEEK_END);
	irLen = ftell(irFile);
	//printf("Total bytes for samples: %d\n", irLen);
	int countIR = irLen / 2;
	//printf("Number of samples (16bits or 2bytes per sample): %d\n", countIR);

	fseek(irFile, 0, SEEK_SET);

	FILE* outputFile = fopen(argv[3], "w");
	if (outputFile == NULL) {
		printf("Could not open output file\n");
		exit(-1);
	}

	short int *irBuffer;
	//	Allocate memory
	irBuffer = (short int *)malloc(irLen + 1);
	if (!irBuffer)
	{
		fprintf(stderr, "Memory error!");
		fclose(irFile);
		return;
	}

	//	Read file contents into buffer
	fread(irBuffer, irLen, 1, irFile);
	fclose(irFile);
	//short int ir_signal[inputCount];
	/*for (int i = 0; i <= countIR; i++)
	{
		ir_signal[i] = ((short int*)irBuffer)[i];
		printf("irBuffer[%d] = %d\n", i, ir_signal[i]);
	}*/

	short int *outputBuffer;
	//	Allocate memory
	int outputLen = (inputLen - sizeof(struct header_file)) + (irLen - sizeof(struct header_file)) - 1;
	outputBuffer = (short int *)malloc(sizeof(struct header_file) + outputLen + 1);
	if (!outputBuffer)
	{
		fprintf(stderr, "Memory error!");
		fclose(outputFile);
		return;
	}

	for (int i = 0; i < sizeof(struct header_file) + outputLen + 1; i++) {
		outputBuffer[i] = 0;
	}

	memcpy(outputBuffer, inputBuffer, sizeof(struct header_file));

	double *input_signal;
	//	Allocate memory
	int inSigLen = inputLen - sizeof(struct header_file);
	input_signal = (short int *)malloc(inSigLen + 1);
	if (!input_signal)
	{
		fprintf(stderr, "Memory error!");
		exit(-1);
	}

	for (int i = 0; i < inSigLen; i++) {
		input_signal[i] = ((double)inputBuffer[i + sizeof(struct header_file)]) / ((double)32768);
	}

	double *ir_signal;
	//	Allocate memory
	int irSLen = irLen - sizeof(struct header_file);
	ir_signal = (short int *)malloc(irSLen + 1);
	if (!ir_signal)
	{
		fprintf(stderr, "Memory error!");
		exit(-1);
	}

	for (int i = 0; i < inSigLen; i++) {
		ir_signal[i] = ((double)inputBuffer[i + sizeof(struct header_file)]) / ((double)32768);
	}

	double *output_signal;
	//	Allocate memory
	int outSigLen = inSigLen + irSLen - 1
	output_signal = (short int *)malloc(outSigLen + 1);
	if (!output_signal)
	{
		fprintf(stderr, "Memory error!");
		exit(-1);
	}

	convolve(input_signal, inSigLen, ir_signal, irSLen,
		output_signal, outSigLen);


	free(inputBuffer);
	free(irBuffer);
	free(outputBuffer);
	fclose(outputFile)
	return 0;
}