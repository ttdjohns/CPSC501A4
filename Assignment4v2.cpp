#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <endian.h>
#include <math.h>


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
	fileLen = ftell(file) - sizeof(struct header_file);
	printf("Total bytes for samples: %d\n", fileLen - sizeof(struct header_file));
	int count = fileLen / 2;
	printf("Number of samples (16bits or 2bytes per sample): %d\n", count);

	fseek(file, sizeof(struct header_file), SEEK_SET);

	char *value = new char[fileLen];
	fread(&value, fileLen, 1, file);
	signals[nextInd] = new float[count];
	int nextShort = 0;
	for (int i = 0; i < fileLen; i += 2) {
		signals[nextInd][nextFloat++] = ((value[i] | 0xff00) & ((value[i + 1] << 8) | 0x00f)) / 32768.0;
	}

	//memset(value, 0, sizeof(float) * samples_count);

	printf("signal initialized\n");

	size[nextInd] = count;

	nextInd++;

	return count;
}

//function prototypes
void convolve(float* x, int N, float* h, int M, float y[], int P);
//void print_vector(char *title, float x[], int N);
void createTestTone(double frequency, double duration,
	int numberOfChannels, char *filename);
void writeWaveFileHeader(int channels, int numberSamples,
	double outputRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);

int main(int argc, char ** argv)
{
	if (argc != 4) {
		printf("Usage: ./Assignment4 [input file] [IR file] [output file]\n");
		return 0;
	}
	signals = new float*[3];

	input_size = readFile(argv[1]);
	//printf("input_size is %lu\n", size[0]);
	bool useImpulse = false;
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
	
	//creating output array
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



// The following was written by  Dr. Lenord Manzara


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

	/*  Do the convolution  */
	/*  Outer loop:  process each input value x[n] in turn  */
	for (n = 0; n < N; n++) {
		/*  Inner loop:  process x[n] with each sample of h[]  */
		for (m = 0; m < M; m++)
			y[n + m] += x[n] * h[m];
	}
}



/*****************************************************************************
*
*    Function:     print_vector
*
*    Description:  Prints the vector out to the screen
*
*    Parameters:   title is a string naming the vector
*                  x[] is the vector to be printed out
*                  N is the number of samples in the vector x[]
*
*****************************************************************************/

void print_vector(char *title, float x[], int N)
{
	int i;

	printf("\n%s\n", title);
	printf("Vector size:  %-d\n", N);
	printf("Sample Number \tSample Value\n");
	for (i = 0; i < N; i++)
		printf("%-d\t\t%f\n", i, x[i]);
}
/******************************************************************************
*
*       function:       createTestTone
*
*       purpose:        Calculates and writes out a sine test tone to file
*
*       arguments:      frequency:  frequency of the test tone in Hz
*                       duration:  length of the test tone in seconds
*                       numberOfChannels:  number of audio channels
*                       filename:  name of the file to create
*
*       internal
*       functions:      writeWaveFileHeader, fwriteShortLSB
*
*       library
*       functions:      ceil, pow, fopen, fprintf, sin, rint, fclose
*
******************************************************************************/

void createTestTone(double frequency, double duration,
	int numberOfChannels, char *filename)
{
	int i;

	/*  Calculate the number of sound samples to create,
		rounding upwards if necessary  */
	int numberOfSamples = (int)ceil(duration * SAMPLE_RATE);

	/*  Calculate the maximum value of a sample  */
	int maximumValue = (int)pow(2.0, (double)BITS_PER_SAMPLE - 1) - 1;

	/*  Open a binary output file stream for writing */
	FILE *outputFileStream = fopen(filename, "wb");
	if (outputFileStream == NULL) {
		fprintf(stderr, "File %s cannot be opened for writing\n", filename);
		return;
	}

	/*  Write the WAVE file header  */
	writeWaveFileHeader(numberOfChannels, numberOfSamples,
		SAMPLE_RATE, outputFileStream);

	/*  Create the sine tone and write it to file  */
	/*  Since the frequency is fixed, the angular frequency
		and increment can be precalculated  */
	double angularFrequency = 2.0 * PI * frequency;
	double increment = angularFrequency / SAMPLE_RATE;
	for (i = 0; i < numberOfSamples; i++) {
		/*  Calculate the sine wave in the range -1.0 to + 1.0  */
		double value = sin(i * increment);

		/*  Convert the value to a 16-bit integer, with the
			range -maximumValue to + maximumValue.  The calculated
			value is rounded to the nearest integer  */
		short int sampleValue = rint(value * maximumValue);

		/*  Write out the sample as a 16-bit (short) integer
			in little-endian format  */
		fwriteShortLSB(sampleValue, outputFileStream);

		/*  If stereo output, duplicate the sample in the right channel  */
		if (numberOfChannels == STEREOPHONIC)
			fwriteShortLSB(sampleValue, outputFileStream);
	}

	/*  Close the output file stream  */
	fclose(outputFileStream);
}



/******************************************************************************
*
*       function:       writeWaveFileHeader
*
*       purpose:        Writes the header in WAVE format to the output file.
*
*       arguments:      channels:  the number of sound output channels
*                       numberSamples:  the number of sound samples
*                       outputRate:  the sample rate
*                       outputFile:  the output file stream to write to
*
*       internal
*       functions:      fwriteIntLSB, fwriteShortLSB
*
*       library
*       functions:      ceil, fputs
*
******************************************************************************/

void writeWaveFileHeader(int channels, int numberSamples,
	double outputRate, FILE *outputFile)
{
	/*  Calculate the total number of bytes for the data chunk  */
	int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;

	/*  Calculate the total number of bytes for the form size  */
	int formSize = 36 + dataChunkSize;

	/*  Calculate the total number of bytes per frame  */
	short int frameSize = channels * BYTES_PER_SAMPLE;

	/*  Calculate the byte rate  */
	int bytesPerSecond = (int)ceil(outputRate * frameSize);

	/*  Write header to file  */
	/*  Form container identifier  */
	fputs("RIFF", outputFile);

	/*  Form size  */
	fwriteIntLSB(formSize, outputFile);

	/*  Form container type  */
	fputs("WAVE", outputFile);

	/*  Format chunk identifier (Note: space after 't' needed)  */
	fputs("fmt ", outputFile);

	/*  Format chunk size (fixed at 16 bytes)  */
	fwriteIntLSB(16, outputFile);

	/*  Compression code:  1 = PCM  */
	fwriteShortLSB(1, outputFile);

	/*  Number of channels  */
	fwriteShortLSB((short)channels, outputFile);

	/*  Output Sample Rate  */
	fwriteIntLSB((int)outputRate, outputFile);

	/*  Bytes per second  */
	fwriteIntLSB(bytesPerSecond, outputFile);

	/*  Block alignment (frame size)  */
	fwriteShortLSB(frameSize, outputFile);

	/*  Bits per sample  */
	fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

	/*  Sound Data chunk identifier  */
	fputs("data", outputFile);

	/*  Chunk size  */
	fwriteIntLSB(dataChunkSize, outputFile);
}



/******************************************************************************
*
*       function:       fwriteIntLSB
*
*       purpose:        Writes a 4-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteIntLSB(int data, FILE *stream)
{
	unsigned char array[4];

	array[3] = (unsigned char)((data >> 24) & 0xFF);
	array[2] = (unsigned char)((data >> 16) & 0xFF);
	array[1] = (unsigned char)((data >> 8) & 0xFF);
	array[0] = (unsigned char)(data & 0xFF);
	return fwrite(array, sizeof(unsigned char), 4, stream);
}



/******************************************************************************
*
*       function:       fwriteShortLSB
*
*       purpose:        Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
	unsigned char array[2];

	array[1] = (unsigned char)((data >> 8) & 0xFF);
	array[0] = (unsigned char)(data & 0xFF);
	return fwrite(array, sizeof(unsigned char), 2, stream);
}