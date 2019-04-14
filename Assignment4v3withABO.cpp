#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <endian.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <time.h> 


using namespace std;


/*  CONSTANTS  ***************************************************************/
#define PI                3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY         440.0

/*  Test tone duration in seconds  */
#define DURATION          2.0				

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2

//function prototypes
void convolve(float* x, unsigned long int N, float* h, unsigned long int M, float y[], unsigned long int P);
//void print_vector(char *title, float x[], int N);
void writeWaveFileHeader(int channels, int numberSamples,
	double outputRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);

struct header_file
{
	char chunk_id;
	unsigned int chunk_size;
	char * format;
	char * subchunk1_id;
	unsigned int subchunk1_size;
	unsigned short int audio_format;
	unsigned short int num_channels;
	unsigned int sample_rate;			// sample_rate denotes the sampling rate.
	unsigned int byte_rate;
	unsigned short int block_align;
	unsigned short int bits_per_sample;
	char * subchunk2_id;
	unsigned long int subchunk2_size;			// subchunk2_size denotes the number of samples.
} inputHead, irHead;

float *input_signal;
unsigned long input_size = 0;
float *ir_signal;
unsigned long ir_size = 0;
struct header_file headers[3];
unsigned long* size;
float ** signals;
int nextInd = 0;

float * readFile(char *inputFile, header_file &headF) {
	cout << "Entered readFile" << endl;
	ifstream file(inputFile, ios::in | ios::binary | ios::ate);

	if (file.is_open()) {
		streampos fileSize = file.tellg();
		file.seekg(0, ios::beg);
		char * input;
		input = new char[fileSize];
		file.read(input, fileSize);
		file.close();

		//Copying from input to provided header file 
		//cout << "memcpy1" << endl;
		memcpy(&headF.chunk_id, &input[0], 4);
		//cout << "memcpy2" << endl;
		memcpy(&headF.chunk_size, &input[4], 4);
		//cout << "memcpy3" << endl;
		memcpy(&headF.format, &input[8], 4);
		//cout << "memcpy4" << endl;
		memcpy(&headF.subchunk1_id, &input[12], 4);
		//cout << "memcpy5" << endl;
		memcpy(&headF.subchunk1_size, &input[16], 4);
		//cout << "memcpy6" << endl;
		memcpy(&headF.audio_format, &input[18], 2);
		//cout << "memcpy7" << endl;
		memcpy(&headF.num_channels, &input[20], 2);
		//cout << "memcpy8" << endl;
		memcpy(&headF.sample_rate, &input[22], 4);

		//cout << "memcpy9" << endl;

		memcpy(&headF.subchunk2_size, &input[40 + headF.subchunk1_size - 16], 4);
		/*headF.subchunk2_size = (headF.subchunk2_size >> 24) |
								((headF.subchunk2_size << 8) & 0x00FF0000) |
								((headF.subchunk2_size >> 8) & 0x0000FF00) |
								(headF.subchunk2_size << 24);//*/


		float * signal = new float[headF.subchunk2_size / 2];

		cout << "data size :" << (headF.subchunk2_size) << endl;

		int nextFloat = 0;
		for (int i = 40 + headF.subchunk1_size - 16 + 4; i < headF.subchunk2_size; i += 2) {
			signal[nextFloat] = (((input[i]) | 0xff00) & ((input[i + 1] << 8) | 0x00ff)) / 32768.0;
			nextFloat++;
		}
		return &signal[0];
	}
	else {
		cout << "failed to open file \n exiting..." << endl;
		exit(-1);
		return NULL;
	}
}

void outputToIntArray(int *ret, unsigned long int s) {
	for (int i = 0; i < s; i++) {
		ret[i] = (int)(signals[2][i] * 32768);
	}
}

void writeOutputToFile(int numberOfSamples, int out[], char *filename)
{
	int i;

	/*  Open a binary output file stream for writing */
	FILE *outputFileStream = fopen(filename, "wb");
	if (outputFileStream == NULL) {
		fprintf(stderr, "File %s cannot be opened for writing\n", filename);
		return;
	}

	/*  Write the WAVE file header  */
	writeWaveFileHeader(MONOPHONIC, numberOfSamples,
		SAMPLE_RATE, outputFileStream);

	for (i = 0; i < numberOfSamples; i++) {

		fwriteShortLSB(out[i], outputFileStream);
	}

	/*  Close the output file stream  */
	fclose(outputFileStream);
}

int main(int argc, char ** argv)
{
	if (argc != 4) {
		printf("Usage: ./Assignment4 [input file] [IR file] [output file]\n");
		return 0;
	}
	cout << "Begin" << endl;
	signals = new float*[3];

	printf("attempting to read %s\n", argv[1]);
	signals[0] = readFile(argv[1], inputHead);

	printf("attempting to read %s\n", argv[2]);
	signals[1] = readFile(argv[2], irHead);

	//creating output array
	cout << "Initilizing output arrays" << endl;
	unsigned long int outSize = (inputHead.subchunk2_size / 2) + (irHead.subchunk2_size / 2) - 1;
	signals[2] = new float[outSize + 1];
	cout << "finished initilizing output arrays" << endl;

	//printf("begining convolution with inputHead.subchunk2_size = %lu, irHead.subchunk2_size = %lu, outSize = %lu\n", inputHead.subchunk2_size, irHead.subchunk2_size, outSize);
	clock_t currentTime;
	currentTime = clock();

	convolve(signals[0], (inputHead.subchunk2_size / 2), signals[1], (irHead.subchunk2_size / 2), signals[2], outSize);

	currentTime = clock() - currentTime;
	printf("Clock cycles taken for convolution: %d\n", currentTime);
	printf("Seconds taken for convolution: %f\n", ((float)currentTime) / CLOCKS_PER_SEC);

	/*
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
	*/
	int *buffer = new int[outSize];
	outputToIntArray(buffer, outSize);

	writeOutputToFile((int)outSize, buffer, argv[3]);

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

void convolve(float x[], unsigned long int N, float h[], unsigned long int M, float y[], unsigned long int P)
{
	unsigned long int n, m;

	/*  Make sure the output buffer is the right size: P = N + M - 1  */
	if (P != (N + M - 1)) {
		printf("Output signal vector is the wrong size\n");
		printf("It is %-d, but should be %-d\n", P, (N + M - 1));
		printf("Aborting convolution\n");
		return;
	}
	//cout << "clearing output buffer" << endl;
	/*  Clear the output buffer y[] to all zero values  */
	for (n = 0; n < P; n++)
		y[n] = 0.0;

	//cout << "doing convolution" << endl;
	/*  Do the convolution  */
	/*  Outer loop:  process each input value x[n] in turn  */
	int num = 0;
	for (n = 0; n < N; n++) {
		/*if (num == 0) {
			printf("working through n = %lu \n", n);
		}
		num = (num + 1) % 100000;//*/
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