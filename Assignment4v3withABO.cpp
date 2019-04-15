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

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

//function prototypes
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

double *input_signal;
unsigned long input_size = 0;
double *ir_signal;
unsigned long ir_size = 0;
struct header_file headers[3];
unsigned long* size;
double ** signals;
int nextInd = 0;

double * readFile(char *inputFile, header_file &headF) {
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


		double * signal = new double[headF.subchunk2_size / 2];

		cout << "data size :" << (headF.subchunk2_size) << endl;

		int nextFloat = 0;
		for (int i = 40 + headF.subchunk1_size - 16 + 4; i < headF.subchunk2_size; i += 2) {
			short int num = (((input[i]) | 0xff00) & ((input[i + 1] << 8) | 0x00ff));
			if (num > 0) 
				signal[nextFloat] = num / 32767.0;
			else 
				signal[nextFloat] = num / 32768.0;
			nextFloat++;
		}
		return &signal[1];
	}
	else {
		cout << "failed to open file \n exiting..." << endl;
		exit(-1);
		return NULL;
	}
}


//Author:		Leonard Manzara
/* The four1 FFT from Numerical Recipes in C,
  nn must be a power of 2
  isign = +1 for an FFT, and -1 for the Inverse FFT.
  array size must be nn*2. 
  This code assumes the array starts
  at index 1, not 0, so subtract 1 when
  calling the routine (see main() below).
  */

//From class handout

void four1(double data[], int nn, int isign)
{

	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;

	for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j + 1], data[i + 1]);
		}
		m = nn;

		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;

	}
	mmax = 2;
	while (n > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j + 1];
				tempi = wr * data[j + 1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}

void outputToIntArray(double signal[], int *ret, unsigned long int s) {
	for (int i = 0; i < s; i++) {
		ret[i] = (int)(signal[i] * 32768);
	}
}

void writeOutputToFile(unsigned long int numberOfSamples, int out[], char *filename)
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

unsigned long int findPow2() {
	unsigned long int nextPow2 = 1;
	while (nextPow2 < (inputHead.subchunk2_size / 2) || nextPow2 < (irHead.subchunk2_size / 2)) {  // optimize this
		nextPow2 *= 2;
	}
	return nextPow2;
}

void padArrays(unsigned long int nextPow2, double* paddedInputArray, double* paddedIRArray) {
	for (int i = 0; i < nextPow2; i++) {
		if (i < (inputHead.subchunk2_size / 2))
			paddedInputArray[i] = signals[0][i];
		else
			paddedInputArray[i] = 0.0;
	}
	for (int i = 0; i < nextPow2; i++) {
		if (i < (irHead.subchunk2_size / 2))
			paddedIRArray[i] = signals[1][i];
		else
			paddedIRArray[i] = 0.0;
	}
}

void prepareComplexArrays(unsigned long int nextPow2, double paddedInputArray[], double paddedIRArray[], double* complexInputArray, double* complexIRArray) {
	unsigned long int c = 0;
	for (unsigned long i = 0; i < nextPow2 * 2; i += 2) {
		complexInputArray[i] = paddedInputArray[c];
		complexInputArray[i + 1] = 0.0;
		c++;
	}

	c = 0;
	for (unsigned long i = 0; i < nextPow2 * 2; i += 2) {
		complexIRArray[i] = paddedIRArray[c];
		complexIRArray[i + 1] = 0.0;
		c++;
	}
}

int main(int argc, char ** argv)
{
	if (argc != 4) {
		printf("Usage: ./Assignment4 [input file] [IR file] [output file]\n");
		return 0;
	}
	cout << "Begin" << endl;
	signals = new double*[3];
	
	clock_t currentTime;
	currentTime = clock();

	printf("attempting to read %s\n", argv[1]);
	signals[0] = readFile(argv[1], inputHead);

	printf("attempting to read %s\n", argv[2]);
	signals[1] = readFile(argv[2], irHead);

	unsigned long int nextPow2 = findPow2();

	// optimize all places with nextPow2 * 2
	//optimize this with next one 
	double * paddedInputArray = new double[nextPow2];
	double * paddedIRArray = new double[nextPow2];
	padArrays(nextPow2, paddedInputArray, paddedIRArray);
	

	// optimize this with the next one
	double * complexInputArray = new double[nextPow2 * 2];
	double * complexIRArray = new double[nextPow2 * 2];

	four1(complexInputArray - 1, nextPow2, 1);
	four1(complexIRArray - 1, nextPow2, 1);

	double * complexOutArray = new double[nextPow2 * 2];
	for (int i = 0; i < (nextPow2 * 2) - 2; i += 2) {
		complexOutArray[i] = (complexInputArray[i] * complexIRArray[i]) - (complexInputArray[i + 1] * complexIRArray[i + 1]);
		complexOutArray[i + 1] = (complexInputArray[i] * complexIRArray[i + 1]) + (complexInputArray[i + 1] * complexIRArray[i]);
	}

	four1(complexOutArray - 1, nextPow2, -1);

	double max = 0.0;
	for (unsigned long int i = 0; i < nextPow2 * 2; i++) {
		complexOutArray[i] = complexOutArray[i] / ((double)nextPow2);
		if (fabs(complexOutArray[i]) > max) {
			max = fabs(complexOutArray[i]);
		}//*/
	}//*/

	for (unsigned long int i = 0; i < nextPow2; i++) {
		complexOutArray[i] = complexOutArray[2 * i] / max;   // only the normalized real part 
	}

	currentTime = clock() - currentTime;
	printf("Clock cycles taken for convolution: %u\n", currentTime);
	printf("Seconds taken for convolution: %f\n", ((float)currentTime) / CLOCKS_PER_SEC);

	unsigned long int outSize = (inputHead.subchunk2_size / 2) + (irHead.subchunk2_size / 2) - 1;
	int *buffer = new int[outSize];
	outputToIntArray(complexOutArray, buffer, outSize);

	writeOutputToFile(outSize, buffer, argv[3]);

	printf("finished\n");
	return 0;
}



// The following was written by  Dr. Lenord Manzara

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