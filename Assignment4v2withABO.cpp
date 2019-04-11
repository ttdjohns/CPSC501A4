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

double *input_signal;
unsigned long input_size = 0;
double *ir_signal;
unsigned long ir_size = 0;
struct header_file headers[3];
unsigned long size[3];
double ** signals;
double ** fftSignals;
int nextInd = 0;

#define PI 3.14159265358979323846264338327950288

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

	signals[nextInd] = new double[count];
	//memset(value, 0, sizeof(float) * samples_count);

	printf("signal initialized\n");

	//Reading data
	for (int i = 0; i < count; i++)
	{
		fread(&value[i], 2, 1, file);
		signals[nextInd][i] = value[i] / 32768.0;
	}
	//printf("signals[200] = %f\n", signals[nextInd][200]);


	size[nextInd] = count;

	nextInd++;

	return count;
}


inline void endian_swap(unsigned int& x)
{
	x = (x >> 24) |
		((x << 8) & 0x00FF0000) |
		((x >> 8) & 0x0000FF00) |
		(x << 24);
}

void swap(double a, double b) {
	double temp = a;
	a = b;
	b = a;
}

void fft(double* data, unsigned long nn)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	// reverse-binary reindexing
	n = nn << 1;
	j = 1;
	printf("entering first for in fft\n");
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			swap(data[j - 1], data[i - 1]);
			swap(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	};
	printf("finished forst for in fft\n");
	// here begins the Danielson-Lanczos section
	mmax = 2;
	while (n > mmax) {
		printf("in fft while with n = %lu and mmax = %lu\n", n, mmax);
		istep = mmax << 1;
		theta = -(2 * M_PI / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j - 1] - wi * data[j];
				tempi = wr * data[j] + wi * data[j - 1];

				data[j - 1] = data[i - 1] - tempr;
				data[j] = data[i] - tempi;
				data[i - 1] += tempr;
				data[i] += tempi;
			}
			wtemp = wr;
			wr += wr * wpr - wi * wpi;
			wi += wi * wpr + wtemp * wpi;
		}
		mmax = istep;
	}
}

void
ifft(double **v, unsigned long n, double **tmp)
{
	if (n > 1) {			/* otherwise, do nothing and return */
		unsigned long k, m;
		double *z, *w, **vo, **ve;
		z = new double[2];
		w = new double[2];
		ve = tmp; vo = tmp + n / 2;
		for (k = 0; k < n / 2; k++) {
			ve[k] = v[2 * k];
			vo[k] = v[2 * k + 1];
		}
		ifft(ve, n << 1, v);		/* FFT on even-indexed elements of v[] */
		ifft(vo, n << 1, v);		/* FFT on odd-indexed elements of v[] */
		for (m = 0; m < n / 2; m++) {
			w[0] = cos(2 * m * PI / (double)n);
			w[1] = sin(2 * m * PI / (double)n);
			z[0] = w[0] * vo[m][0] - w[1] * vo[m][1];	/* Re(w*vo[m]) */
			z[1] = w[0] * vo[m][1] + w[1] * vo[m][0];	/* Im(w*vo[m]) */
			v[m][0] = ve[m][0] + z[0];
			v[m][1] = ve[m][1] + z[1];
			v[m + n / 2][0] = ve[m][0] - z[0];
			v[m + n / 2][1] = ve[m][1] - z[1];
		}
	}
	return;
}

void writeHeader(FILE* file) {
	unsigned short int *siBuff = new unsigned short int[4];
	unsigned int *iBuff = new unsigned int[4];

	fwrite(headers[0].chunk_id, 4, 1, file);
	iBuff[0] = htole32(headers[0].chunk_size);
	fwrite(iBuff, 4, 1, file);
	fwrite(headers[0].format, 4, 1, file);
	fwrite(headers[0].subchunk1_id, 4, 1, file);
	iBuff[0] = htole32(headers[0].subchunk1_size);
	fwrite(iBuff, 4, 1, file);
	siBuff[0] = htole16(headers[0].audio_format);
	fwrite(siBuff, 2, 1, file);
	siBuff[0] = htole16(headers[0].num_channels);
	fwrite(siBuff, 2, 1, file);
	iBuff[0] = htole32(headers[0].sample_rate);
	fwrite(iBuff, 4, 1, file);
	iBuff[0] = htole32(headers[0].byte_rate);
	fwrite(iBuff, 4, 1, file);
	siBuff[0] = htole16(headers[0].block_align);
	fwrite(siBuff, 2, 1, file);
	siBuff[0] = htole16(headers[0].bits_per_sample);
	fwrite(siBuff, 2, 1, file);
	fwrite(headers[0].subchunk2_id, 4, 1, file);
	iBuff[0] = htole32(headers[0].subchunk2_size);
	fwrite(iBuff, 4, 1, file);

	return;
}

int main(int argc, char ** argv)
{
	if (argc != 4) {
		printf("Usage: ./Assignment4 [input file] [IR file] [output file]\n");
		return 0;
	}
	signals = new double*[3];

	input_size = readFile(argv[1]);
	//printf("input_size is %lu\n", size[0]);
	bool useImpulse = true;
	if (!useImpulse) {
		ir_size = readFile(argv[2]);
	} else {
		signals[1] = new double[1];
		signals[1][0] = 1.0;
		size[1] = 1;//*/
	}

	/*size[2] = size[0] + size[1] - 1;
	signals[2] = new float[size[2]];*/

	unsigned long nn = 1;
	while (nn < 2 * size[0]) {
		nn >> 1;
	}


	fftSignals = new double*[2];
	//convolve(signals[0], size[0], signals[1], size[1], signals[2], size[2]);
	for (int i = 0; i < 2; i++) {
		fftSignals[i] = new double[nn];
	}

	// put the real part in every second index
	for (int i = 0; i < nn; i++) {
		for (int k = 0; k < 2; k++) {					// optimize this part
			if ((i % 2) == 0) {							// optimize this part too 
				fftSignals[k][i] = signals[k][i << 1];
			}
			else {
				fftSignals[k][i] = 0;
			}
		}
	}
	printf("applying fft\n");
	nn = nn << 1;
	fft(fftSignals[0], nn);
	fft(fftSignals[1], nn);
	printf("fft completed \n");
	nn = nn >> 1;
	for (unsigned long i = 0; i < nn; i++) {
		fftSignals[0][i] += fftSignals[1][i];
	}

	nn = nn << 1;
	double ** v = new double*[nn];
	double ** tmp = new double*[nn];
	for (int i = 0; i < nn; i++) {
		v[i] = new double[2];
		tmp[i] = new double[2];
		tmp[i][0] = 0;
		tmp[i][1] = 0;
		v[i][0] = fftSignals[0][(2 * i)];
		v[i][1] = fftSignals[0][(2 * i) + 1];
	}

	printf("applying ifft \n");
	ifft(v, nn, tmp);

	printf("finished the ifft\n");
	short int *buffer = new short int[sizeof(struct header_file) / 2 + size[0]];
	printf("initialized the ouput buffer\n");

	//file = fopen(argv[1], "rb");
	//fread(&buffer, sizeof(struct header_file), 1, file);
	//fclose(file);
	printf("changing the size of the header\n");
	headers[0].subchunk2_size = size[0];
	//memcpy(buffer, &headers[0], sizeof(struct header_file));


	printf("converting the signal back into short ints\n");
	for (unsigned long i = 0; i < size[0]; i++) {
		buffer[i] = htole16((short int)(v[i][0] * 32768));
	}
	//memcpy(buffer + 40, &size[2], 4);
	
	FILE *file;
	file = fopen(argv[3], "w");
	writeHeader(file);
	fwrite(&buffer, 2, input_size, file);
	fclose(file);

	printf("finished\n");
	return 0;
}
