#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Program to read wavfile
int readFile(char *name, short int *signal);

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

short int * input_signal;
int input_size = 0;
short int * ir_signal;
int ir_size = 0;

int main(int argc, char ** argv)
{
	if (argc != 4) {
		printf("Usage: ./Assignment4 [input file] [IR file] [output file]\n");
		return 0;
	}
	input_size = readFile(argv[1], input_signal);
	printf("input_size is %lu\n", input_size);

	ir_size = readFile(argv[2], ir_signal);
	printf("ir_size is %lu\n", ir_size);

	return 0;
}

int readFile(char *name, short int *signal)
{
	FILE *file;
	short int *buffer;
	unsigned long fileLen;
	struct header_file meta;
	//	Open file
	file = fopen(name, "rb");
	if (!file)
	{
		fprintf(stderr, "Unable to open file %s", name);
		exit(-1);
	}
	//	Read header
	fread(&meta, sizeof(meta), 1, file);
	printf(" Size of Header file : %d  bytes\n", sizeof(struct header_file));
	printf(" Sampling rate of the input wave file : %d Hz \n", meta.sample_rate);
	printf(" Bits per sample in wave file : %d \n", meta.bits_per_sample);

	//	Get file length
	fseek(file, 0, SEEK_END);
	fileLen = ftell(file) - sizeof(struct header_file);
	printf("Total bytes for samples: %d\n", fileLen);
	int count = fileLen / 2;
	printf("Number of samples (16bits or 2bytes per sample): %d\n", count);

	fseek(file, sizeof(struct header_file), SEEK_SET);

	//	Allocate memory
	buffer = (short int *)malloc(fileLen + 1);
	if (!buffer)
	{
		fprintf(stderr, "Memory error!");
		fclose(file);
		exit(-1);
	}

	//	Read file contents into buffer
	fread(buffer, fileLen, 1, file);
	fclose(file);
	short int in_signal[count];
	for (int i = 0; i < count; i++)
	{
		in_signal[i] = ((short int*)buffer)[i];
		//printf("buffer[%d] = %d\n", i, input_signal[i]);
	}
	printf("completed creating the input signal\n");
	signal = in_signal;

	free(buffer);
	return count;
}
