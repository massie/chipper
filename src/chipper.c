#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "version.h"
#include "chipper.h"

void print_usage(char *program)
{
	fprintf(stderr,
		"Usage %s -i <fasta input> -o <fastq output> -m <model>\n",
		program);
}

int main(int argc, char *argv[])
{
	extern char *optarg;
	extern int optind;
	int opt;
	char *input_fasta_file = "-";
	char *output_fastq_file = "-";
	char *model_file = "chipper.model";

	while ((opt = getopt(argc, argv, "i:o:m:hv")) != -1) {
		switch (opt) {
		case 'i':
			input_fasta_file = optarg;
			break;
		case 'o':
			output_fastq_file = optarg;
			break;
		case 'm':
			model_file = optarg;
			break;
		case 'v':
			printf("chipper version: %d.%d.%d\n",
			       CHIPPER_VERSION_MAJOR, CHIPPER_VERSION_MINOR,
			       CHIPPER_VERSION_PATCH);
			printf("Best cutoff is %.2f (MCC= %.2f)\n",
			       BEST_CUTOFF_VALUE, BEST_CUTOFF_MCC);
			exit(EXIT_SUCCESS);
		case 'h':
			print_usage(argv[0]);
			exit(EXIT_SUCCESS);
		default:
			print_usage(argv[0]);
			exit(EXIT_FAILURE);
		}
	}

	return predict_cleavage(input_fasta_file, output_fastq_file,
				model_file, yes_no);
}
