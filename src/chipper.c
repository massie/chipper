#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "project.h"
#include "chipper.h"
#include "argtable3.h"

void freetable_exit(void **argtable, int exitcode)
{
	arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
	exit(exitcode);
}

const char *get_filename(struct arg_file *arg, const char *default_name)
{
	if (arg->count == 1) {
		return arg->filename[0];
	} else {
		return default_name;
	}
}

const char *get_default_model()
{
	const char *prefix = CHIPPER_INSTALL_PREFIX;
	const char *suffix = "share/chipper/lr.model";
	void *dest = malloc(strlen(prefix) + strlen(suffix) + 2);
	sprintf(dest, "%s/%s", prefix, suffix);
	return dest;
}

int main(int argc, char *argv[])
{
	struct arg_lit *help, *version, *probabilities;
	struct arg_file *input, *output, *model;
	struct arg_dbl *cutoff;
	struct arg_end *end;
	int nerrors;
	char *progname = "chipper";
	const char *default_model = get_default_model();
	void *default_model_help =
	    malloc(strlen(default_model) + strlen("(default:)") + 1);
	sprintf(default_model_help, "(default:%s)", default_model);

	void *argtable[] = {
		help =
		    arg_litn("h", "help", 0, 1, "Display this help and exit"),
		version = arg_litn("v", "version", 0, 1,
				   "Display version info and exit"),
		input = arg_filen("i", "input", "<fasta file>", 0, 1,
				  "FASTA file with protein(s) (default: stdin)"),
		model = arg_filen("m", "model", "<liblinear model>", 0, 1,
				  default_model_help),
		output = arg_filen("o", "output", "<fastq file>", 0, 1,
				   "FASTQ file with predictions (default: stdout)"),
		cutoff =
		    arg_dbln("c", "cutoff", "<cutoff>", 0, 1,
			     "Cutoff in range (0,1) for probability models (default: cutoff with highest MCC)"),
		probabilities = arg_litn("p", "probs", 0, 1,
					 "Output probabilities to FASTQ"),
		end = arg_end(20),
	};

	nerrors = arg_parse(argc, argv, argtable);

	if (help->count > 0) {
		printf("Usage: %s", progname);
		arg_print_syntax(stdout, argtable, "\n");
		arg_print_glossary(stdout, argtable, " %-25s %s\n");
		freetable_exit(argtable, EXIT_SUCCESS);
	}

	if (version->count > 0) {
		printf("chipper version: %d.%d.%d\n",
		       CHIPPER_VERSION_MAJOR, CHIPPER_VERSION_MINOR,
		       CHIPPER_VERSION_PATCH);
		printf("Best cutoff is %.2f (MCC= %.2f)\n",
		       BEST_CUTOFF_VALUE, BEST_CUTOFF_MCC);
		freetable_exit(argtable, EXIT_SUCCESS);
	}

	if (nerrors > 0) {
		/* Display the error details contained in the arg_end struct. */
		arg_print_errors(stdout, end, progname);
		printf("Try '%s --help' for more information.\n", progname);
		freetable_exit(argtable, EXIT_FAILURE);
	}

	return predict_cleavage(get_filename(input, NULL),
				get_filename(output, NULL),
				get_filename(model, default_model),
				probabilities->count,
				cutoff->count == 1, cutoff->dval[0]);
}
