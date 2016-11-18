#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <zlib.h>
#include <math.h>
#include <assert.h>
#include "project.h"
#include "chipper.h"
#include "linear.h"
#include "output.h"
#include "kseq.h"
#include "util.h"

KSEQ_INIT(gzFile, gzread)

void init_features(struct feature_node *features)
{
	int i;
	for (i = 0; i < NUM_FEATURES; i++) {
		/* note: liblinear uses 1-based indexing */
		features[i].index = i + 1;
	}
	/* Mark the last position */
	features[NUM_FEATURES].index = -1;
}

void
peptide_to_features(const char *peptide, int peptide_len,
		    struct feature_node *features, int features_len)
{
	int i, rval;
	assert(features_len == peptide_len * NUM_AA_PROPERTIES + 1);
	for (i = 0; i < peptide_len; i++) {
		rval =
		    get_aa_properties(peptide[i],
				      features + i * NUM_AA_PROPERTIES);
		if (rval == 0) {
			fprintf(stderr, "Illegal amino acid='%c'", peptide[i]);
			exit(EXIT_FAILURE);
		}
	}
}

struct chipper_output *new_chipper_output(const char *filename,
					  chipper_output_format format,
					  int is_probability_model,
					  int output_probabilities,
					  int cutoff_provided, double cutoff)
{
	double active_cutoff = cutoff_provided ? cutoff : BEST_CUTOFF_VALUE;

	if (!is_probability_model && output_probabilities) {
		fprintf(stderr,
			"Model does not support outputting probabilities. Exiting.\n");
		exit(EXIT_FAILURE);
	}

	switch (format) {
	case FASTQ_OUTPUT:
		return new_fastq_output(filename, output_probabilities,
					active_cutoff);
	case NETCHOP_OUTPUT:
		return new_netchop_output(filename, active_cutoff, 0);
	case NETCHOP_SHORT_OUTPUT:
		return new_netchop_output(filename, active_cutoff, 1);
	default:
		fprintf(stderr, "Unknown format requested. Exiting.\n");
		exit(EXIT_FAILURE);
	}
}

int
predict_cleavage(const char *fasta_filename, const char *output_filename,
		 chipper_output_format output_format, const char *model_file,
		 int output_probabilities, int cutoff_provided, double cutoff)
{
	gzFile input;
	kseq_t *seq;
	int i, len, sample_index;
	struct model *data_model;
	/* The feature vector has a terminator node '+1' */
	struct feature_node features[NUM_FEATURES + 1];
	double *predictions;
	int predictions_len = 0;
	double probabilities_true_false[2];
	int is_probability_model;
	struct chipper_output *output = NULL;
	int cleavage_site;

	predictions_len = 256;
	predictions = malloc(predictions_len * sizeof(double));
	if (predictions == NULL) {
		fprintf(stderr,
			"Unable to allocate memory for predictions. Exiting.\n");
		exit(EXIT_FAILURE);
	}

	/* Initialize feature vector */
	init_features(features);

	/* Load the model */
	data_model = load_model(model_file);
	if (data_model == NULL) {
		fprintf(stderr, "Unable to load model file '%s'\n", model_file);
		exit(EXIT_FAILURE);
	}

	input = gzFile_create(fasta_filename, "r");
	seq = kseq_init(input);
	is_probability_model = check_probability_model(data_model);
	output =
	    new_chipper_output(output_filename, output_format,
			       is_probability_model, output_probabilities,
			       cutoff_provided, cutoff);

	while ((len = kseq_read(seq)) >= 0) {
		if (strchr(seq->seq.s, 'U') || strchr(seq->seq.s, 'u')) {
			fprintf(stderr,
				"Ignoring protein with Selenocysteine. Not supported (yet).\n");
			continue;
		}
		/* Make sure we have enough space for predictions */
		if (predictions_len < seq->seq.l) {
			predictions =
			    realloc(predictions, seq->seq.l * sizeof(double));
			if (predictions == NULL) {
				fprintf(stderr,
					"Unable to allocate memory for prediction data. Exiting.\n");
				exit(EXIT_FAILURE);
			}
			predictions_len = seq->seq.l;
		}
		/* Initialize all position to no information */
		for (i = 0; i < seq->seq.l; i++) {
			predictions[i] = NAN;
		}

		for (sample_index = GENERATED_SAMPLE_LEN;
		     sample_index <= seq->seq.l; sample_index++) {
			peptide_to_features(seq->seq.s + sample_index -
					    GENERATED_SAMPLE_LEN,
					    GENERATED_SAMPLE_LEN, features,
					    ARRAY_LEN(features));
			cleavage_site =
			    sample_index - GENERATED_SAMPLE_LEN / 2 - 1;
			if (is_probability_model) {
				predict_probability(data_model, features,
						    probabilities_true_false);
				predictions[cleavage_site] =
				    probabilities_true_false[0];
			} else {
				predictions[cleavage_site] =
				    predict(data_model, features);
			}
		}

		output->write_record_cb(output, seq->name.s, seq->name.l,
					seq->comment.s, seq->comment.l,
					seq->seq.s, seq->seq.l, predictions);
	}
	output->close_cb(output);
	kseq_destroy(seq);
	gzclose(input);
	return 0;
}
