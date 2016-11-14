#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <zlib.h>
#include <math.h>
#include "kseq.h"
#include "version.h"
#include "chipper.h"
#include "linear.h"

#define CHECK(x) do { \
  int retval = (x);\
  if (retval <= 0) {\
    fprintf(stderr, "Runtime error: %s returned %d at %s:%d", #x, retval, __FILE__, __LINE__);\
	exit(EXIT_FAILURE);\
  }\
} while (0)

KSEQ_INIT(gzFile, gzread)
gzFile gzFile_create(const char *filename, char *mode)
{
	FILE *stream = NULL;
	if (filename == NULL) {
		if (strcmp(mode, "r") == 0) {
			stream = stdin;
		} else if (strcmp(mode, "w") == 0) {
			stream = stdout;
		} else {
			fprintf(stderr, "Unknown file mode='%s'\n", mode);
			exit(EXIT_FAILURE);
		}
	} else {
		stream = fopen(filename, mode);
		if (!stream) {
			fprintf(stderr, "file:%s, error: %s", filename,
				strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
	return gzdopen(fileno(stream), mode);
}

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

char prediction_to_char(double prediction)
{
	return (char)((int)floor(prediction * 10) + 48);
}

int
predict_cleavage(const char *fasta_input, const char *fastq_output,
		 const char *model_file, int output_probabilities)
{
	const int sample_half_len = GENERATED_SAMPLE_LEN / 2;
	gzFile input, output;
	kseq_t *seq;
	int i, len, cut_index;
	char *preds = NULL;
	int preds_len = 0;
	struct model *data_model;
	/* The feature vector has a terminator node '+1' */
	struct feature_node features[NUM_FEATURES + 1];
	double prediction;
	double probabilities_true_false[2];
	char amino_acid;
	int rval;

	/* Initialize feature vector */
	init_features(features);

	/* Load the model */
	data_model = load_model(model_file);
	if (data_model == NULL) {
		fprintf(stderr, "Unable to load model file '%s'\n", model_file);
		exit(EXIT_FAILURE);
	}

	input = gzFile_create(fasta_input, "r");
	seq = kseq_init(input);
	output = gzFile_create(fastq_output, "w");

	while ((len = kseq_read(seq)) >= 0) {
		if (seq->qual.l) {
			fprintf(stderr,
				"Input file is a FASTQ file. Exiting.\nqual: %s\n",
				seq->qual.s);
			exit(EXIT_FAILURE);
		}
		if (strchr(seq->seq.s, 'U') || strchr(seq->seq.s, 'u')) {
			fprintf(stderr,
				"Ignoring protein with Selenocysteine. Not supported (yet).\n");
			continue;
		}
		/* Echo back the name and comment */
		CHECK(gzwrite(output, "@", 1));
		CHECK(gzwrite(output, seq->name.s, seq->name.l));
		if (seq->comment.l) {
			CHECK(gzwrite(output, " ", 1));
			CHECK(gzwrite(output, seq->comment.s, seq->comment.l));
		}
		CHECK(gzwrite(output, "\n", 1));

		/* Make sure we have enough space for predictions */
		if (preds_len < seq->seq.l) {
			preds = realloc(preds, seq->seq.l);
			if (preds == NULL) {
				fprintf(stderr,
					"Unable to allocate memory for prediction data. Exiting.\n");
				exit(EXIT_FAILURE);
			}
			preds_len = seq->seq.l;
		}
		/* Initialize all position to no information */
		memset(preds, '!', seq->seq.l);

		for (cut_index = sample_half_len;
		     cut_index < seq->seq.l - sample_half_len; cut_index++) {
			for (i = 0; i < GENERATED_SAMPLE_LEN; i++) {
				amino_acid =
				    seq->seq.s[cut_index + i - sample_half_len];
				rval =
				    get_aa_properties(amino_acid,
						      features +
						      i * NUM_AA_PROPERTIES);
				if (rval == 0) {
					fprintf(stderr,
						"Illegal amino acid='%c'.\n",
						amino_acid);
					exit(EXIT_FAILURE);
				}
			}

			if (output_probabilities) {
				prediction =
				    predict_probability(data_model, features,
							probabilities_true_false);
				preds[cut_index] =
				    prediction_to_char(probabilities_true_false
						       [0]);
			} else {
				prediction = predict(data_model, features);
				preds[cut_index] =
				    prediction_to_char(prediction);
			}
		}

		CHECK(gzwrite(output, seq->seq.s, seq->seq.l));
		CHECK(gzwrite(output, "\n+", 2));
		CHECK(gzwrite(output, seq->name.s, seq->name.l));
		if (seq->comment.l) {
			CHECK(gzwrite(output, " ", 1));
			CHECK(gzwrite(output, seq->comment.s, seq->comment.l));
		}
		CHECK(gzwrite(output, "\n", 1));
		CHECK(gzwrite(output, preds, seq->seq.l));
		CHECK(gzwrite(output, "\n", 1));
	}
	kseq_destroy(seq);
	gzclose(output);
	gzclose(input);
	return 0;
}
