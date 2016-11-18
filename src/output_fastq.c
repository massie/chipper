#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include "output.h"
#include "util.h"

struct chipper_fastq_output {
	struct chipper_output output;
	gzFile file;
	int output_probabilities;
};

/* Visible for testing */
char prediction_to_char(double prediction)
{
	if (isnan(prediction)) {
		return '!';
	} else {
		return (char)((int)floor(prediction * 10) + 48);
	}
}

static void fastq_write_record(struct chipper_output *co, char *name,
			       int name_len, char *comment, int comment_len,
			       char *seq, int seq_len, double *predictions)
{
	struct chipper_fastq_output *fq =
	    container_of(co, struct chipper_fastq_output, output);
	int i;

	/* Echo back the name and comment */
	CHECK(gzwrite(fq->file, "@", 1));
	CHECK(gzwrite(fq->file, name, name_len));
	if (comment_len) {
		CHECK(gzwrite(fq->file, " ", 1));
		CHECK(gzwrite(fq->file, comment, comment_len));
	}
	CHECK(gzwrite(fq->file, "\n", 1));

	CHECK(gzwrite(fq->file, seq, seq_len));
	CHECK(gzwrite(fq->file, "\n+", 2));
	CHECK(gzwrite(fq->file, name, name_len));
	if (comment_len) {
		CHECK(gzwrite(fq->file, " ", 1));
		CHECK(gzwrite(fq->file, comment, comment_len));
	}
	CHECK(gzwrite(fq->file, "\n", 1));
	for (i = 0; i < seq_len; i++) {
		if (fq->output_probabilities) {
			CHECK(gzprintf
			      (fq->file, "%c",
			       prediction_to_char(predictions[i])));
		} else {
			CHECK(gzprintf
			      (fq->file, "%c",
			       isnan(predictions[i]) ? '!' :
			       predictions[i] >= co->cutoff ? '+' : '-'));
		}
	}
	CHECK(gzwrite(fq->file, "\n", 1));
}

static void fastq_close(struct chipper_output *co)
{
	struct chipper_fastq_output *fq =
	    container_of(co, struct chipper_fastq_output, output);
	gzclose(fq->file);
}

struct chipper_output *new_fastq_output(const char *output_path,
					int output_probabilities, double cutoff)
{
	struct chipper_fastq_output *fq =
	    malloc(sizeof(struct chipper_fastq_output));
	if (fq == NULL) {
		fprintf(stderr, "Can't allocate memory for fastq output\n");
		exit(EXIT_FAILURE);
	}
	fq->file = gzFile_create(output_path, "w");
	fq->output_probabilities = output_probabilities;

	struct chipper_output *out = &fq->output;
	out->format = FASTQ_OUTPUT;
	out->write_record_cb = fastq_write_record;
	out->close_cb = fastq_close;
	out->cutoff = cutoff;
	return out;
}
