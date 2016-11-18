#ifndef OUTPUT_H_DEFINED
#define OUTPUT_H_DEFINED

enum chipper_output_format {
	FASTQ_OUTPUT,
	NETCHOP_OUTPUT,
	NETCHOP_SHORT_OUTPUT
};
typedef enum chipper_output_format chipper_output_format;

struct chipper_output {
	chipper_output_format format;
	double cutoff;
	void (*write_record_cb) (struct chipper_output * output, char *name,
				 int name_len, char *comment, int comment_len,
				 char *seq, int seq_len, double *predictions);
	void (*close_cb) (struct chipper_output * output);
};

struct chipper_output *new_fastq_output(const char *output_path,
					int output_probabilities,
					double cutoff);

struct chipper_output *new_netchop_output(const char *output_path,
					  double cutoff, int short_output);

#endif
