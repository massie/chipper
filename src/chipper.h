
/*** THIS IS AUTO-GENERATED SOURCE! DO NOT EDIT! ***/

#ifndef CHIPPER_H_DEFINED
#define CHIPPER_H_DEFINED

#include "linear.h"

#define NUM_AA_PROPERTIES 50
#define TRAINING_SET_SIZE 49159
#define BEST_CUTOFF_VALUE 0.36
#define BEST_CUTOFF_MCC 0.604
#define GENERATED_SAMPLE_LEN 28
#define NUM_FEATURES (GENERATED_SAMPLE_LEN * NUM_AA_PROPERTIES)

typedef enum output_mode_t { yes_no, probability } output_mode;

int get_aa_properties(char aa, struct feature_node *features);
int predict_cleavage(char *fasta_input, char *fastq_output, char *model_file,
		     output_mode mode);

#endif
