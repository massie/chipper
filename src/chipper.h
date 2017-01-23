
/*** THIS IS AUTO-GENERATED SOURCE! DO NOT EDIT! ***/

#ifndef CHIPPER_H_DEFINED
#define CHIPPER_H_DEFINED

#include "linear.h"
#include "output.h"

#define NUM_AA_PROPERTIES 50
#define TRAINING_SET_SIZE 52332
#define BEST_CUTOFF_VALUE 0.36
#define BEST_CUTOFF_MCC 0.620
#define GENERATED_SAMPLE_LEN 20
#define NUM_FEATURES (GENERATED_SAMPLE_LEN * NUM_AA_PROPERTIES)

int get_aa_properties(char aa, struct feature_node *features);
int predict_cleavage(const char *fasta_input, const char *output,
		     chipper_output_format format, const char *model_file,
		     int output_probabilities, int cutoff_provided,
		     double cutoff);

#endif
