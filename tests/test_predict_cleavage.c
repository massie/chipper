#include "chipper.h"
#include "minunit.h"
#include "util.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int tests_run = 0;

char prediction_to_char(double prediction);

static char *test_probability_to_char() {
  mu_assert("prediction_to_char(0.0) failed", prediction_to_char(0.0) == '0');
  mu_assert("prediction_to_char(0.1) failed", prediction_to_char(0.1) == '1');
  mu_assert("prediction_to_char(0.11) failed", prediction_to_char(0.11) == '1');
  mu_assert("prediction_to_char(0.2) failed", prediction_to_char(0.2) == '2');
  mu_assert("prediction_to_char(0.99) failed", prediction_to_char(0.99) == '9');
  mu_assert("prediction_to_char(1.0) failed", prediction_to_char(1.0) == ':');
  return 0;
}

/* See predict_cleavage.c */
void init_features(struct feature_node *features);
void peptide_to_features(const char *seq, int seq_len,
                         struct feature_node *features, int features_len);

static char *test_pos_neg_samples() {
  const char *positive_sample = "SDAIGPREQWEPVFQNGKMA";
  const char *negative_sample = "LVGHSDAIGPREQWEPVFQN";
  struct feature_node
      positive_features[GENERATED_SAMPLE_LEN * NUM_AA_PROPERTIES + 1];
  struct feature_node
      negative_features[GENERATED_SAMPLE_LEN * NUM_AA_PROPERTIES + 1];
  struct model *data_model = load_model("../models/lr.model");
  double probabilities_true_false[2];
  double positive_probability, negative_probability;

  if (data_model == NULL) {
    fprintf(stderr, "Unable to load model\n");
    exit(EXIT_FAILURE);
  }

  init_features(positive_features);
  init_features(negative_features);
  peptide_to_features(positive_sample, GENERATED_SAMPLE_LEN, positive_features,
                      ARRAY_LEN(positive_features));
  peptide_to_features(negative_sample, GENERATED_SAMPLE_LEN, negative_features,
                      ARRAY_LEN(negative_features));

  predict_probability(data_model, positive_features, probabilities_true_false);
  positive_probability = probabilities_true_false[0];
  predict_probability(data_model, negative_features, probabilities_true_false);
  negative_probability = probabilities_true_false[0];

  printf("pos=%f neg=%f\n", positive_probability, negative_probability);
  assert(positive_probability > negative_probability);
  return 0;
}

static char *all_tests() {
  mu_run_test(test_probability_to_char);
  mu_run_test(test_pos_neg_samples);
  return 0;
}

int main(int argc, char **argv) {
  char *result = all_tests();
  if (result != 0) {
    printf("%s\n", result);
  } else {
    printf("ALL TESTS PASSED\n");
  }
  printf("Tests run: %d\n", tests_run);

  return result != 0;
}
