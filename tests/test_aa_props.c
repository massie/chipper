#include "chipper.h"
#include "minunit.h"
#include <stdio.h>

int tests_run = 0;

struct feature_node test_features[NUM_AA_PROPERTIES + 1];

static char *test_non_alpha() {
  mu_assert("didn't return zero on non-alphabet character",
            get_aa_properties('1', test_features) == 0);
  return 0;
}

static char *test_alpha_non_aa() {
  mu_assert("didn't return zero on alpha character that is not an amino acid",
            get_aa_properties('B', test_features) == 0);
  return 0;
}

static char *test_valid_aa() {
  mu_assert("didn't return number of properties with valid amino acid",
            get_aa_properties('A', test_features) == NUM_AA_PROPERTIES);
  return 0;
}

static char *test_lowercase() {
  mu_assert("lowercase amino acid notation failed",
            get_aa_properties('a', test_features) ==
                get_aa_properties('A', test_features));
  return 0;
}

static char *all_tests() {
  mu_run_test(test_non_alpha);
  mu_run_test(test_alpha_non_aa);
  mu_run_test(test_valid_aa);
  mu_run_test(test_lowercase);
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
