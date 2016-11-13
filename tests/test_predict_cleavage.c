#include <stdio.h>
#include "minunit.h"

int tests_run = 0;

char prediction_to_char(double prediction);

static char *test_probability_to_char()
{
	mu_assert("prediction_to_char(0.0) failed",
		  prediction_to_char(0.0) == '0');
	mu_assert("prediction_to_char(0.1) failed",
		  prediction_to_char(0.1) == '1');
	mu_assert("prediction_to_char(0.11) failed",
		  prediction_to_char(0.11) == '1');
	mu_assert("prediction_to_char(0.2) failed",
		  prediction_to_char(0.2) == '2');
	mu_assert("prediction_to_char(0.99) failed",
		  prediction_to_char(0.99) == '9');
	mu_assert("prediction_to_char(1.0) failed",
		  prediction_to_char(1.0) == ':');
	return 0;
}

static char *all_tests()
{
	mu_run_test(test_probability_to_char);
	return 0;
}

int main(int argc, char **argv)
{
	char *result = all_tests();
	if (result != 0) {
		printf("%s\n", result);
	} else {
		printf("ALL TESTS PASSED\n");
	}
	printf("Tests run: %d\n", tests_run);

	return result != 0;
}
