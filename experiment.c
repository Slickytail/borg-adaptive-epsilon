/**
 * An experiment that runs the BORG MOEA on a COCO suite.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/wait.h>

#include "coco.h"
#include "borg.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

/**
 * The number of threads to use for the experiment.
 * Problems in a suite will be divided among threads and each thread will run solve problems synchronously.
 */
static const int THREADS = 4;

/**
 * The maximal budget for evaluations done by an optimization algorithm equals dimension * BUDGET_MULTIPLIER.
 * Increase the budget multiplier value gradually to see how it affects the runtime.
 */
static const unsigned int BUDGET_MULTIPLIER = 1e3;

/**
 * The number of evaluations to spend on estimating the epsilon values.
 */
static const unsigned int EPSILON_EVALUATIONS = BUDGET_MULTIPLIER / 10;

/**
 * Scale multiplier for epsilon values.
 */
static const float EPSILON_SCALE = 1e-5;

/**
 * The random seed used for the experiment.
 */
static const unsigned int RANDOM_SEED = 0xc0c0;

/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);

/**
 * A pointer to the problem to be optimized (needed in order to simplify the interface between the optimization
 * algorithm and the COCO platform).
 */
static coco_problem_t *PROBLEM;


/**
 * Evaluates the function and constraints.
 * Borg expects a single function that does both.
 */
static void evaluate_problem(double* vars, double* objs, double* consts) {
    int nvars = (int) coco_problem_get_dimension(PROBLEM);
    int nints = (int) coco_problem_get_number_of_integer_variables(PROBLEM);
    int nconsts = (int) coco_problem_get_number_of_constraints(PROBLEM);

    // Truncate integer variables
    // do we have to clone vars?
    for (int i = 0; i < nints; i++)
        vars[i] = floor(vars[i]);

    coco_evaluate_function(PROBLEM, vars, objs);

    if (nconsts > 0) {
        coco_evaluate_constraint(PROBLEM, vars, consts);
        // constraints in coco are satisfied if their values are <= 0.
        // constraints in the borg implementtion are satisfied if their values are = 0.
        for (int j = 0; j < nconsts; j++)
            consts[j] = max(consts[j], 0);
    }

}

/**
 * Compares two doubles for qsort.
 */
static int compare_doubles(const void *a, const void *b) {
    double arg1 = *(const double *) a;
    double arg2 = *(const double *) b;

    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

void experiment(const char *suite_name,
                const char *suite_options,
                const char *observer_name,
                const char *observer_options,
                coco_random_state_t *random_generator);

int main(void) {

    /* Change the log level to "warning" to get less output */
    coco_set_log_level("warning");

    printf("Running the Borg MOEA Experiment\n");
    fflush(stdout);

    // bbob-biobj has 55 functions.
    const int N_FUNCS = 55;
    for (int i = 0; i < THREADS; i++) {
        if (fork() == 0) {
            // sprintf the thread number into the result_folder for the observer options
            char observer_options[100];
            sprintf(observer_options, "result_folder: borg_adaptive/thread_%d algorithm_name: BorgMOEA", i);
            // put the function range into the suite options
            // apparently supposed to be 1-indexed, lol
            int lower_bound = (N_FUNCS * i) / THREADS + 1;
            int upper_bound = (N_FUNCS * (i + 1)) / THREADS;
            char suite_options[100];
            sprintf(suite_options, "function_indices: %i-%i", lower_bound, upper_bound);

            // initialize the random number generator
            coco_random_state_t *random_generator = coco_random_new(RANDOM_SEED + i);
            // call the experiment
            printf("Starting thread %d.\n", i);
            experiment("bbob-biobj", suite_options, "bbob-biobj", observer_options, random_generator);

            printf("Thread %d finished.\n", i);
            exit(0);
        }
        // Parent thread can just do nothing, the children will take care of the experiment
    }
    // Wait until all children are done
    while(wait(NULL) > 0);

    printf("Done!\n");
    fflush(stdout);

    return 0;
}

/**
 * Runs the Borg MOEA on a suite of problems.
 * @param suite_name Name of the suite (e.g. "bbob" or "bbob-biobj").
 * @param suite_options Options of the suite (e.g. "dimensions: 2,3,5,10,20 instance_indices: 1-5").
 * @param observer_name Name of the observer matching with the chosen suite (e.g. "bbob-biobj"
 * when using the "bbob-biobj-ext" suite).
 * @param observer_options Options of the observer (e.g. "result_folder: folder_name")
 */
void experiment(const char *suite_name,
                const char *suite_options,
                const char *observer_name,
                const char *observer_options,
                coco_random_state_t *random_generator) {

    coco_suite_t *suite;
    coco_observer_t *observer;

    /* Initialize the suite and observer. */
    suite = coco_suite(suite_name, "", suite_options);
    observer = coco_observer(observer_name, observer_options);

    /* Iterate over all problems in the suite */
    while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {

        // Get the number of variables and objectives. 
        // Constraints is apparently zero for all coco suites, but we'll get it too
        int dimension = (int) coco_problem_get_dimension(PROBLEM);
        int objectives = (int) coco_problem_get_number_of_objectives(PROBLEM);
        int constraints = (int) coco_problem_get_number_of_constraints(PROBLEM);
        int int_vars = (int) coco_problem_get_number_of_integer_variables(PROBLEM);

        // Initialize the BORG problem
        BORG_Problem bproblem = BORG_Problem_create(dimension, objectives, constraints, evaluate_problem);

        // Set bounds for variables in BORG
        double* lower = coco_problem_get_smallest_values_of_interest(PROBLEM);
        double* upper = coco_problem_get_largest_values_of_interest(PROBLEM);
        for (int i = 0; i < dimension; i++) {
            // integer variables need to have their upper bound increased by 0.99
            // then, when we evaluate the coco problem, we'll take the floor. 
            // this way each integer has an equal amount of space in the continuous interval.
            double this_upper = upper[i];
            if (i < int_vars)
                this_upper += 0.99;
            BORG_Problem_set_bounds(bproblem, i, lower[i], this_upper);
        }

        // Estimate epsilon values
        // We do this by a random search to estimate the scale of each objective
        int n_epsilon_evals = dimension * EPSILON_EVALUATIONS;
        double* random_values = (double*) malloc(objectives * n_epsilon_evals * sizeof(double));
        double* vars = (double*) malloc(dimension * sizeof(double));
        double* objs = (double*) malloc(objectives * sizeof(double));

        for (int i = 0; i < n_epsilon_evals; i++) {
            // Pick a random point in the search space
            for (int j = 0; j < dimension; j++) {
                vars[j] = lower[j] + (upper[j] - lower[j]) * coco_random_uniform(random_generator);
            }
            coco_evaluate_function(PROBLEM, vars, objs);
            // copy the objectives into the random_values array
            // such that the first N values are the first objective, etc.
            for (int j = 0; j < objectives; j++) {
                random_values[j * n_epsilon_evals + i] = objs[j];
            }
        }
        // Sort each objective
        // Note that this is O(n log n) for each objective
        // and the median-of-medians algorithm is O(n) for each objective
        // but we don't care because we're only doing it once per problem
        for (int i = 0; i < objectives; i++) {
            double* this_obj = random_values + i * n_epsilon_evals;
            qsort(this_obj, n_epsilon_evals, sizeof(double), compare_doubles);

            // Set the epsilon values.
            // BORG requires that epsilon values be explicitly set.
            // The use of the epsilon values is that the algorithm will periodly restart
            // if it hasn't improved by at least epsilon in the last N iterations.
            // and epsilon-dominance is used to determine if a solution is a "duplicate".
    //
            // We'll set epsilon to be the 1st quartile of the objective values minus the min, times a constant
            double quart = this_obj[n_epsilon_evals / 4];
            double min = this_obj[0];
            BORG_Problem_set_epsilon(bproblem, i, (quart - min) * EPSILON_SCALE);
        }
        // free memory used in epsilon estimation
        free(vars);
        free(objs);
        free(random_values);

        // Initialize the Borg algorithm
        // create operators with default settings
        BORG_Operator pm = BORG_Operator_create("PM", 1, 1, 2, BORG_Operator_PM);
        BORG_Operator_set_parameter(pm, 0, 1.0 / dimension);
        BORG_Operator_set_parameter(pm, 1, 20.0);

        BORG_Operator sbx = BORG_Operator_create("SBX", 2, 2, 2, BORG_Operator_SBX);
        BORG_Operator_set_parameter(sbx, 0, 1.0);
        BORG_Operator_set_parameter(sbx, 1, 15.0);
        BORG_Operator_set_mutation(sbx, pm);

        BORG_Operator de = BORG_Operator_create("DE", 4, 1, 2, BORG_Operator_DE);
        BORG_Operator_set_parameter(de, 0, 0.1);
        BORG_Operator_set_parameter(de, 1, 0.5);
        BORG_Operator_set_mutation(de, pm);

        BORG_Operator um = BORG_Operator_create("UM", 1, 1, 1, BORG_Operator_UM);
        BORG_Operator_set_parameter(um, 0, 1.0 / dimension);

        BORG_Operator spx = BORG_Operator_create("SPX", 10, 2, 1, BORG_Operator_SPX);
        BORG_Operator_set_parameter(spx, 0, 3.0);

        BORG_Operator pcx = BORG_Operator_create("PCX", 10, 2, 2, BORG_Operator_PCX);
        BORG_Operator_set_parameter(pcx, 0, 0.1);
        BORG_Operator_set_parameter(pcx, 1, 0.1);

        BORG_Operator undx = BORG_Operator_create("UNDX", 10, 2, 2, BORG_Operator_UNDX);
        BORG_Operator_set_parameter(undx, 0, 0.5);
        BORG_Operator_set_parameter(undx, 1, 0.35);

        // Create the algorithm with specified operators
        // Also reduce the initial population size
        BORG_Algorithm algorithm = BORG_Algorithm_create(bproblem, 6);
        //BORG_Algorithm_set_initial_population_size(algorithm, dimension * 10);
        BORG_Algorithm_set_operator(algorithm, 0, sbx);
        BORG_Algorithm_set_operator(algorithm, 1, de);
        BORG_Algorithm_set_operator(algorithm, 2, pcx);
        BORG_Algorithm_set_operator(algorithm, 3, spx);
        BORG_Algorithm_set_operator(algorithm, 4, undx);
        BORG_Algorithm_set_operator(algorithm, 5, um);

        // Inject the initial solution
        double* x_init = coco_allocate_vector(dimension);
        coco_problem_get_initial_solution(PROBLEM, x_init);

        BORG_Solution initial_solution = BORG_Solution_create(bproblem);
        BORG_Solution_set_variables(initial_solution, x_init);
        // evaluate = 1, nfe = 1. 
        BORG_Algorithm_inject(algorithm, initial_solution, 1, 1);

        int evaluations_budget = dimension * BUDGET_MULTIPLIER;
        while (BORG_Algorithm_get_nfe(algorithm) < evaluations_budget) {
            BORG_Algorithm_step(algorithm);
        }

        // don't bother getting the result
        // BORG_Archive result = BORG_Algorithm_get_result(algorithm);

        /* Warn if the MOEA didn't use all its evaluations or an unexpected thing happened */
        long evaluations_done = coco_problem_get_evaluations(PROBLEM);
        if (evaluations_done < evaluations_budget) {
            printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n",
                   (unsigned long) evaluations_done, (unsigned long) evaluations_budget);
        }

        // COCO will keep track of the results, so we can immediately free BORG's memory.
        BORG_Operator_destroy(sbx);
        BORG_Operator_destroy(de);
        BORG_Operator_destroy(pm);
        BORG_Operator_destroy(um);
        BORG_Operator_destroy(spx);
        BORG_Operator_destroy(pcx);
        BORG_Operator_destroy(undx);
        BORG_Algorithm_destroy(algorithm);
        BORG_Problem_destroy(bproblem);
        // Does destroying the algorithm also destroy the initial solution?
        //BORG_Solution_destroy(initial_solution);
        coco_free_memory(x_init);

    }

    coco_observer_free(observer);
    coco_suite_free(suite);

}

