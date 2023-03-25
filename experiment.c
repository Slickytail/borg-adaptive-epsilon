/**
 * An experiment that runs the BORG MOEA on a COCO suite.
 */
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/wait.h>

#include "coco.h"
#include "borg.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

/**
 * Whether to print verbose output.
 */
static int VERBOSE = 0;

/**
 * The number of threads to use for the experiment.
 * Problems in a suite will be divided among threads and each thread will run solve problems synchronously.
 */
static int THREADS = 4;

/**
 * The maximal budget for evaluations done by an optimization algorithm equals dimension * BUDGET_MULTIPLIER.
 */
static long BUDGET_MULTIPLIER = 1e3;


/**
 * Epsilon scale multiplier for initial adaptation.
 */
static double ETA = 1e-4;


/**
 * Whether to inject the initial COCO solution
 */
static int INJECT_INITIAL_COCO = 0;


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
                const char *observer_options);

void adapt_epsilons(BORG_Algorithm algorithm,
                    BORG_Problem problem,
                    double multiplier,
                    int decrease_only);

void adapt_epsilons_hypervolume(
                    BORG_Algorithm algorithm,
                    BORG_Problem problem,
                    double multiplier,
                    int decrease_only);

int main(int argc, char **argv) {
    // Parse arguments
    int c;
    static struct option long_options[] = {
        {"verbose", no_argument, &VERBOSE, 1},
        {"inject", no_argument, &INJECT_INITIAL_COCO, 1},
        {"threads", required_argument, 0, 't'},
        {"budget", required_argument, 0, 'b'},
        {"eta", required_argument, 0, 'e'},
        {0, 0, 0, 0}
    };
    while (1) {
        int option_index = 0;
        c = getopt_long(argc, argv, "t:b:e:",
                        long_options, &option_index);
        if (c == -1) break;

        // switch to handle argument values
        // strtol sets errno if it fails
        errno = 0;
        switch (c) {
            case 0:
                break;
            case 't':
                THREADS = (int) strtol(optarg, NULL, 10);
                if (errno != 0 || THREADS < 1) {
                    printf("Invalid thread count: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'b':
                // we parse as a long double to correctly parse scientific notation
                BUDGET_MULTIPLIER = (long) strtold(optarg, NULL);
                if (errno != 0 || BUDGET_MULTIPLIER < 1) {
                    printf("Invalid budget: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'e':
                ETA = strtod(optarg, NULL);
                if (errno != 0 || ETA < 0) {
                    printf("Invalid eta: %s\n", optarg);
                    exit(1);
                }
                break;
            default:
                abort();
        }
    }
    // End of argument parsing

    // info log level is more verbose than warning
    if (VERBOSE)
        coco_set_log_level("info");
    else
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

            // call the experiment
            printf("Starting thread %d.\n", i);
            experiment("bbob-biobj", suite_options, "bbob-biobj", observer_options);

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
                const char *observer_options) {

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
        // Set default epsilon values
        // If we try to initialize the algorithm before the epsilons are set, borg will complain
        for (int i = 0; i < objectives; i++) {
            BORG_Problem_set_epsilon(bproblem, i, 1e-4);
        }


        // Initialize the Borg algorithm.
        // This initialization routine and parameters is taken from BORG_Algorithm_run
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
        BORG_Algorithm algorithm = BORG_Algorithm_create(bproblem, 6);
        // Set the initial population size
        //BORG_Algorithm_set_initial_population_size(algorithm, 10 * dimension);
        BORG_Algorithm_set_operator(algorithm, 0, sbx);
        BORG_Algorithm_set_operator(algorithm, 1, de);
        BORG_Algorithm_set_operator(algorithm, 2, pcx);
        BORG_Algorithm_set_operator(algorithm, 3, spx);
        BORG_Algorithm_set_operator(algorithm, 4, undx);
        BORG_Algorithm_set_operator(algorithm, 5, um);

        // Estimate epsilon values
        // We run the first iteration of the algorithm, which will cause it to generate an initial population.
        BORG_Algorithm_step(algorithm);
        // Which method to adapt?
        adapt_epsilons(algorithm, bproblem, ETA, 0);

        if (INJECT_INITIAL_COCO) {
            // Inject the initial solution from COCO
            double* x_init = coco_allocate_vector(dimension);
            coco_problem_get_initial_solution(PROBLEM, x_init);

            BORG_Solution initial_solution = BORG_Solution_create(bproblem);
            BORG_Solution_set_variables(initial_solution, x_init);
            // evaluate = 1, nfe = 1. 
            BORG_Algorithm_inject(algorithm, initial_solution, 1, 1);
            // free the memory storing the coco initial solution
            // this is safe because borg has already copied it into its initial solution
            coco_free_memory(x_init);
        }

        // Loop the main algorithm iterate
        long evaluations_budget = dimension * BUDGET_MULTIPLIER;
        while (BORG_Algorithm_get_nfe(algorithm) < evaluations_budget) {
            BORG_Algorithm_step(algorithm);
        }

        // don't bother getting the result
        // BORG_Archive result = BORG_Algorithm_get_result(algorithm);

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

    }

    coco_observer_free(observer);
    coco_suite_free(suite);

}

void adapt_epsilons(BORG_Algorithm algorithm, BORG_Problem problem, double multiplier, int decrease_only) {
    // TODO: I've modified the BORG header to allow direct access to some data structures.
    // It should be possible to use the borg api to do everything in here, but it's definitely more annoying.
    BORG_Population population = algorithm->population;
    int popsize = population->size;
    double* initial_population_objective = (double *) malloc(popsize * sizeof(double));
    for (int obj = 0; obj < problem->numberOfObjectives; obj++) {
        // We copy the borg population to a buffer that we can sort
        for (int i = 0; i < popsize; i++)
            initial_population_objective[i] = population->members[i]->objectives[obj];
        // Sort the objective values
        qsort(initial_population_objective, popsize, sizeof(double), compare_doubles);

        // We'll set epsilon to be the 1st quartile of the objective values minus the min, times a constant
        double quart = initial_population_objective[popsize / 4];
        double min = initial_population_objective[0];
        double new_epsilon = (quart - min) * multiplier;
        // If decrease_only is set, only apply the update if the new epsilon is smaller
        // This ensures that the new epsilon-dominance relation is weaker than the old one,
        // so it remains satisfied.
        if ((new_epsilon < problem->epsilons[obj]) || !decrease_only)
            BORG_Problem_set_epsilon(problem, obj, new_epsilon);
    }
}

void adapt_epsilons_hypervolume(BORG_Algorithm algorithm, BORG_Problem problem, double multiplier, int decrease_only) {
    // Compute the max and min for the pareto front from the BORG archive
    int nobj = problem->numberOfObjectives;
    double* max = (double *) malloc(nobj * sizeof(double));
    double* min = (double *) malloc(nobj * sizeof(double));
    for (int i = 0; i < nobj; i++) {
        max[i] = -INFINITY;
        min[i] = INFINITY;
    }
    BORG_Archive archive = BORG_Algorithm_get_result(algorithm);
    int size = BORG_Archive_get_size(archive);
    for (int i = 0; i < size; i++) {
        BORG_Solution solution = BORG_Archive_get(archive, i);
        for (int j = 0; j < nobj; j++) {
            double obj = BORG_Solution_get_objective(solution, j);
            if (obj > max[j])
                max[j] = obj;
            if (obj < min[j])
                min[j] = obj;
        }
    }
    for (int i = 0; i < nobj; i++) {
        double new_epsilon = (max[i] - min[i]) * multiplier;
        if ((new_epsilon < problem->epsilons[i]) || !decrease_only)
            BORG_Problem_set_epsilon(problem, i, new_epsilon);
    }
}
