#ifndef H_GLOBAL
#define H_GLOBAL
#include <bits/stdc++.h>
using namespace std;


///SYSTEM TAGS ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
#define EPSILON 1e-10
#define N_CATEGORIES 3
#define CATEGORY_1 1
#define CATEGORY_2 2
#define CATEGORY_BOTH 0
#define GLOBAL 1
#define DIARIA 2

//v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner, v_both_snack;
//encoded times by each day of the individual
#define N_TIMES 6
#define N_OPT_DAY 8

//Options per time
//Breakfast
#define BREAKFAST_1 0
#define MORNING_SNACK_1 1
#define STARTER_1 2
#define STARTER_2 3
#define MAIN_COURSE_1 4
#define MAIN_COURSE_2 5
#define EVENING_SNACK_1 6
#define DINNER_1 7

///Times per day
#define BREAKFAST 0
#define MORNING_SNACK 1
#define STARTER 2
#define MAIN_COURSE 3
#define EVENING_SNACK 4
#define DINNER 5


#define PAIR_BASED_CROSSOVER 1
#define UNIFORM_CROSSOVER 2
#define UNIFORM2_CROSSOVER 3

//////////////////////////////////////////////////
///USER PARAMETERS TAGS
//These parameters should be fitted as is desired
//////////////////////////////////////////////////

//The lowest number of generations required to apply a memetic algorithm
//this numbers is calculated based on the time required to find a feasible solution
#define LOWEST_NUMBER_OF_GENERATIONS 10

//This factor extends time required to attain a feasible solution of the local search for each generation in the memetic algorithm,
//in other words it is multiplied by the lowest time in a local search to attain a feasible solution
//note that as this value is longer a high intensification for each local search might be applied
//and might be a good parameter to tunne for the most challenging instances
#define FACTOR_TIME_PER_LOCAL_SEARCH 10
//#define FACTOR_TIME_PER_LOCAL_SEARCH 100 //NOTE test this!!

//This defines the minimum days allowed for repetition
#define DAYS_FAVORITE 7*3
#define DAYS_NO_FAVORITE 7*4

//Iterations applied inside the local search
#define MAX_ITE_ONE_DAY_LS 200

///////////////////////////////////////////////////////
//Weights that are taken into account in the model NOTE check 1e6...
#define W_FEASIBILITY 1e20
//Penalization for unfulfilled constraint per day ->Daily constraints
#define WEIGHT_DAY 1.0e6

#endif
