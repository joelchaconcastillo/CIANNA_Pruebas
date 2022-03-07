#ifndef H_GLOBAL
#define H_GLOBAL
#include <bits/stdc++.h>
using namespace std;
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
/***
 * USER PARAMETERS
 * */


#define DAYS_FAVORITE 7*3
#define DAYS_NO_FAVORITE 7*4

//Iterations applied inside the local search
#define ITERATIONS_LS 100
#define MAX_ITE_ONE_DAY_LS 200

//Weights that are taken into account in the model
//Penalization for unfulfilled constraint per day ->Daily constraints
#define WEIGHT_DAY 1.0e6  
//Penalization for unfunfilled constraint globally (all the days) -->Global constraints
#define W_VAR_GLOBAL 100
//Penalization rised by assigning a different category
#define W_VAR_GLOBAL_CAT 1

#define W_VAR_BREAKFAST 1
#define W_VAR_MORNING_SNACK 1
#define W_VAR_STARTER 1
#define W_VAR_MAIN_COURSE 1
#define W_VAR_EVENING_SNACK 1
#define W_VAR_DINNER 1


#endif
