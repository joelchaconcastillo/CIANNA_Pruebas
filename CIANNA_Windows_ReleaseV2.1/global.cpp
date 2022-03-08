#include <bits/stdc++.h>
#include "global.h"
#define CONFIGURATION1 0
#define CONFIGURATION2 1

using namespace std;


//population size
int g_N=10;
//xover probability
double g_pc=1;
//mutation probability
double g_pm=0.01;
//Stopping criterion set by time in minutes
double g_finalTime=1;


/*
  Each of six times has several options of dishes 

    TIME          : {OPTIONS}
  =================================================
  BREAKFAST       : {BREAKFAST_1}
  MORNING_SNACK   : {MORNING_SNACK_1}
  STARTER         : {STARTER_1, STARTER_2}
  MAIN_COURSE     : {MAIN_COURSE_1, MAIN_COURSE_2}
  EVENING_SNACK   : {EVENING_SNACK_1}
  DINNER          : {DINNER_1}
  =================================================

  Each day is encoded by two configurations which are:
     configuration1 = {BREAKFAST_1, MORNING_SNACK_1, STARTER_1, MAIN_COURSE_1, EVENING_SNACK_1, DINNER_1}
     configuration2 = {BREAKFAST_1, MORNING_SNACK_1, STARTER_2, MAIN_COURSE_2, EVENING_SNACK_1, DINNER_1}
  A day is composed by six times: breakfast, morning_snack, starter, main_course, evening_snack and dinner.


  g_Idtime2Configs : stores for each time the available id's configuration, note that STARTERS and MAIN_COURSES have two options and remaining times only one.
                   In other words it maps from a time option to its configuration.
    TIME          : {CONFIGURATIONS}
  =================================================
  BREAKFAST_1       : {CONFIGURATION1, CONFIGURATION2}
  MORNING_SNACK_1   : {CONFIGURATION1, CONFIGURATION2}
  STARTER_1         : {CONFIGURATION1}
  STARTER_2         : {CONFIGURATION2}
  MAIN_COURSE_1     : {CONFIGURATION1}
  MAIN_COURSE_2     : {CONFIGURATION2}
  EVENING_SNACK_1   : {CONFIGURATION1, CONFIGURATION2}
  DINNER_1          : {CONFIGURATION1, CONFIGURATION2}
  =================================================

  g_time2Ids     : stores for each time the available id's options, similarly than with g_times2Configs  starters and main_courses have two options.
     TIME          : {OPTIONS}
  =================================================
  BREAKFAST       : {BREAKFAST_1}
  MORNING_SNACK   : {MORNING_SNACK_1}
  STARTER         : {STARTER_1, STARTER_2}
  MAIN_COURSE     : {MAIN_COURSE_1, MAIN_COURSE_2}
  EVENING_SNACK   : {EVENING_SNACK_1}
  DINNER          : {DINNER_1}
  =================================================
 
  g_timesIdPerConf : given the required configurations (default two), for each configuration fixs the six times

  TIME         : {CONFIGURATION1}         : {CONFIGURATION2}
  =================================================
  BREAKFAST       : {BREAKFAST_1}         : {BREAKFAST_1}
  MORNING_SNACK   : {MORNING_SNACK_1}     : {MORNING_SNACK_1}
  STARTER         : {STARTER_1}           : {STARTER_2}
  MAIN_COURSE     : {MAIN_COURSE_1}       : {MAIN_COURSE_2}
  EVENING_SNACK   : {EVENING_SNACK_1}     : {EVENING_SNACK_1}
  DINNER          : {DINNER_1}            : {DINNER_1}
  =================================================


 */
vector<vector<int>> g_Idtime2Configs={{CONFIGURATION1,CONFIGURATION2},{CONFIGURATION1,CONFIGURATION2},{CONFIGURATION1},{CONFIGURATION2},{CONFIGURATION1},{CONFIGURATION2},{CONFIGURATION1,CONFIGURATION2},{CONFIGURATION1,CONFIGURATION2}};

vector<vector<int>> g_time2Ids={{BREAKFAST_1}, {MORNING_SNACK_1}, {STARTER_1, STARTER_2}, {MAIN_COURSE_1, MAIN_COURSE_2}, {EVENING_SNACK_1}, {DINNER_1}};

vector<vector<int> > g_timesIdPerConf ={{BREAKFAST_1, MORNING_SNACK_1, STARTER_1, MAIN_COURSE_1, EVENING_SNACK_1, DINNER_1},{BREAKFAST_1, MORNING_SNACK_1, STARTER_2, MAIN_COURSE_2, EVENING_SNACK_1, DINNER_1}};

int g_crossoverType=PAIR_BASED_CROSSOVER;
/**
 * Weight applied to penalize variability
 * */
//This rank is added to each normalized score in variability terms, in this way there is a rank between times, this avoid tradeoff between times regarding variability
//A higher rank gives more relevance to the time since the problem is of minimization, note that it begins with 7 since there a different kind of varibility which is by categories.
//TIMES				   ||  rank
//========================================
//MORNING_SNACK, EVENING_SNACK     ||  7
//BREAKFAST, DINNER		   ||  8
//STARTER                          ||  9
//MAIN_COURSE                      ||  10
///////////////////////////////////////
vector<double> g_rank={8, 7, 9, 10, 7, 8};

//This is applied in a similar way than the previous rank, instead this rank is lower resulting in a less relevance principle.
//TIMES				   ||  rank
//========================================
//MORNING_SNACK, EVENING_SNACK     ||  1
//BREAKFAST, DINNER		   ||  2
//STARTER                          ||  3
//MAIN_COURSE                      ||  4
///////////////////////////////////////
vector<double> g_rank_category={2, 1, 3, 4, 1, 2};

