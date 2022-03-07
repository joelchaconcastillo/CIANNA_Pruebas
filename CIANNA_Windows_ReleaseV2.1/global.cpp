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
vector<double> g_weights={W_VAR_BREAKFAST, W_VAR_MORNING_SNACK, W_VAR_STARTER, W_VAR_MAIN_COURSE, W_VAR_EVENING_SNACK, W_VAR_DINNER};
