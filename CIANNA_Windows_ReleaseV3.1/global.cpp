#include <bits/stdc++.h>
#include "global.h"
#define CONFIGURATION1 0
#define CONFIGURATION2 1

using namespace std;


//population size
int g_N=10;
//xover probability
double g_pc=0.2;
//mutation probability
double g_pm=0.01;
//Stopping criterion set by time in minutes
double g_finalTime=1;
double g_ITV=0.5;


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
//TIMES				   ||  l va || u val|
//========================================
//MORNING_SNACK, EVENING_SNACK     ||  1e3 1e4
//BREAKFAST, DINNER		   ||  1e4 1e5
//STARTER                          ||  1e5 1e6
//MAIN_COURSE                      ||  1e6 1e7
///////////////////////////////////////
//
vector<double> g_l_rank={1e4, 1e3, 1e5, 1e6, 1e3, 1e4};
vector<double> g_u_rank={1e5, 1e4, 1e6, 1e7, 1e4, 1e5};

//This is applied in a similar way than the previous rank, instead this rank is lower resulting in a less relevance principle.
//TIMES				   ||  Lowest val || upper val
//========================================
//MORNING_SNACK, EVENING_SNACK     ||  *0   1   make sure that there is not a tradeoff between weights
//BREAKFAST, DINNER		   ||  *1 10
//STARTER                          ||  *10 100
//MAIN_COURSE                      ||  *100 1000
///////////////////////////////////////
vector<double> g_l_rank_category={1, 0, 10, 1e2, 0, 1};
vector<double> g_u_rank_category={10, 1, 1e2, 1e3, 1, 10};

