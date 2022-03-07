#include <bits/stdc++.h>
using namespace std;
#include "MA.h"
#include "global.h"
extern int g_N;
extern double g_pc;
extern double g_pm;
extern double g_pm;
extern double g_finalTime;
/*
   Input arguments:
   1) file of Dishes information
   2) file of Constraintment nutrients
   3) Days
   4) Time
   5) out file
*/
int main(int argc, char **argv){
	srand(time(NULL));

	g_finalTime = atof(argv[4])*60;//25 * 60;
	//loading the input data...
	MPP STP;
        STP.load_data(argc, argv);
	MPP::MPP_problem = &STP;
        //everithing seem's to be OK, thus starting the algorithm...
	MA ma(g_N, g_pc, g_pm, g_finalTime);
	//minimum amount of local searches that should be applied
	//if there is not enough time then it applies just local searches
	int minimumLS=10*g_N; //ten times the population size
	ma.run(minimumLS);
}
