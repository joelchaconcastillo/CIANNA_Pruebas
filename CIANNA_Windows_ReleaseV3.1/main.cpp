#include <bits/stdc++.h>
using namespace std;
#include "MA.h"
#include "global.h"
extern int g_N;
extern double g_pc;
extern double g_pm;
extern double g_finalTime;
extern double g_ITV;
/*
   Input arguments:
   1) file of Dishes information
   2) file of Constraintment nutrients
   3) Days
   4) Time
   5) out file
*/
int main(int argc, char **argv){
	int seed = 11; //atoi(argv[6]);
	string saveFile="";//=string(argv[7]);
	g_pc=1.0;
	g_ITV=0.4;
//	g_pc=atof(argv[8]);
//	g_ITV=atof(argv[9]);
//	srand(time(NULL));
	srand(seed);
	g_finalTime = atof(argv[4])*60;//25 * 60;
	//loading information for the problem
	//note that MPP is a derived class of MPP_Problem.
	MPP STP;
        STP.load_data(argc, argv);
//	g_ITV = g_ITV*STP.nDias;
	MPP::MPP_problem = &STP;
        //everithing seem's to be OK, thus starting the algorithm...
	MA ma(g_N, g_pc, g_pm, g_finalTime, saveFile);
	//minimum amount of local searches that should be applied
	//if there is not enough time then it applies just local searches
	int minimumLS=LOWEST_NUMBER_OF_GENERATIONS*g_N; //ten times the population size
	ma.run(minimumLS);
}
