#ifndef __MA_H__
#define __MA_H__
#include "MPP.h"

struct ExtendedIndividual {
        MPP ind;
	int dist;
};
class MA {
	public:
		MA(int N_, double pc_, double pm_, double finalTime_);
		void run(int minimumLS);
		void runMemetic();
	private:
		//Parameters of MA
		int N;//Population Size
		double pc;//crossover probability
		double pm;//mutation probability
		double finalTime;//Seconds

		//Basic procedures of MA
		void initPopulation();
		void initDI();
		void selectParents();
		void crossover();
		void mutation();
		void localSearch();
		void replacement();
		int getRandomInteger0_N(int n){ 
                   return (int) ((n + 1.0)*rand()/(RAND_MAX+1.0));
                }

                double generateRandomDouble0_Max(double maxValue){
	               return (double)(rand()) / RAND_MAX * maxValue;
                }
         public:
		//Internal attributes of MA
		vector< ExtendedIndividual * > population; 
		vector< ExtendedIndividual * > parents;
		vector< ExtendedIndividual * > offspring;
		double initialTime;
		double DI;
		double timeLS;
};
#endif
