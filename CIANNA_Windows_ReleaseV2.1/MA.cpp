/*
 * Developers: Carlos Segura, Joel Chacón
 * Last modification: 03/05/22
 * Description:
 *  The MA class refers to "Memetic Algorithm" which consists of a genetic algorithm and a local search procedure.
 * */

#include <sys/time.h>
#include <iostream>
#include <signal.h>
#include "MA.h"
#include "global.h"

using namespace std;
/*
 *  This constructor recieves
 *  N_: Population size
 *  pc_: crossover probability
 *  pm_: mutation probability
 *  finalTime_: total time of execution
 * */
MA::MA(int N_, double pc_, double pm_, double finalTime_){
	if (N_ % 2){ cerr << "El tam. de poblacion debe ser par" << endl; exit(-1); }
	N = N_;
	pc = pc_;
	pm = pm_;
	finalTime = finalTime_;
	struct timeval currentTime; 
	gettimeofday(&currentTime, NULL);
	initialTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
}
/*
 *This function randomly initializes the population with size N
 * */
void MA::initPopulation(){
	for (int i = 0; i < N; i++){
		ExtendedIndividual *ei = new ExtendedIndividual();
		ei->ind.init();
		population.push_back(ei);	
	}
}
/*
 *In this functions is carried out a binary tournament selection,
 *nota that it considers minimization.
 * */
void MA::selectParents(){
	parents.clear();
	for (int i = 0; i < N; i++){
		int first = getRandomInteger0_N(N - 1);
		int second = getRandomInteger0_N(N - 1);
		if (population[first]->ind.fitness <= population[second]->ind.fitness){
			parents.push_back(population[first]);
		} else {
			parents.push_back(population[second]);
		}
	}
}
/*
 *This function applyes a crossover which is defined in the class Individual,
 *in particular the dependent crossover is taken into account.
 *Note that "ExtendedIndividual" is  just a structure of one Individual and a distance.
 * */
void MA::crossover(){
	for (int i = 0; i < parents.size(); i++){
		ExtendedIndividual *ei = new ExtendedIndividual();
		*ei = *parents[i];
		offspring.push_back(ei);
	}
	for (int i = 0; i < offspring.size(); i+=2){
		if (generateRandomDouble0_Max(1) <= pc){
			offspring[i]->ind.dependentCrossover(offspring[i+1]->ind);
		}
	}
}
/*
 *In this function is applied a mutation which should be defined in the Individual class.
 * */
void MA::mutation(){
	for (int i = 0; i < offspring.size(); i++){
		offspring[i]->ind.dependentMutation(pm);
	}
}
/*
 *This function applies a local search for each offspring
 * */
void MA::localSearch(){
	for (int i = 0; i < offspring.size(); i++){
		struct timeval currentTime; 
		gettimeofday(&currentTime, NULL);
		double elapsedTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
		double availableTime = min(finalTime-(elapsedTime-initialTime), timeLS);
		offspring[i]->ind.localSearch(availableTime);
	}
}
/*
 *This function applies the Best-Non-Penalized (BNP) procedure,
 *that explicitly considers diversity in the decision variable space.
 * */
void MA::replacement(){
	vector < ExtendedIndividual* > all;
	
	//Join population and offspring multisets
	for (int i = 0; i < population.size(); i++){
		all.push_back(population[i]);
		all.back()->dist = INT_MAX;
	}
	population.clear();
	for (int i = 0; i < offspring.size(); i++){
		all.push_back(offspring[i]);
		all.back()->dist = INT_MAX;
	}
	offspring.clear();
	
	//Select best solution, which has the lowest fitness
	int indexBest = 0;
	for (int i = 1; i < all.size(); i++){
		if (all[i]->ind.fitness < all[indexBest]->ind.fitness){
			indexBest = i;
		}
	}
	population.push_back(all[indexBest]);
	all[indexBest] = all.back(); 
	all.pop_back();

	struct timeval currentTime; 
	gettimeofday(&currentTime, NULL);
	double elapsedTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
	elapsedTime -= initialTime;

	double D = DI - DI * elapsedTime / finalTime;
	while(population.size() != N){
		//Update distances
		for (int i = 0; i < all.size(); i++){
			all[i]->dist = min(all[i]->dist, all[i]->ind.getDistance(population.back()->ind));
		}
		//Select best option
		indexBest = 0;
		for (int i = 1; i < all.size(); i++){
			bool betterInDist =	(all[i]->dist > all[indexBest]->dist);
			bool eqInDist = (all[i]->dist == all[indexBest]->dist);
			bool betterInFit = (all[i]->ind.fitness < all[indexBest]->ind.fitness);
			bool eqInFit = (all[i]->ind.fitness == all[indexBest]->ind.fitness);
			if (all[indexBest]->dist < D){//Do not fulfill distance requirement
				if ((betterInDist) || (eqInDist && betterInFit)){
					indexBest = i;
				}
			} else {
				if (all[i]->dist >= D){
					if ((betterInFit) || (eqInFit && betterInDist)){
						indexBest = i;
					}
				}
			}
		}
		//Insert best option
		population.push_back(all[indexBest]);
		all[indexBest] = all.back();
		all.pop_back();
	}
	//Release memory
	for (int i = 0; i < all.size(); i++){
		delete(all[i]);
	}
}
/*
 *The initial distance factor is considered as the average distance between individuals.
 * */
void MA::initDI(){
	double meanDistance = 0;
	for (int i = 0; i < population.size(); i++){
		for (int j = i + 1; j < population.size(); j++){
			meanDistance += population[i]->ind.getDistance(population[j]->ind);
		}
	}
	meanDistance /= (population.size() * (population.size() - 1)) / 2;
	DI = meanDistance * 1;
}
void MA::runMemetic(){
	initPopulation();
	initDI();
	int generation = 0;
	struct timeval currentTime; 
	gettimeofday(&currentTime, NULL);
	double elapsedTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
	elapsedTime -= initialTime;
	while(elapsedTime < finalTime){//Infinitas generaciones
		int minDistance = INT_MAX;
		for (int i = 0; i < population.size(); i++){
			for (int j = i + 1; j < population.size(); j++){
				minDistance = min(minDistance, population[i]->ind.getDistance(population[j]->ind));
			}
		}
		selectParents();
		crossover();
		mutation();
		localSearch();
		replacement();
		generation++;
		gettimeofday(&currentTime, NULL);
		elapsedTime = ((double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6)-initialTime;
	}
	int indexBest = 0;
	for (int i = 1; i < population.size(); i++)
        {
	   if(population[indexBest]->ind.fitness > population[i]->ind.fitness)
	   indexBest = i;
        }
	population[indexBest]->ind.exportcsv();

}
/*
 *This function 
 * */
void MA::run(int minimumLS){
        //get the time needed to attain a feasible solution..	
 	ExtendedIndividual *ei = new ExtendedIndividual();
	//init individual
	ei->ind.init();
	/*
	 *This local search measures the total time needed to calculate a feasible solution
	 *if time allows it then is applied a memetic algorithm otherwise a local search.
	 * */
  	//double ls_time = ei->ind.localSearch_testing_time(finalTime);
        bool isMeasuringTime=false;
  	double ls_time = ei->ind.localSearch(finalTime, isMeasuringTime);
        if(ls_time*minimumLS > finalTime)
	{
	  ei->ind.localSearch(finalTime-ls_time);
	  ei->ind.exportcsv();
	}
	else{
	  finalTime -=ls_time;
	  timeLS = 10*ls_time;
	  runMemetic();
	}
	delete ei;

}

