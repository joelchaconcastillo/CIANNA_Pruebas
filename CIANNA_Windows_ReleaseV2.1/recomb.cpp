#include "global.h"
#include "MPP.h"
extern int g_crossoverType;
void MPP::dependentCrossover(MPP &i2){
   if (g_crossoverType == PAIR_BASED_CROSSOVER){
	pairBasedCrossover(i2);
    } else if (g_crossoverType == UNIFORM_CROSSOVER){
	uniformCrossover(i2);
    } else if (g_crossoverType == UNIFORM2_CROSSOVER){
	uniform2Crossover(i2);
    }
      else
    {
      cout << "Operador de cruce desconocido "<<endl;
      exit(EXIT_FAILURE);
    }
}
void MPP::uniformCrossover(MPP &i2)
{
   for (int i = 0; i < nDias; i++)
   {
      if (rand() > (RAND_MAX / 2))
	for(int j = 0; j < N_OPT_DAY; j++)
	   swap(x_var[i*N_OPT_DAY+ j], i2.x_var[i*N_OPT_DAY + j]);
   }
}
void MPP::uniform2Crossover(MPP &i2){
   for (int i = 0; i < (int)x_var.size(); i++)
   {
	if (rand() > (RAND_MAX / 2))
		swap(x_var[i], i2.x_var[i]);
   }
}
void MPP::pairBasedCrossover(MPP &i2)
{
	vector<Food> pendingI1, pendingI2;
	map<Food, int> f1;
	for (int i = 0; i < nDias; i++){
		Food f;
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = x_var[i*N_OPT_DAY+j];
		f1[f]++;
	}
	for (int i = 0; i < nDias; i++){
		Food f;

		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = i2.x_var[i*N_OPT_DAY+j];
		if (f1.count(f)){//Comida en ambos
 			for(int j = 0; j < N_OPT_DAY; j++)
			{
			  x_var[i*N_OPT_DAY + j] = f.p[j];
			  i2.x_var[i*N_OPT_DAY + j] = f.p[j];
			}
			f1[f]--;
			if (f1[f] == 0){
				f1.erase(f);
			}
		} else {
			pendingI2.push_back(f);
		}
	}
	for (map<Food, int>::iterator it = f1.begin(); it != f1.end(); it++){
		for (int j = 0; j < it->second; j++){
			pendingI1.push_back(it->first);
		}
	}
	if (pendingI1.size() != pendingI2.size()){ cerr << "Error interno. PendingI1 != PendingI2" << endl; exit(-1); }
	random_shuffle(pendingI1.begin(), pendingI1.end());
	int next = nDias - pendingI1.size();
	for (int i = 0; i < pendingI1.size(); i++){
		Food f1 = pendingI1[i];
		Food f2 = pendingI2[i];
		if (rand() < RAND_MAX / 2.0){
			swap(f1, f2);
		}
 		for(int j = 0; j < N_OPT_DAY; j++)
	        {
	           x_var[i*N_OPT_DAY + j] = f1.p[j];
	           i2.x_var[i*N_OPT_DAY + j] =f2.p[j];
		}
		next++;
	}
}
void MPP::dependentMutation(double pm){
}
