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
	map<Food, vector<int> > f1;
	vector<int> freeDays1, freeDays2;
	for (int i = 0; i < nDias; i++){
		Food f(N_OPT_DAY);
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = x_var[i*N_OPT_DAY+j], f.sP.insert(f.p[j]);
		f1[f].push_back(i);
	}
	for (int i = 0; i < nDias; i++){
	    Food f(N_OPT_DAY);
	    for(int j = 0; j < N_OPT_DAY; j++){
	       f.p[j] = i2.x_var[i*N_OPT_DAY+j], f.sP.insert(f.p[j]);
	    }
//	    if(f1.find(f)!=f1.end()){
	    if(f1.count(f)!=0){
 	     //  for(int j = 0; j < N_OPT_DAY; j++){
	     //    x_var[currentDay*N_OPT_DAY + j] = f.p[j];
	     //    i2.x_var[currentDay*N_OPT_DAY + j] = f.p[j];
	     //  }
//	       f1[f]--;
	       f1[f].pop_back();
	    if(f1[f].empty()){
	       f1.erase(f);
	     }
	    }else {
		pendingI2.push_back(f);
		freeDays2.push_back(i);
	    }
	}
	for(auto it:f1){
	   for(auto it2:it.second){
		pendingI1.push_back(it.first);
		freeDays1.push_back(it2);
	   }	   
	}
//	for (map<Food, int>::iterator it = f1.begin(); it != f1.end(); it++){
//		for (int j = 0; j < it->second; j++){
//			pendingI1.push_back(it->first);
//		}
//	}
	if (pendingI1.size() != pendingI2.size()){ cerr << "Error interno. PendingI1 != PendingI2" << endl; exit(-1); }
        
	random_shuffle(pendingI1.begin(), pendingI1.end());
	random_shuffle(pendingI2.begin(), pendingI2.end());
	for (int i = 0; i < pendingI1.size(); i++){
		Food f1 = pendingI1[i];
		Food f2 = pendingI2[i];
		if (rand()%2){
			swap(f1, f2);
		}
 		for(int j = 0; j < N_OPT_DAY; j++)
	        {
	           x_var[freeDays1[i]*N_OPT_DAY + j] = f1.p[j];
	           i2.x_var[freeDays2[i]*N_OPT_DAY + j] =f2.p[j];
		}
	}
///    vector<int>perm1((int)x_var.size()), perm2((int)x_var.size());
///    for(int i = 0; i < nDias; i++) perm1[i]=i;
///    for(int i = 0; i < nDias; i++) perm2[i]=i;
///    random_shuffle(perm1.begin(), perm1.end());
///    random_shuffle(perm2.begin(), perm2.end());
///    auto tmp1=x_var;
///    auto tmp2=i2.x_var;
///    for(int i = 0; i < nDias; i++){
///       for(int j = 0; j < N_OPT_DAY; j++){
///	   x_var[i*N_OPT_DAY+j]=tmp1[perm1[i]*N_OPT_DAY+j];
///	   i2.x_var[i*N_OPT_DAY+j]=tmp2[perm2[i]*N_OPT_DAY+j];
///       }
///    }
}
void MPP::dependentMutation(double pm){
}
