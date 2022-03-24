/*
 * Developers: Joel Chacón Castillo, Carlos Segura
 * Last modification: 03/05/22
 * The local search variants are defined in "localsearches.cpp"
 * The parsing of .xls files is defined in "parsingFiles.cpp"
 * The model that is taken into account as the fitness is defined in "Model.cpp"
 * */
#include <sys/time.h>
#include <string>
#include "MPP.h"
#include "global.h"
using namespace std;
MPP_Problem* MPP::MPP_problem;
extern vector<vector<int>> g_Idtime2Configs;
extern vector<vector<int> > g_timesIdPerConf;
extern vector<vector<int>> g_time2Ids;
extern vector<double> g_weights;
extern int g_crossoverType;

ostream & operator << (ostream &out, const vector<double> &data)
{
   for(int i = 0; i < data.size(); i++){
	   out << data[i]<< " ";
   }
   return out;
}

void MPP::init_neighbours(){
     for (int i = 0; i < (int)nDias; i++){
	for (int j = 0; j <N_OPT_DAY; j++){
	   for (int k = 0; k < (int) v_opt_dishes[j].size(); k++){
	     Neighbor n;
	     n.variable = i * N_OPT_DAY + j;
	     n.newValue = k;
	     neighbors.push_back(n);
	   }
        }
     }
     for(int i = 0; i < (int)nDias; i++){
      for(int j = i+1; j < (int)nDias; j++){
         neighbors_swap.push_back({i, j});
      }
     }

//       for(int i = 0; i < (int)nDias; i++){
//      for(int j = i+1; j < (int)nDias; j++){
//          for(int ii= j+1; ii <nDias; ii++){
//             for(int jj=ii; jj<nDias; jj++){
//                neighbors_swap2.push_back({i,j,ii, jj});
//             }
//          }
//      }
//     }


     for (int j = 0; j <N_OPT_DAY; j++)
     {
       for (int k = 0; k < (int) v_opt_dishes[j].size(); k++)
       {
	     Neighbor n;
	     n.variable = j;
	     n.newValue = k;
	     oneDayneighbors.push_back(n);
       }
     }

}
void MPP::init(){

   v_opt_dishes=MPP_problem->v_opt_dishes;
   v_constraints=MPP_problem->v_constraints;
   dic_nut_id = MPP_problem->dic_nut_id;
   v_constraint_global=MPP_problem->v_constraint_global;
   v_constraint_day=MPP_problem->v_constraint_day;
   out_filename=MPP_problem->out_filename;
   max_description_id=MPP_problem->max_description_id;
   nDias=MPP_problem->nDias;
   num_nutr=MPP_problem->num_nutr; 
   x_var.assign(N_OPT_DAY*(int)nDias, 0);
   //obj_values : {feasibility_day, feasibility_global, {variability_dishes}, {variability_categories}}
   objFeasibility.assign(2, 0);
   objVariability.assign(2*N_TIMES, pair<double, double>());
   for (int i = 0; i < (int)nDias; i++){
      perturb_day(x_var, i);
   }
   evaluate();
   //Definition of neighbourhood in Local search
   init_neighbours(); 

}
/*
 *This function takes each day and reassign it randomly
 * */
void MPP::restart(){
	for (int i = 0; i < (int)nDias; i++){
	   perturb_day(x_var, i);
	}
	evaluate();
}
/*
 *This function reassign just one item (dish) randomly
  note that this new item is different from those that are already selected.
 * */
void MPP::perturb_opt(vector<int> &data, int day, int which){
	   
	   data[day * N_OPT_DAY + which] = random_dish(which);
	  
         //  set<int> id_Day;
	 //  for(int opt = 0; opt < N_OPT_DAY; opt++)
	 //    id_Day.insert(v_opt_dishes[opt][data[day*N_OPT_DAY+opt]].description);
	 //  int rd = random_dish(which);
         //  int id = v_opt_dishes[which][rd].description;
	 //  while(id_Day.find(id) != id_Day.end())//it forces to different dishes in a day
	 //  { 
	 //     rd = random_dish(which);
	 //     id = v_opt_dishes[which][rd].description;
         //  }
	 //  data[day * N_OPT_DAY + which] = rd;
}
/*
 *
 * */
void MPP::perturb_day(vector<int> &data, int day){ 		
	set<int> dish;
	for(int k = 0; k < N_OPT_DAY; k++) 
 	{
	   int rd = random_dish(k);
	   int id = v_opt_dishes[k][rd].description;
	   while(dish.find(id) != dish.end())
	   { 
	      rd = random_dish(k);
	      id = v_opt_dishes[k][rd].description;
	   }
	   dish.insert(id);
	   data[day*N_OPT_DAY + k] = rd;
	
	}
}
/*
 *
 * */
void MPP::update_dcn_pair(int diff, pair<int, int> &p_dcn){ 
   if(diff < p_dcn.first){
      p_dcn = make_pair(diff, 1);
   }
   else if(diff == p_dcn.first){
	   p_dcn.second++; 
   }
}
/*
 *
 * */
int MPP::getDistance(MPP &ind2) {
///  map<int, pair<int, int > >counter;
///       for(int j = 0; j < x_var.size(); j++){
///           int i=x_var[j];
///           counter[i] = {counter[i].first+1, (int) v_opt_dishes[j%N_OPT_DAY].size()};
///       }
///       for(int j = 0; j < ind2.x_var.size(); j++){
///           int i=ind2.x_var[j];
///           counter[i] = {counter[i].first-1, (int) v_opt_dishes[j%N_OPT_DAY].size()};
///       }
///       double dist2 = 0;
///          for(auto i:counter) dist2 +=fabs(i.second.first)/(double)i.second.second;
///        return dist2;


	map<Food, int> f1;
	int dist = 0;
	for (int i = 0; i < nDias; i++){
		Food f(N_OPT_DAY);
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = x_var[i*N_OPT_DAY+j], f.sP.insert(f.p[j]);
		f1[f]++;
	}
	for (int i = 0; i < nDias; i++){
		Food f(N_OPT_DAY);
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = ind2.x_var[i*N_OPT_DAY+j], f.sP.insert(f.p[j]);
		if (f1.count(f)){
			f1[f]--;
			if (f1[f] == 0){
				f1.erase(f);
			}
		} else {
			dist++;
		}
	}
	return dist;
}
void MPP::swap_days(vector<int> &data, int day1, int day2)
{
   for(int ii = 0; ii < g_Idtime2Configs.size(); ii++)
   {
      swap(data[day1*N_OPT_DAY + ii], data[day2*N_OPT_DAY + ii] );
   }
}

void MPP::print(ostream &os) const {
	for (int i = 0; i < x_var.size(); i++){
		os << x_var[i] << " ";
	}
	os << fitness <<endl;
}








