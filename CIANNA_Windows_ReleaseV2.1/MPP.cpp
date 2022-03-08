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
           set<int> id_Day;
	   for(int opt = 0; opt < N_OPT_DAY; opt++)
	     id_Day.insert(v_opt_dishes[opt][data[day*N_OPT_DAY+opt]].description);
	   int rd = random_dish(which);
           int id = v_opt_dishes[which][rd].description;
	   while(id_Day.find(id) != id_Day.end())//it forces to different dishes in a day
	   { 
	      rd = random_dish(which);
	      id = v_opt_dishes[which][rd].description;
           }
	   data[day * N_OPT_DAY + which] = rd;
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
 *This function computes the distance to closest neightbour, note that it also takes into account the repetition of dcn
 *Therefore if there are several dcn's values equal it also count the number of repetition.
 *The input is {distance to nearest 'day' of a dish, number of occurrences}
 *With this metric the best vaues the maximum
 * */
double MPP::normalize_dcn(pair<int, int> data_dcn){ 
	return ((double)data_dcn.first + ( 1.0 - ((double)data_dcn.second/(double)nDias)));
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
	map<Food, int> f1;
	int dist = 0;
	for (int i = 0; i < nDias; i++){
		Food f(N_OPT_DAY);
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = x_var[i*N_OPT_DAY+j];
		f1[f]++;
	}
	for (int i = 0; i < nDias; i++){
		Food f(N_OPT_DAY);
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = ind2.x_var[i*N_OPT_DAY+j];
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
void MPP::inc_eval(struct Solution_LS &current, Neighbor &new_neighbor, vector<double> &new_objs, vector<pair<double, double> > &newVariability){
    int num_nutr = (int)v_constraints.size();
    int day =  new_neighbor.variable/N_OPT_DAY;
    int opt = new_neighbor.variable%N_OPT_DAY;
    double new_partial_infeasibility_day = 0.0, original_partial_infeasibility_day = 0.0 ;
    double new_partial_infeasibility_global = 0.0, original_partial_infeasibility_global = 0.0 ;
   for(int b = 0; b < g_Idtime2Configs[opt].size(); b++)
   {
    int a = g_Idtime2Configs[opt][b];
    for(unsigned int j = 0; j < num_nutr; j++)
    {
       	//update sumatory of nutriments....
	double new_nut_value = (-v_opt_dishes[opt][current.x_var[new_neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][new_neighbor.newValue].v_nutrient_value[j]);
          double minv = v_constraints[j].min;
          double maxv = v_constraints[j].max;
	  double middle = (maxv+minv)*0.5;
       if(v_constraints[j].type == DIARIA)
       {
	  double nut = current.nutriment_per_day[a][day][j] + new_nut_value;
  	  double original_nut = current.nutriment_per_day[a][day][j];
	  if( nut < minv)new_partial_infeasibility_day+= ((minv- nut)/middle)*((minv - nut)/middle)*WEIGHT_DAY;
	  else if (nut > maxv) new_partial_infeasibility_day+=((nut - maxv)/middle)*((nut - maxv)/middle)*WEIGHT_DAY;
	  if( original_nut  < minv)original_partial_infeasibility_day += ((minv - original_nut)/middle)*((minv - original_nut)/middle)*WEIGHT_DAY;
	  else if (original_nut > maxv) original_partial_infeasibility_day +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle)*WEIGHT_DAY;
       }
       else if(v_constraints[j].type == GLOBAL)
       {
          double nut = current.globalPlan[a][j] + new_nut_value;
          double original_nut = current.globalPlan[a][j];
          if(nut < minv) new_partial_infeasibility_global+= ((minv - nut)/(middle))*((minv - nut)/(middle));
          else if (nut > maxv) new_partial_infeasibility_global+=((nut - maxv)/middle)*((nut - maxv)/middle);
          if(original_nut < minv) original_partial_infeasibility_global += ((minv - original_nut)/(middle))*((minv- original_nut)/(middle));
          else if ( original_nut > maxv) original_partial_infeasibility_global +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle);
       }
     }
   }
   new_objs[0]  = current.objFeasibility[0] - original_partial_infeasibility_day + new_partial_infeasibility_day;
   new_objs[1]  = current.objFeasibility[1] - original_partial_infeasibility_global + new_partial_infeasibility_global;

   if(new_objs[0] != current.objFeasibility[0]) return; //kind of branch procedure....
   if(new_objs[1] != current.objFeasibility[1]) return;
   //variability... this code-part will be optimized...
   int tmp = current.x_var[new_neighbor.variable];
   current.x_var[new_neighbor.variable] = new_neighbor.newValue;
   calculateVariability(current.x_var, newVariability);
   current.x_var[new_neighbor.variable] = tmp;
}
void MPP::update_inc(struct Solution_LS &current, Neighbor &neighbor, vector<double> &newFeas, vector<pair<double, double> > &newVar){
    int num_nutr = (int)v_constraints.size();
    int day =  neighbor.variable/N_OPT_DAY;
    int opt = neighbor.variable%N_OPT_DAY;
    struct infoDishes &dish_in = v_opt_dishes[opt][neighbor.newValue], &dish_out = v_opt_dishes[opt][current.x_var[day*N_OPT_DAY + opt]];
   for(int b = 0; b < g_Idtime2Configs[opt].size(); b++)
   {
    int a = g_Idtime2Configs[opt][b];
    for(unsigned int j = 0; j < num_nutr; j++)//check new neighbor..
    {
       	//update sumatory of nutriments....
	double diff = (-v_opt_dishes[opt][current.x_var[neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][neighbor.newValue].v_nutrient_value[j]);
        current.nutriment_per_day[a][day][j] += diff;
	current.globalPlan[a][j] += diff;
    }
   }
   current.x_var[neighbor.variable]= neighbor.newValue;
   current.objFeasibility = newFeas;
   current.objVariability = newVar;
   current.uniq_per_day[day][dish_in.description]++, current.uniq_per_day[day][dish_out.description]--;
}
void MPP::init_incremental_evaluation(struct Solution_LS &current){
        //feasibility information
	int num_nutr = (int)v_constraints.size();
   	current.uniq_per_day.assign(nDias+1, vector<int> (max_description_id+1, 0));
        current.globalPlan.assign((int)g_timesIdPerConf.size(),vector<double> (num_nutr, 0.0 ));
	current.nutriment_per_day.assign( (int)g_timesIdPerConf.size(), vector<vector<double> > (nDias, vector<double> (num_nutr, 0)));
	current.objFeasibility.assign((int)objFeasibility.size(), 0.0);
	current.objVariability.assign((int)objVariability.size(), pair<double, double>());
        for(int a = 0; a < g_timesIdPerConf.size(); a++){
 	  for(int j = 0; j < num_nutr; j++){
	   for(int i = 0; i < nDias; i++){
		 for(int b = 0; b < g_timesIdPerConf[a].size(); b++){
		   int k = g_timesIdPerConf[a][b];
	    	   current.nutriment_per_day[a][i][j] += v_opt_dishes[k][current.x_var[i*N_OPT_DAY + k]].v_nutrient_value[j];
		 }
	        current.globalPlan[a][j] += current.nutriment_per_day[a][i][j];
	        if( v_constraints[j].type == DIARIA){
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = current.nutriment_per_day[a][i][j];
	           if( nut < minv) current.objFeasibility[0] += ((minv - nut)/middle)*((minv - nut)/middle)*WEIGHT_DAY;
	           else if (nut > maxv) current.objFeasibility[0] +=((nut - maxv)/middle)*((nut - maxv)/middle)*WEIGHT_DAY;
                }
	   }
	   if( v_constraints[j].type == GLOBAL){
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = current.globalPlan[a][j];
	           if( nut < minv) current.objFeasibility[1] += ((minv - nut)/middle)*((minv - nut)/middle);
	           else if (nut > maxv) current.objFeasibility[1] += ((nut - maxv)/middle)*((nut - maxv)/middle);
            }
	  }
	}
        calculateVariability(current.x_var, current.objVariability);
}

void MPP::init_inc_eval_day(vector<int> &current_solution, vector<vector<double>> &nut_info)
{
   nut_info.assign((int)g_timesIdPerConf.size(), vector<double>(num_nutr, 0.0));
   for(int a = 0; a < g_timesIdPerConf.size(); a++)
   {
     for(int i = 0; i < v_constraint_day.size(); i++)
     {
       int j =  v_constraint_day[i];
       for(int b = 0; b < g_timesIdPerConf[a].size(); b++)
       {
         int k = g_timesIdPerConf[a][b];
         nut_info[a][j] += v_opt_dishes[k][current_solution[k]].v_nutrient_value[j];
       }
     }
   }

}
double MPP::inc_eval_day(Neighbor &new_neighbor, vector<vector<double>> &nut_info, vector<int> &current_solution, double current_feasibility)
{
   double new_partial_infeasibility_day = 0.0, original_partial_infeasibility_day = 0.0 ;
   int opt = new_neighbor.variable;
   for(int b = 0; b < g_Idtime2Configs[opt].size(); b++)
   {
       int a = g_Idtime2Configs[opt][b];
       for(int i = 0; i < v_constraint_day.size(); i++)
       {
          int j =  v_constraint_day[i];
         //update sumatory of nutriments....
         double new_nut_value = (-v_opt_dishes[opt][current_solution[opt]].v_nutrient_value[j] + v_opt_dishes[opt][new_neighbor.newValue].v_nutrient_value[j]);
         double minv = v_constraints[j].min;
         double maxv = v_constraints[j].max;
         double middle = (maxv+minv)*0.5;
         double nut = nut_info[a][j] + new_nut_value;
         double original_nut = nut_info[a][j];

         if( nut < minv)new_partial_infeasibility_day+= ((minv- nut)/middle)*((minv - nut)/middle);
         else if (nut > maxv) new_partial_infeasibility_day+=((nut - maxv)/middle)*((nut - maxv)/middle);
         if( original_nut  < minv)original_partial_infeasibility_day += ((minv - original_nut)/middle)*((minv - original_nut)/middle);
         else if (original_nut > maxv) original_partial_infeasibility_day +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle);
       }
   }
  return current_feasibility -original_partial_infeasibility_day + new_partial_infeasibility_day;
}
void MPP::update_inc_day(vector<vector<double>> &nut_info, Neighbor &new_neighbor, vector<int> &current_solution)
{
   int opt = new_neighbor.variable;
   for(int b = 0; b < g_Idtime2Configs[opt].size(); b++)
   {
       int a = g_Idtime2Configs[opt][b];
       for(int i = 0; i < v_constraint_day.size(); i++)
       {
          int j =  v_constraint_day[i];
         nut_info[a][j] = nut_info[a][j]-v_opt_dishes[opt][current_solution[new_neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][new_neighbor.newValue].v_nutrient_value[j];
       }
   }
   current_solution[opt] = new_neighbor.newValue;
}

