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
#include "utils.h"
using namespace std;
MPP_Problem* MPP::MPP_problem;
ostream & operator << (ostream &out, const vector<double> &data)
{
   for(int i = 0; i < data.size(); i++){
	   out << data[i]<< " ";
   }
   return out;
}


void MPP::init(){
   x_var.assign(N_OPT_DAY*nDias, 0);
   obj_values.assign(N_TIMES+2, 0);
   for (int i = 0; i < nDias; i++){
      perturb_day(x_var, i);
   }
   evaluate();
}
/*
 *This function takes each day and reassign it randomly
 * */
void MPP::restart(){
	for (int i = 0; i < nDias; i++){
	   perturb_day(x_var, i);
	}
	evaluate();
}

void MPP::perturb_opt(vector<int> &data, int day, int which){
           set<int> id_Day;
	   for(int opt = 0; opt < N_OPT_DAY; opt++)
	     id_Day.insert(MPP_problem->v_opt_dishes[opt][data[day*N_OPT_DAY+opt]].description);
	   int rd = MPP_problem->random_dish(which);
           int id = MPP_problem->v_opt_dishes[which][rd].description;
	   while(id_Day.find(id) != id_Day.end())//it forces to different dishes in a day
	   { 
	      rd = MPP_problem->random_dish(which);
	      id = MPP_problem->v_opt_dishes[which][rd].description;
           }
	   data[day * N_OPT_DAY + which] = rd;
}
void MPP::perturb_day(vector<int> &data, int day){ 		
	set<int> dish;
	for(int k = 0; k < N_OPT_DAY; k++) 
 	{
	   int rd = MPP_problem->random_dish(k);
	   int id = MPP_problem->v_opt_dishes[k][rd].description;
	   while(dish.find(id) != dish.end())
	   { 
	      rd = MPP_problem->random_dish(k);
	      id = MPP_problem->v_opt_dishes[k][rd].description;
	   }
	   dish.insert(id);
	   data[day*N_OPT_DAY + k] = rd;
	
	}
}
double MPP::f(pair<int, int> data_dcn){ 
	return ((double)data_dcn.first + ( 1.0 - ((double)data_dcn.second/(double)nDias)));
}
void MPP::update_dcn_pair(int diff, pair<int, int> &p_dcn){ 
   if(diff < p_dcn.first){
      p_dcn = make_pair(diff, 1);
   }
   else if(diff == p_dcn.first){
	   p_dcn.second++; 
   }
}


void MPP::dependentCrossover(MPP &i2){
   if (crossoverType == PAIR_BASED_CROSSOVER){
	pairBasedCrossover(i2);
    } else if (crossoverType == UNIFORM_CROSSOVER){
	uniformCrossover(i2);
    } else if (crossoverType == UNIFORM2_CROSSOVER){
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
int MPP::getDistance(MPP &ind2) {
	map<Food, int> f1;
	int dist = 0;
	for (int i = 0; i < nDias; i++){
		Food f;
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = x_var[i*N_OPT_DAY+j];
		f1[f]++;
	}
	for (int i = 0; i < nDias; i++){
		Food f;
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
   for(int ii = 0; ii < MPP_problem->opt_conf.size(); ii++)
   {
      swap(data[day1*N_OPT_DAY + ii], data[day2*N_OPT_DAY + ii] );
   }
}
bool MPP::comp_objs(vector<double> &variability_v1, vector<double> &variability_v2)
{
  //feasibility
  if( variability_v1[0] < variability_v2[0]) return true;
  else if( variability_v1[0] > variability_v2[0]) return false;
  if( variability_v1[1] < variability_v2[1]) return true;
  else if( variability_v1[1] > variability_v2[1]) return false;

  //variability tchebycheff approach
  double max1 = -nDias, max2=-nDias;
  double ref = nDias*W_VAR_GLOBAL;
  for(int time_opt = 0; time_opt < N_TIMES; time_opt++)
  {
     double v1 = variability_v1[time_opt+2];
     double v2 = variability_v2[time_opt+2];
     max1 = max(max1,fabs(v1-(ref))/MPP_problem->weights[time_opt]);
     max2 = max(max2,fabs(v2-(ref))/MPP_problem->weights[time_opt]);
  }
 if( max1 < max2) return true;
 return false;
}
void MPP::print(ostream &os) const {
	for (int i = 0; i < x_var.size(); i++){
		os << x_var[i] << " ";
	}
	os << fitness <<endl;
}
void MPP::dependentMutation(double pm){
}
void MPP::evaluate(vector<int> &sol, vector<double> &objs){
  calculateFeasibilityDegree(sol, objs[0], objs[1]);
  calculateVariability(sol, objs);
}
void MPP::inc_eval(struct Solution_LS &current, Neighbor &new_neighbor, vector<double> &new_objs)
{
    int num_nutr = (int)MPP_problem->v_constraints.size();
    vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
    vector<int> &v_constraint_global = MPP_problem->v_constraint_global, &v_constraint_day = MPP_problem->v_constraint_day;
    vector<vector<infoDishes> > &v_opt_dishes = (MPP_problem->v_opt_dishes);
    int day =  new_neighbor.variable/N_OPT_DAY;
    int opt = new_neighbor.variable%N_OPT_DAY;
    double new_partial_infeasibility_day = 0.0, original_partial_infeasibility_day = 0.0 ;
    double new_partial_infeasibility_global = 0.0, original_partial_infeasibility_global = 0.0 ;
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {
    int a = MPP_problem->opt_conf[opt][b];
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
   new_objs[0]  = current.obj_values[0] - original_partial_infeasibility_day + new_partial_infeasibility_day;
   new_objs[1]  = current.obj_values[1] - original_partial_infeasibility_global + new_partial_infeasibility_global;

   if(new_objs[0] != current.obj_values[0]) return; //kind of branch procedure....
   if(new_objs[1] != current.obj_values[1]) return;
   //variability... this code-part will be optimized...
   int tmp = current.x_var[new_neighbor.variable];
   current.x_var[new_neighbor.variable] = new_neighbor.newValue;
   calculateVariability(current.x_var, new_objs);
   current.x_var[new_neighbor.variable] = tmp;
}
void MPP::update_inc(struct Solution_LS &current, Neighbor &neighbor, vector<double> &new_objs)
{
    int num_nutr = (int)MPP_problem->v_constraints.size();
    vector<vector<infoDishes> > &v_opt_dishes = (MPP_problem->v_opt_dishes);
    int day =  neighbor.variable/N_OPT_DAY;
    int opt = neighbor.variable%N_OPT_DAY;
    struct infoDishes &dish_in = MPP_problem->v_opt_dishes[opt][neighbor.newValue], &dish_out = MPP_problem->v_opt_dishes[opt][current.x_var[day*N_OPT_DAY + opt]];
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {
    int a = MPP_problem->opt_conf[opt][b];
    for(unsigned int j = 0; j < num_nutr; j++)//check new neighbor..
    {
       	//update sumatory of nutriments....
	double diff = (-v_opt_dishes[opt][current.x_var[neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][neighbor.newValue].v_nutrient_value[j]);
        current.nutriment_per_day[a][day][j] += diff;
	current.globalPlan[a][j] += diff;
    }
   }
   current.x_var[neighbor.variable]= neighbor.newValue;
   current.obj_values = new_objs;
   current.uniq_per_day[day][dish_in.description]++, current.uniq_per_day[day][dish_out.description]--;
}
void MPP::init_incremental_evaluation(struct Solution_LS &current)
{
        //feasibility information
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
        vector<int> &v_constraint_global = MPP_problem->v_constraint_global, &v_constraint_day = MPP_problem->v_constraint_day;
        int &max_description_id = MPP_problem->max_description_id;
        vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   	current.uniq_per_day.assign(nDias+1, vector<int> (MPP_problem->max_description_id+1, 0));
        current.globalPlan.assign((int)MPP_problem->conf_day.size(),vector<double> (num_nutr, 0.0 ));
	current.nutriment_per_day.assign( (int)MPP_problem->conf_day.size(), vector<vector<double> > (nDias, vector<double> (num_nutr, 0)));
	current.obj_values.assign(N_TIMES+2, 0.0);
        for(int a = 0; a < MPP_problem->conf_day.size(); a++)
	{
 	  for(int j = 0; j < num_nutr; j++)
	  {
	   for(int i = 0; i < nDias; i++)
	   {
		 for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
		 {
		   int k = MPP_problem->conf_day[a][b];
	    	   current.nutriment_per_day[a][i][j] += v_opt_dishes[k][current.x_var[i*N_OPT_DAY + k]].v_nutrient_value[j];
		 }
	        current.globalPlan[a][j] += current.nutriment_per_day[a][i][j];
	        if( v_constraints[j].type == DIARIA)
                {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = current.nutriment_per_day[a][i][j];
	           if( nut < minv) current.obj_values[0] += ((minv - nut)/middle)*((minv - nut)/middle)*WEIGHT_DAY;
	           else if (nut > maxv) current.obj_values[0] +=((nut - maxv)/middle)*((nut - maxv)/middle)*WEIGHT_DAY;
                }
	   }
	   if( v_constraints[j].type == GLOBAL)
           {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = current.globalPlan[a][j];
	           if( nut < minv) current.obj_values[1] += ((minv - nut)/middle)*((minv - nut)/middle);
	           else if (nut > maxv) current.obj_values[1] += ((nut - maxv)/middle)*((nut - maxv)/middle);
            }
	  }
	}
        calculateVariability(current.x_var, current.obj_values);
}
void MPP::calculateFeasibilityDegree(vector<int> &sol, double &feas_day, double &feas_global){
	feas_day = 0.0;
	feas_global = 0.0;
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
	badDaysFeas.clear();
	double heaviestValue = 0;
	heaviestNut = -1;
        for(int a = 0; a < MPP_problem->conf_day.size(); a++)
        {
           for(int i = 0; i < num_nutr; i++)
           {
              double minv = v_constraints[i].min;
              double maxv = v_constraints[i].max;
	      double middle = (maxv+minv)*0.5;
              double globalNutr = 0.0;
	      for(int j = 0; j < nDias; j++)
	      {
	         double dayNutr = 0.0;
		 for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
		 {
		   int k = MPP_problem->conf_day[a][b];
	   	   dayNutr += MPP_problem->v_opt_dishes[k][sol[j*N_OPT_DAY+k]].v_nutrient_value[i];
		 }
	   	globalNutr += dayNutr;
	   	if(v_constraints[i].type == DIARIA)
             	{
	               if(dayNutr  < minv) feas_day += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle)*WEIGHT_DAY, badDaysFeas.insert(j);
	               else if (dayNutr > maxv) feas_day +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle)*WEIGHT_DAY, badDaysFeas.insert(j);
           	}
	      }
	      if(v_constraints[i].type == GLOBAL)
              {
	        if(globalNutr < minv)
	        {
	          double v = ((minv - globalNutr)/middle)*((minv - globalNutr)/middle);
	          feas_global += v;
	          if( v >  heaviestValue)
	          {
	            heaviestValue = v;
	            heaviestNut = i;
	            heaviestType = -1;
	          }
	        }
	        else if (globalNutr > maxv)
	        {
	          double v =((globalNutr - maxv)/middle)*((globalNutr - maxv)/middle);
	          feas_global +=v;
	          if( v >  heaviestValue)
	          {
	            heaviestValue = v;
	            heaviestNut = i;
	            heaviestType = 1;
	          }
	        }
              }
           }
	}
}
void MPP::feasibility_day(vector<int> &best_solution, double &best_feasibility)
{
   vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   //get daily contribution
   best_feasibility = 0.0;
   for(int a = 0; a < MPP_problem->conf_day.size(); a++)
   {
   for(int ij = 0; ij < v_constraint_day.size(); ij++)
     {
        double dayNutr = 0;
        int i =  v_constraint_day[ij];
        double minv = v_constraints[i].min;
        double maxv = v_constraints[i].max;
	double middle = (maxv+minv)*0.5;
	for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
	{
	  int k = MPP_problem->conf_day[a][b];
	  dayNutr += MPP_problem->v_opt_dishes[k][best_solution[k]].v_nutrient_value[i];
        }
	if(dayNutr  < minv) best_feasibility += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle);
	else if (dayNutr > maxv) best_feasibility +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle);
     }
   }
}
void MPP::init_inc_eval_day(vector<int> &current_solution, vector<vector<double>> &nut_info)
{
   int num_nutr = (int)MPP_problem->v_constraints.size();
   vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;

   nut_info.assign((int)MPP_problem->conf_day.size(), vector<double>(num_nutr, 0.0));
   for(int a = 0; a < MPP_problem->conf_day.size(); a++)
   {
     for(int i = 0; i < v_constraint_day.size(); i++)
     {
       int j =  v_constraint_day[i];
       for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
       {
         int k = MPP_problem->conf_day[a][b];
         nut_info[a][j] += v_opt_dishes[k][current_solution[k]].v_nutrient_value[j];
       }
     }
   }

}
double MPP::inc_eval_day(Neighbor &new_neighbor, vector<vector<double>> &nut_info, vector<int> &current_solution, double current_feasibility)
{
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   double new_partial_infeasibility_day = 0.0, original_partial_infeasibility_day = 0.0 ;
   int opt = new_neighbor.variable;
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {
       int a = MPP_problem->opt_conf[opt][b];
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
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   int opt = new_neighbor.variable;
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {
       int a = MPP_problem->opt_conf[opt][b];
       for(int i = 0; i < v_constraint_day.size(); i++)
       {
          int j =  v_constraint_day[i];
         nut_info[a][j] = nut_info[a][j]-v_opt_dishes[opt][current_solution[new_neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][new_neighbor.newValue].v_nutrient_value[j];
       }
   }
   current_solution[opt] = new_neighbor.newValue;
}
void MPP::calculateVariability(vector<int> &sol, vector<double> &objs_var)
{
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   double variability_global = 0.0, variability_cat_day=0.0;
   int &max_description_id = MPP_problem->max_description_id;
   double v_global = 0, v_global_id = 0, v_global_cat = 0;
   badDaysVar.clear();
   vector< vector<int> > last_day_seen(N_TIMES, vector<int>(max_description_id+1, -1)), last_day_seen_cat(N_TIMES, vector<int>(N_TIMES+1, -1));
   vector<pair<int, int>> min_dcn(N_TIMES, make_pair(nDias+1, 0)), min_dcn_cat(N_TIMES, make_pair(N_CATEGORIES+1, 0));
   int min_dcn_day = nDias+1;
   for(int d = 0; d < nDias; d++)
   {
	for(int i = 0; i < MPP_problem->unique_opt_time.size(); i++)
	{
	   for(int j = 0; j < MPP_problem->unique_opt_time[i].size(); j++)
           {
	  	int opt = MPP_problem->unique_opt_time[i][j];
		struct infoDishes &dish = v_opt_dishes[opt][sol[d*N_OPT_DAY + opt]];
           	if(last_day_seen[i][dish.description] != -1)
           	{
           	   int diff = d-last_day_seen[i][dish.description];
	   	   diff = min((dish.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE, diff);
	   	   update_dcn_pair(diff, min_dcn[i]);
		   if( min_dcn_day == diff)
		   {
		     badDaysVar.insert(d);
		     badDaysVar.insert(last_day_seen[i][dish.description]);
		   }
		   else if(min_dcn_day > diff)
		   {
		     min_dcn_day = diff;
		     badDaysVar.clear();
		     badDaysVar.insert(d);
		     badDaysVar.insert(last_day_seen[i][dish.description]);
		   }
           	}
           	last_day_seen[i][dish.description] = d;
	 	if(dish.category == CATEGORY_BOTH)
		{
		   for(int c = 1; c <=2; c++)
		   {
	   	      if(last_day_seen_cat[i][c] != -1)
           	      {
           	         int diff = d-last_day_seen_cat[i][c];
	   	         update_dcn_pair(diff, min_dcn_cat[i]);
           	      }
           	      last_day_seen_cat[i][c] = d;
		   }
		}
		else
		{
		   if(last_day_seen_cat[i][dish.category] != -1)
           	   {
           	      int diff = d-last_day_seen_cat[i][dish.category];
	   	      update_dcn_pair(diff, min_dcn_cat[i]);
           	   }
           	   last_day_seen_cat[i][dish.category] = d;
		}
	   }
	}
   }
  for(int i = 0; i < N_TIMES; i++)
  {
    objs_var[i+2] = W_VAR_GLOBAL*f(min_dcn[i]) + W_VAR_GLOBAL_CAT*f(min_dcn_cat[i]);
  }
}
