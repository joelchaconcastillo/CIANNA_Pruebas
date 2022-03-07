/*
 * Developers: Joel Chacón Castillo, Carlos Segura
 * Last modification: 03/05/22
 * This file defines the model that is taken into account as the fitness of each individuals given a combination of dishes
 * */
#include "global.h"
#include "MPP.h"
extern vector<vector<int> > g_timesIdPerConf;
extern vector<vector<int>> g_time2Ids;
extern vector<double> g_weights;

void MPP::evaluate(){
  calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
  calculateVariability(x_var, obj_values);
  fitness = obj_values[0]+obj_values[1];
  double ref = nDias*W_VAR_GLOBAL;
  double max1 = -ref, sum1=0.0;
  for(int a = 2; a < obj_values.size(); a++)
  {
     double v1 = obj_values[a];
	sum1 += fabs(v1-ref);
	max1 = max(max1,fabs(v1-ref)/g_weights[a-2]);
  }
  fitness = 1e4*(obj_values[0]+obj_values[1]) + max1;
}
void MPP::feasibility_day(vector<int> &best_solution, double &best_feasibility){
   //get daily contribution
   best_feasibility = 0.0;
   for(int a = 0; a < g_timesIdPerConf.size(); a++)
   {
   for(int ij = 0; ij < v_constraint_day.size(); ij++)
     {
        double dayNutr = 0;
        int i =  v_constraint_day[ij];
        double minv = v_constraints[i].min;
        double maxv = v_constraints[i].max;
	double middle = (maxv+minv)*0.5;
	for(int b = 0; b < g_timesIdPerConf[a].size(); b++)
	{
	  int k = g_timesIdPerConf[a][b];
	  dayNutr += v_opt_dishes[k][best_solution[k]].v_nutrient_value[i];
        }
	if(dayNutr  < minv) best_feasibility += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle);
	else if (dayNutr > maxv) best_feasibility +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle);
     }
   }
}
void MPP::calculateFeasibilityDegree(vector<int> &sol, double &feas_day, double &feas_global){
	feas_day = 0.0;
	feas_global = 0.0;
	int num_nutr = (int)v_constraints.size();
	badDaysFeas.clear();
	double heaviestValue = 0;
	heaviestNut = -1;
        for(int a = 0; a < g_timesIdPerConf.size(); a++)
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
		 for(int b = 0; b < g_timesIdPerConf[a].size(); b++)
		 {
		   int k = g_timesIdPerConf[a][b];
	   	   dayNutr += v_opt_dishes[k][sol[j*N_OPT_DAY+k]].v_nutrient_value[i];
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
void MPP::calculateVariability(vector<int> &sol, vector<double> &objs_var)
{
   double variability_global = 0.0, variability_cat_day=0.0;
   double v_global = 0, v_global_id = 0, v_global_cat = 0;
   badDaysVar.clear();
   vector< vector<int> > last_day_seen(N_TIMES, vector<int>(max_description_id+1, -1)), last_day_seen_cat(N_TIMES, vector<int>(N_TIMES+1, -1));
   vector<pair<int, int>> min_dcn(N_TIMES, make_pair(nDias+1, 0)), min_dcn_cat(N_TIMES, make_pair(N_CATEGORIES+1, 0));
   int min_dcn_day = nDias+1;
   for(int d = 0; d < nDias; d++){
	for(int i = 0; i < g_time2Ids.size(); i++)
	{
	   for(int j = 0; j < g_time2Ids[i].size(); j++)
           {
	  	int opt = g_time2Ids[i][j];
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
  for(int i = 0; i < N_TIMES; i++){
    objs_var[i+2] = W_VAR_GLOBAL*f(min_dcn[i]) + W_VAR_GLOBAL_CAT*f(min_dcn_cat[i]);
  }
}
