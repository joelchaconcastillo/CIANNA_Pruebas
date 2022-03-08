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
extern vector<double> g_rank;
extern vector<double> g_rank_category;

double MPP::modelation(vector<double> &objFeas, vector<pair<double, double>> &objVar){
  double score = 0;
  double sum=0;
  for(int i = 0; i < N_TIMES; i++){
	  //////dcn by foods
	  double dcni=objVar[i].first;//dcn score max
	  double repi=1.0-(objVar[i].second/((double) nDias)); //repetitions of this dcn score max
	  double normalized = (1.0-((dcni+repi)/((double) nDias)))+g_rank[i];//min..
	  //dcn by categories..
	  double dcnCat=objVar[i+N_TIMES].first;
	  double repCat=(1.0-objVar[i+N_TIMES].second/((double)nDias));
	  double normalizedCat = (1.0-((dcnCat+repCat)/((double) nDias)))+g_rank_category[i];
	  sum += normalized+normalizedCat;
  }
  score=W_FEASIBILITY*(objFeas[0]+objFeas[1]) + sum;
  return score;
}
void MPP::evaluate(vector<int> &sol, vector<double> &objFeas, vector<pair<double, double> > &objVar){
  calculateFeasibilityDegree(sol, objFeas[0], objFeas[1]);
  calculateVariability(sol, objVar);
  fitness = modelation(objFeas, objVar);
}
/**
 *This function defines the mix the problem's information 
 *into a fitness function
 * */
void MPP::evaluate(){
  evaluate(x_var, objFeasibility, objVariability);
}
bool MPP::comp_objs(vector<double> &objFeas1, vector<pair<double, double> > &objVar1, vector<double> &objFeas2, vector<pair<double, double> > &objVar2){
  //feasibility
  if( objFeas1[0] < objFeas2[0]) return true;
  else if( objFeas1[0] > objFeas2[0]) return false;

  if( objFeas1[1] < objFeas2[1]) return true;
  else if( objFeas1[1] >  objFeas2[1]) return false;

  double score1=modelation(objFeas1, objVar1), score2=modelation(objFeas2, objVar2);
  return (score1<score2);
}
/*
 *In this functions is calculated the feasibility per day
 * */
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
	//reset variables
	feas_day = 0.0;
	feas_global = 0.0;
	badDaysFeas.clear();
	double heaviestValue = 0;
	heaviestNut = -1;
	for(auto &id_config:g_timesIdPerConf){
           for(int i = 0; i < num_nutr; i++){
              double minv = v_constraints[i].min;
              double maxv = v_constraints[i].max;
	      double middle = (maxv+minv)*0.5;
              double globalNutr = 0.0;
	      for(int j = 0; j < nDias; j++){
	         double dayNutr = 0.0;
		 for(auto &id_time:id_config){
	   	   dayNutr += v_opt_dishes[id_time][sol[j*N_OPT_DAY+id_time]].v_nutrient_value[i];
		 }
	   	globalNutr += dayNutr;
	   	if(v_constraints[i].type == DIARIA){
	               if(dayNutr  < minv) feas_day += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle)*WEIGHT_DAY, badDaysFeas.insert(j);
	               else if (dayNutr > maxv) feas_day +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle)*WEIGHT_DAY, badDaysFeas.insert(j);
           	}
	      }
	      if(v_constraints[i].type == GLOBAL){
	        if(globalNutr < minv){
	          double v = ((minv - globalNutr)/middle)*((minv - globalNutr)/middle);
	          feas_global += v;
	          if( v >  heaviestValue){
	            heaviestValue = v;
	            heaviestNut = i;
	            heaviestType = -1;
	          }
	        }
	        else if (globalNutr > maxv){
	          double v =((globalNutr - maxv)/middle)*((globalNutr - maxv)/middle);
	          feas_global +=v;
	          if( v >  heaviestValue){
	            heaviestValue = v;
	            heaviestNut = i;
	            heaviestType = 1;
	          }
	        }
              }
           }
	}
}
void MPP::calculateVariability(vector<int> &sol, vector<pair<double, double> > &objVar){
   badDaysVar.clear();
   vector< vector<int> > last_day_seen(N_TIMES, vector<int>(max_description_id+1, -1)), last_day_seen_cat(N_TIMES, vector<int>(N_TIMES+1, -1));
   vector<pair<int, int>> min_dcn(N_TIMES, make_pair(nDias+1, 0)), min_dcn_cat(N_TIMES, make_pair(N_CATEGORIES+1, 0));
   int min_dcn_day = nDias+1;
   for(int d = 0; d < nDias; d++){
	for(int i = 0; i < g_time2Ids.size(); i++){
	   for(auto &opt:g_time2Ids[i]){
		struct infoDishes &dish = v_opt_dishes[opt][sol[d*N_OPT_DAY + opt]];
           	if(last_day_seen[i][dish.description] != -1){
           	   int diff = d-last_day_seen[i][dish.description];
		   //just penalize distances below a threshold
		   //in this way repetition are possible
		   //for instance favorite days can be repeated every "DAYS_FAVORITE"
	   	   diff = min((dish.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE, diff);
                   diff=min(diff, 3);
	   	   update_dcn_pair(diff, min_dcn[i]);
		   if( min_dcn_day == diff){
		     badDaysVar.insert(d);
		     badDaysVar.insert(last_day_seen[i][dish.description]);
		   }
		   else if(min_dcn_day > diff){
		     min_dcn_day = diff;
		     badDaysVar.clear();
		     badDaysVar.insert(d);
		     badDaysVar.insert(last_day_seen[i][dish.description]);
		   }
           	}
           	last_day_seen[i][dish.description] = d;
	 	if(dish.category == CATEGORY_BOTH){
		   for(int c = 1; c <=2; c++){
	   	      if(last_day_seen_cat[i][c] != -1){
           	         int diff = d-last_day_seen_cat[i][c];
	   	         update_dcn_pair(diff, min_dcn_cat[i]);
           	      }
           	      last_day_seen_cat[i][c] = d;
		   }
		}
		else{
		   if(last_day_seen_cat[i][dish.category] != -1){
           	      int diff = d-last_day_seen_cat[i][dish.category];
	   	      update_dcn_pair(diff, min_dcn_cat[i]);
           	   }
           	   last_day_seen_cat[i][dish.category] = d;
		}
	   }
	}
  }
  for(int i = 0; i < N_TIMES; i++){
	  objVar[i]= min_dcn[i];
	  objVar[i+N_TIMES]=min_dcn_cat[i];
  }
}
