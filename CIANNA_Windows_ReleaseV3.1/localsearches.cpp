#include <chrono>
#include <sys/time.h>
#include "global.h"
#include "MPP.h"
using namespace std;
using namespace std::chrono;
extern vector<vector<int>> g_Idtime2Configs;
extern vector<vector<int> > g_timesIdPerConf;
extern vector<vector<int>> g_time2Ids;
extern int g_crossoverType;
extern vector<double> g_l_rank;
extern vector<double> g_u_rank;
extern vector<double> g_l_rank_category;
extern vector<double> g_u_rank_category;


/*
 *This is the main function that integrates two local searches:
   First_Improvement_Hill_Climbing: Tries convergency into a feasible solution and as second objective it tries to improve in terms of variability.
   First_Improvement_Hill_Climbing_swap: Aims to maximize the variability between dishes
 * */
double MPP::localSearch(double finalTime, bool isMeasuringTime){
     calculateFeasibilityDegree(x_var, objFeasibility[0], objFeasibility[1]);
     calculateVariability(x_var, objVariability);
     vector<int> bestIndividual = x_var;
     vector<double> bestFeas = objFeasibility;
     vector<pair<double, double> > bestVar=objVariability;
     //load incremental evaluation values...
     struct timeval currentTime;
     gettimeofday(&currentTime, NULL);
     double initialTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
     double elapsedTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
     elapsedTime -= initialTime;
     while(elapsedTime < finalTime){
        First_Improvement_Hill_Climbing(neighbors, x_var, objFeasibility, objVariability);
        First_Improvement_Hill_Climbing_swap(neighbors_swap, x_var, objFeasibility, objVariability);
        if(comp_objs(objFeasibility, objVariability, bestFeas, bestVar)){
	  bestFeas=objFeasibility;
	  bestVar=objVariability;
          bestIndividual = x_var;
//	  cout <<setprecision(10)<< modelation(bestFeas, bestVar)<<endl;
        }
        else{
	       	x_var = bestIndividual;
	}

        if(isMeasuringTime && bestFeas[0] == 0.0 && bestFeas[1] == 0.0){
	      	break;
	}
        //These two functions are necessary to fix numerical errors raised by the incremental evaluations..
        calculateFeasibilityDegree(x_var, objFeasibility[0], objFeasibility[1]);
        calculateVariability(x_var, objVariability);

        if(badDaysFeas.empty()){
	   int selectedDay = -1;
	   if (heaviestNut != -1){
		vector< pair<double, int> > infoNut;
		for (int i = 0; i < nDias; i++){
		   double total = 0.0;
		   for(int k = 0; k < N_OPT_DAY; k++) total +=  v_opt_dishes[k][x_var[i*N_OPT_DAY + k]].v_nutrient_value[heaviestNut];
		   infoNut.push_back(make_pair(total, i));
		}
		sort(infoNut.begin(), infoNut.end());
		if (heaviestType == 1){
		       	reverse(infoNut.begin(), infoNut.end());
		}
		int nbestday = rand()%min(nDias ,  6);
		selectedDay = infoNut[nbestday].second;
		perturb_day(x_var, selectedDay);
	        oneDaylocalSearch(x_var, selectedDay);
	   }
	   else{
		 selectedDay = rand() % nDias;
		 perturb_day(x_var, selectedDay);
	         oneDaylocalSearch(x_var, selectedDay);
	   }
	 }
	 else{
	   for (auto it = badDaysFeas.begin(); it != badDaysFeas.end(); it++){
	      perturb_day(x_var, *it);
	      oneDaylocalSearch(x_var, *it);
	   }
 	 }
	gettimeofday(&currentTime, NULL);
	elapsedTime = ((double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6)-initialTime;
     }

     x_var = bestIndividual;
     evaluate();
     return elapsedTime;
}

/*
 *In this function is applied the first improvement hill climbing by randomly swapping days,
 *this is better to increase variability between dishes, although swapping can be integrated into crossover, 
 *some preliminary test indicates that to converge faste is better a local search instead into the crossover.
 * */
void MPP::First_Improvement_Hill_Climbing_swap(vector<Neighbor_swap> &neighbors, vector<int> &best_sol, vector<double> &bestFeas, vector<pair<double, double> > &bestVar){
  evaluate(best_sol, bestFeas, bestVar);
  bool improved= true;
  vector<int> current_sol = best_sol;
  vector<pair<double, double> > currentVar = bestVar;
  vector<double> currentFeas=bestFeas;
  while(improved){
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for(int i = 0; i < neighbors.size(); i++){
	swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
	calculateVariability(current_sol, currentVar);
	if( comp_objs(currentFeas, currentVar, bestFeas, bestVar)){
            improved = true;
	    bestFeas=currentFeas;
	    bestVar=currentVar;
	    swap_days(best_sol, neighbors[i].day1, neighbors[i].day2);
	}
	else
	    swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
     }
  }
  evaluate(best_sol, bestFeas, bestVar);
}
/*
 *This is the "First improvement hill climbing regarding first feasibility and second to the variability
 * */
void  MPP::First_Improvement_Hill_Climbing(vector<Neighbor> &neighbors, vector<int> &currentSol, vector<double> &currentFeas, vector<pair<double, double> > &currentVar){

   evaluate(currentSol, currentFeas, currentVar); 
   //incremental evaluation values...
   struct Solution_LS best;
   best.objFeasibility=currentFeas;
   best.objVariability=currentVar;
   best.x_var = currentSol;
   vector<double> newFeas = currentFeas;
   vector<pair<double, double> > newVar=currentVar;
   init_incremental_evaluation(best);
   bool improved = true;
   while(improved){
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for (int i = 0; i < neighbors.size(); i++){
	int day = neighbors[i].variable/N_OPT_DAY, opt =  neighbors[i].variable%N_OPT_DAY;
        struct infoDishes &dish_in = v_opt_dishes[opt][neighbors[i].newValue], &dish_out = v_opt_dishes[opt][best.x_var[day*N_OPT_DAY + opt]];
	//force to always have an unique dish per day..
	if(best.uniq_per_day[day][dish_in.description]>0) continue;
        //incremental evaluation...
	inc_eval(best, neighbors[i], newFeas, newVar);
        if(comp_objs(newFeas, newVar, best.objFeasibility, best.objVariability)){
	   improved = true;
	   update_inc(best, neighbors[i], newFeas, newVar);
	}
      }
    }
    currentSol= best.x_var;
    currentFeas = best.objFeasibility;
    currentVar = best.objVariability;
    evaluate(currentSol, currentFeas, currentVar); //check again.. the overall values....//avoid numerical error of the incremental evaluation sum..
}
void MPP::oneDaylocalSearch(vector<int> &solution, int day){
     vector<int> best_solution(solution.begin()+day*N_OPT_DAY, solution.begin() + (day+1)*N_OPT_DAY );
     double best_feasibility = 0.0;
     feasibility_day(best_solution, best_feasibility);
     double current_feasibility = best_feasibility;
     vector<int> current_solution = best_solution;
     auto start = high_resolution_clock::now();
     for (int i = 0; i < MAX_ITE_ONE_DAY_LS; i++){
 	First_Improvement_Hill_Climbing_Day(oneDayneighbors, current_solution, current_feasibility);
 	if( current_feasibility < best_feasibility){
	   best_feasibility = current_feasibility;
	   best_solution = current_solution;
	}
	else{
	   current_solution = best_solution;
	   current_feasibility = best_feasibility;
	}
	//restart....
	int which = rand() % N_OPT_DAY;
	perturb_opt(current_solution, 0, which);
	feasibility_day(current_solution, current_feasibility);
	if(current_feasibility == 0.0) break;
     }
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);
    for(int i = 0; i  < N_OPT_DAY; i++) solution[day*N_OPT_DAY + i] = best_solution[i];
}
void  MPP::First_Improvement_Hill_Climbing_Day(vector<Neighbor> &oneDayneighbors, vector<int> &best_solution, double &best_feasibility){
   //incremental evaluation values...
   double current_feasibility = best_feasibility;
   vector<bool> uniq_per_day(max_description_id, false);

   vector<vector<double>> nut_info;
   for(int i = 0 ; i  < N_OPT_DAY; i++) uniq_per_day[v_opt_dishes[i][best_solution[i]].description] = true;
   init_inc_eval_day(best_solution, nut_info);
   bool improved = true;
   while(improved)
   {
     improved = false;
     random_shuffle(oneDayneighbors.begin(), oneDayneighbors.end());
     for (int i = 0; i < oneDayneighbors.size(); i++)
     {
	int day = oneDayneighbors[i].variable/N_OPT_DAY, opt =  oneDayneighbors[i].variable;
        struct infoDishes &dish_in = v_opt_dishes[opt][oneDayneighbors[i].newValue];
        struct infoDishes &dish_out = v_opt_dishes[opt][best_solution[opt]];
	//this is to avoid repeated dishes per day
	if(uniq_per_day[dish_in.description]) continue;
        //incremental evaluation...
	double new_feasibility = inc_eval_day(oneDayneighbors[i], nut_info , best_solution, best_feasibility);
	if( new_feasibility < best_feasibility)
	{
	   improved = true;
	   update_inc_day(nut_info, oneDayneighbors[i], best_solution);
	   best_feasibility = new_feasibility;
	   best_feasibility = max(0.0, best_feasibility);
	   uniq_per_day[dish_in.description] = true;
	   uniq_per_day[dish_out.description] = false;
	}
      }
    }
}
/**
 *This section is devoted to the incremental evaluation regarding the feasibility
 *therefore some parts of the sum are not longer recomputed.
 * */
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
/***
 *
 *This section is dedicated to the incremental evaluation per day
 *this coding optimization gives a boosting to the time!
 *
 * */
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
double MPP::inc_eval_day(Neighbor &new_neighbor, vector<vector<double>> &nut_info, vector<int> &current_solution, double current_feasibility){
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
void MPP::update_inc_day(vector<vector<double>> &nut_info, Neighbor &new_neighbor, vector<int> &current_solution){
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
void MPP::Greedy_Improvement(vector<int> &x_var, vector<double> &objFeasibility, vector<pair<double, double> > &objVariability){
   vector<int> new_x_var((int)x_var.size(), 0);
   auto newFeas=objFeasibility;
   auto newVar=objVariability;

   //find the heaviest time
   pair<long long, int> min_time(-INT_MAX, 0);
   for(int i = 0; i < g_u_rank.size(); i++) min_time=max(min_time, {g_u_rank[i], i});
   int idHeaviest=min_time.second;
   map<int, vector<int> > itemDay;
   map<int, int> rep;
   int maxCount = -INT_MAX;
   for(int i = 0; i <nDias; i++){
     int item = x_var[i*N_OPT_DAY+idHeaviest];
     itemDay[item].push_back(i);
     rep[item]++;
     maxCount=max(maxCount, rep[item]);
   }
   vector< pair<int, int> > items;
   for(auto ii:rep) items.push_back({-ii.second, ii.first});
   sort(items.begin(), items.end());
   vector<vector<int> > buckets(maxCount, vector<int>());
   for(int i = 0, k=0; i < items.size(); i++){
      for(int j = 0; j < (int)-items[i].first; j++){
	  buckets[k].push_back(items[i].second);
	  k = (k+1)%maxCount;
      }
   }
   for(int i = 0, toDay=0; i < buckets.size(); i++){
      for(int j = 0; j < buckets[i].size(); j++){
	  int fromDay = itemDay[buckets[i][j]].back();
	  itemDay[buckets[i][j]].pop_back();
	  for(int k = 0; k < N_OPT_DAY; k++){
	     new_x_var[toDay*N_OPT_DAY+k]=x_var[fromDay*N_OPT_DAY+k];
	  }
	  toDay++;
      }
   }
   evaluate(new_x_var, newFeas, newVar);
   auto newfit = modelation(newFeas, newVar);
   fitness=modelation(objFeasibility, objVariability);
   if(newfit < fitness){
     x_var=new_x_var; 
     objFeasibility=newFeas;
     objVariability=newVar;
     fitness=newfit;
   }
}
void MPP::First_Improvement_Hill_Climbing_swap2(vector<Neighbor_swap> &neighbors, vector<int> &best_sol, vector<double> &bestFeas, vector<pair<double, double> > &bestVar){

  bool improved= true;
  vector<int> current_sol = best_sol;
  vector<pair<double, double> > currentVar = bestVar;
  vector<double> currentFeas=bestFeas;
  while(improved){
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for(int i = 0; i < neighbors.size(); i++){
        swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
        swap_days(current_sol, neighbors[i].day3, neighbors[i].day4);
        calculateVariability(current_sol, currentVar);
        if( comp_objs(currentFeas, currentVar, bestFeas, bestVar)){
            improved = true;
//	    cout <<setprecision(10) <<"MEJORA! "<<modelation(bestFeas, bestVar) << " ";
            bestFeas=currentFeas;
            bestVar=currentVar;
//	    cout <<setprecision(10) <<" "<<modelation(bestFeas, bestVar) << endl;
            swap_days(best_sol, neighbors[i].day1, neighbors[i].day2);
            swap_days(best_sol, neighbors[i].day3, neighbors[i].day4);
        }
        else{
           swap_days(current_sol, neighbors[i].day3, neighbors[i].day4);
           swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
	}
     }
  }
}
void MPP::Improvement(vector<int> &best_sol, vector<double> &bestFeas, vector<pair<double, double> > &bestVar){
   int K=3;
   calculateVariability(best_sol, bestVar);
   auto candidateX=best_sol;
   auto candidateVar=bestVar;
   auto xb=best_sol;
   double fitb=modelation(bestFeas, bestVar);
   vector<int> permutationDays;
   bool improved=true;
   for(int i=0;i<nDias; i++)permutationDays.push_back(i);
//   while(improved)
   {
	   improved=false;
   random_shuffle(permutationDays.begin(), permutationDays.end());
      for(auto idDay:permutationDays){
        for(int i = max(0, idDay-K); i < min(idDay+K, nDias); i++){
      	  vector<int> perm;
              for(int ii = i; ii < min(i+K, nDias); ii++)perm.push_back(ii);
   	   do{
   	      for(int j=0; j < perm.size(); j++){
   	        for(int jj = 0; jj < N_OPT_DAY; jj++){
   		   candidateX[(i+j)*N_OPT_DAY+jj]=best_sol[perm[j]*N_OPT_DAY+jj];
   	        }
   	      }
   	      calculateVariability(candidateX, candidateVar);
   	      auto nf=modelation(bestFeas,candidateVar);
   	      if( nf<fitb){
   		     fitb=nf;
   		     xb=candidateX;
   		     improved=true;
   		     //cout << nf<<endl;
   	      } 
   	   }while(next_permutation(perm.begin(), perm.end()));
           best_sol=xb;
           candidateX=xb;
           calculateVariability(best_sol, bestVar);
           candidateVar=bestVar;
           fitb=modelation(bestFeas, bestVar);
	}
      }
    }
}

