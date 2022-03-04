#include <chrono>
#include <sys/time.h>
#include "MPP.h"
using namespace std;
using namespace std::chrono;

void MPP::localSearch(double finalTime){

     vector<Neighbor> neighbors;
     for (int i = 0; i < nDias; i++)
     {
	for (int j = 0; j <N_OPT_DAY; j++)
	{
	   for (int k = 0; k < (int) MPP_problem->v_opt_dishes[j].size(); k++)
	   {
	     Neighbor n;
	     n.variable = i * N_OPT_DAY + j;
	     n.newValue = k;
	     neighbors.push_back(n);
	   }
        }
     }
     vector<Neighbor_swap> neighbors_swap;
     for(int i = 0; i < nDias; i++)
     {
      for(int j = i+1; j < nDias; j++)
      {
         neighbors_swap.push_back({i, j});
      }
     }

     calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
     calculateVariability(x_var, obj_values);

     vector<int> bestIndividual = x_var;
     vector<double> best_objs = obj_values;
     //load incremental evaluation values...
     struct timeval currentTime;
     gettimeofday(&currentTime, NULL);
     double initialTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
     double elapsedTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
     elapsedTime -= initialTime;
     while(elapsedTime < finalTime )
     {
        calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
        calculateVariability(x_var, obj_values);
        First_Improvement_Hill_Climbing(neighbors, x_var, obj_values);
        First_Improvement_Hill_Climbing_swap(neighbors_swap, x_var, obj_values);
        if(comp_objs(obj_values, best_objs))
        {
	  best_objs = obj_values;
          bestIndividual = x_var;
//	cout << best_objs<<endl;
        }
        else x_var = bestIndividual;

        calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
        calculateVariability(x_var, obj_values);
        if(badDaysFeas.empty())
        {
	   int selectedDay = -1;
	   if (heaviestNut != -1)
           {
		vector< pair<double, int> > infoNut;
		for (int i = 0; i < nDias; i++)
		{
		   double total = 0.0;
		   for(int k = 0; k < N_OPT_DAY; k++) total +=  MPP_problem->v_opt_dishes[k][x_var[i*N_OPT_DAY + k]].v_nutrient_value[heaviestNut];
		   infoNut.push_back(make_pair(total, i));
		}
		sort(infoNut.begin(), infoNut.end());
		if (heaviestType == 1) reverse(infoNut.begin(), infoNut.end());
		int nbestday = rand()%min(nDias ,  6);
		selectedDay = infoNut[nbestday].second;
	   }
	   else
	   {
		if( rand()%2 && !badDaysVar.empty())
		{
	 	   vector<int> tmp;
	           for(auto it = badDaysVar.begin(); it != badDaysVar.end(); it++) tmp.push_back(*it);
		   random_shuffle(tmp.begin(), tmp.end());
		   selectedDay = tmp[0];
		}
		else
		 selectedDay = rand() % nDias;
	   }
	   perturb_day(x_var, selectedDay);
	   oneDaylocalSearch(x_var, selectedDay);
	 }
	 else
	 {
	   for (auto it = badDaysFeas.begin(); it != badDaysFeas.end(); it++)
	   {
	      perturb_day(x_var, *it);
	      oneDaylocalSearch(x_var, *it);
	   }
 	 }
	gettimeofday(&currentTime, NULL);
	elapsedTime = ((double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6)-initialTime;
     }
     x_var = bestIndividual;
//	cout << obj_values<<endl;
     evaluate();
}
void MPP::First_Improvement_Hill_Climbing_swap(vector<Neighbor_swap> &neighbors, vector<int> &best_sol, vector<double> &best_objs)
{
  bool improved= true;
  vector<int> current_sol = best_sol;
  vector<double> current_objs = best_objs;
  while(improved)
  {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for(int i = 0; i < neighbors.size(); i++)
     {
	swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
	calculateVariability(current_sol, current_objs);
	if( comp_objs(current_objs, best_objs))
	{
            improved = true;
	    best_objs= current_objs;
	    swap_days(best_sol, neighbors[i].day1, neighbors[i].day2);
	}
	else
	    swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
     }
  }
}
void  MPP::First_Improvement_Hill_Climbing(vector<Neighbor> &neighbors, vector<int> &best_sol, vector<double> &best_objs)
{
   //incremental evaluation values...
   struct Solution_LS best;
   best.obj_values = best_objs;
   best.x_var = best_sol;
   vector<double> new_objs = best_objs;
   init_incremental_evaluation(best);
   bool improved = true;
   while(improved)
   {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for (int i = 0; i < neighbors.size(); i++)
     {
	int day = neighbors[i].variable/N_OPT_DAY, opt =  neighbors[i].variable%N_OPT_DAY;
        struct infoDishes &dish_in = MPP_problem->v_opt_dishes[opt][neighbors[i].newValue], &dish_out = MPP_problem->v_opt_dishes[opt][best.x_var[day*N_OPT_DAY + opt]];
	if(best.uniq_per_day[day][dish_in.description]>0) continue;
        //incremental evaluation...
	inc_eval(best, neighbors[i], new_objs);
        if(comp_objs(new_objs, best.obj_values))
	{
	   improved = true;
	   update_inc(best, neighbors[i], new_objs);
	}
      }
    }
    best_sol = best.x_var;
    best_objs = best.obj_values;
    evaluate(best_sol, best_objs); //check again.. the overall values....//avoid numerical error of the incremental evaluation sum..
}
void MPP::oneDaylocalSearch(vector<int> &solution, int day) {
     vector<Neighbor> neighbors;
     int num_nutr = (int)MPP_problem->v_constraints.size();
     vector<int> best_solution(solution.begin()+day*N_OPT_DAY, solution.begin() + (day+1)*N_OPT_DAY );

     for (int j = 0; j <N_OPT_DAY; j++)
     {
       for (int k = 0; k < (int) MPP_problem->v_opt_dishes[j].size(); k++)
       {
	     Neighbor n;
	     n.variable = j;
	     n.newValue = k;
	     neighbors.push_back(n);
       }
     }
     double best_feasibility = 0.0;
     feasibility_day(best_solution, best_feasibility);
//     cout << "entra day " << best_feasibility<<endl;
     double current_feasibility = best_feasibility;
     vector<int> current_solution = best_solution;
     auto start = high_resolution_clock::now();
     for (int i = 0; i < 2000; i++)
     {
	//first improvement hc
 	First_Improvement_Hill_Climbing_Day(neighbors, current_solution, current_feasibility);
	///save the best....
 	if( current_feasibility < best_feasibility)
	{
	//   cout << current_feasibility<<endl;
	   best_feasibility = current_feasibility;
	   best_solution = current_solution;
	}
	else
	{
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
//     cout << duration.count()/1.0e6 << endl;
//     cout << "sale--- " << best_feasibility<<endl;
     //exit(0);
    for(int i = 0; i  < N_OPT_DAY; i++) solution[day*N_OPT_DAY + i] = best_solution[i];
}
void  MPP::First_Improvement_Hill_Climbing_Day(vector<Neighbor> &neighbors, vector<int> &best_solution, double &best_feasibility)
{
   //incremental evaluation values...
   double current_feasibility = best_feasibility;
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   int &max_description_id = MPP_problem->max_description_id;
   vector<bool> uniq_per_day(max_description_id, false);

   int num_nutr = (int)MPP_problem->v_constraints.size();
   vector<vector<double>> nut_info;
   for(int i = 0 ; i  < N_OPT_DAY; i++) uniq_per_day[v_opt_dishes[i][best_solution[i]].description] = true;
   init_inc_eval_day(best_solution, nut_info);
   bool improved = true;
   while(improved)
   {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for (int i = 0; i < neighbors.size(); i++)
     {
	int day = neighbors[i].variable/N_OPT_DAY, opt =  neighbors[i].variable;
        struct infoDishes &dish_in = MPP_problem->v_opt_dishes[opt][neighbors[i].newValue];
        struct infoDishes &dish_out = MPP_problem->v_opt_dishes[opt][best_solution[opt]];
	if(uniq_per_day[dish_in.description]) continue;
        //incremental evaluation...
	double new_feasibility = inc_eval_day(neighbors[i], nut_info , best_solution, best_feasibility);
	if( new_feasibility < best_feasibility)
	{
	   improved = true;
	   update_inc_day(nut_info, neighbors[i], best_solution);
	   best_feasibility = new_feasibility;
	   best_feasibility = max(0.0, best_feasibility);
	  // feasibility_day(best_solution, new_feasibility);
	  // cout << fabs(new_feasibility-best_feasibility)<<endl;
	   uniq_per_day[dish_in.description] = true;
	   uniq_per_day[dish_out.description] = false;
	}
      }
    }
}
double MPP::localSearch_testing_time(double finalTime) {

     vector<Neighbor> neighbors;
     for (int i = 0; i < nDias; i++)
     {
	for (int j = 0; j <N_OPT_DAY; j++)
	{
	   for (int k = 0; k < (int) MPP_problem->v_opt_dishes[j].size(); k++)
	   {
	     Neighbor n;
	     n.variable = i * N_OPT_DAY + j;
	     n.newValue = k;
	     neighbors.push_back(n);
	   }
        }
     }

     vector<Neighbor_swap> neighbors_swap;
     for(int i = 0; i < nDias; i++)
     {
      for(int j = i+1; j < nDias; j++)
      {
         neighbors_swap.push_back({i, j});
      }
     }

     calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
     calculateVariability(x_var, obj_values);

     vector<int> bestIndividual = x_var;
     vector<double> best_objs = obj_values;
     //load incremental evaluation values...
 //    cout << "entra--- " << obj_values<<endl;
     struct timeval currentTime;
     gettimeofday(&currentTime, NULL);
     double initialTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
     double elapsedTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6;
     elapsedTime -= initialTime;
     while(elapsedTime < finalTime )
     {
        calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
        calculateVariability(x_var, obj_values);
        First_Improvement_Hill_Climbing(neighbors, x_var, obj_values);
        First_Improvement_Hill_Climbing_swap(neighbors_swap, x_var, obj_values);
        if(comp_objs(obj_values, best_objs))
        {
	  best_objs = obj_values;
          bestIndividual = x_var;
//	  cout << obj_values<< " "<< " -- " << badDaysFeas.size() << endl;
        }
        else x_var = bestIndividual;
        if(best_objs[0] == 0.0 && best_objs[1] == 0.0)
	break;
        calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
        calculateVariability(x_var, obj_values);

	if(badDaysFeas.empty())
        {
	   int selectedDay = -1;
	   if (heaviestNut != -1)
           {
		vector< pair<double, int> > infoNut;
		for (int i = 0; i < nDias; i++)
		{
		   double total = 0.0;
		   for(int k = 0; k < N_OPT_DAY; k++) total +=  MPP_problem->v_opt_dishes[k][x_var[i*N_OPT_DAY + k]].v_nutrient_value[heaviestNut];
		   infoNut.push_back(make_pair(total, i));
		}
		sort(infoNut.begin(), infoNut.end());
		if (heaviestType == 1) reverse(infoNut.begin(), infoNut.end());
		int nbestday = rand()%min(nDias ,  6);
		selectedDay = infoNut[nbestday].second;
	   }
	   else
	   {
		if( rand()%2 && !badDaysVar.empty())
		{
	 	   vector<int> tmp;
	           for(auto it = badDaysVar.begin(); it != badDaysVar.end(); it++) tmp.push_back(*it);
		   random_shuffle(tmp.begin(), tmp.end());
		   selectedDay = tmp[0];
		}
		else
		 selectedDay = rand() % nDias;
	   }
	   perturb_day(x_var, selectedDay);
	   oneDaylocalSearch(x_var, selectedDay);
	 }
	 else
	 {
	   for (auto it = badDaysFeas.begin(); it != badDaysFeas.end(); it++)
	   {
	      perturb_day(x_var, *it);
	      oneDaylocalSearch(x_var, *it);
	   }
 	 }
	gettimeofday(&currentTime, NULL);
	elapsedTime = ((double) (currentTime.tv_sec) + (double) (currentTime.tv_usec)/1.0e6)-initialTime;
     }
     x_var = bestIndividual;
     evaluate();
//     cout << "sale--- " << obj_values <<" fitness.."<< fitness <<endl;
    return elapsedTime;
}
