#ifndef __MPP_H__
#define __MPP_H__
#include "MPP_Problem.h"

using namespace std;
struct Food {
   vector<int> p;
   multiset<int> sP;
   int n;
   Food(){}
   Food(int _n){
	this->n=_n;
	p.assign(n,0);
   }
   bool operator<(const Food &i2) const {
	 return sP < i2.sP;
   //   for(int i = 0; i < this->n; i++){
   //     if(p[i] < i2.p[i]) return true;
   //     else if(p[i] > i2.p[i]) return false;
   //   }
   //   return false; //dishes are equal.
   }
};


struct Neighbor{
   int variable;
   int newValue;
};
struct Neighbor_swap{
   int day1, day2, day3, day4;
};
struct Solution_LS{
   set<int> badDays;
   int heaviestNut, heaviestType;
   vector<double> objFeasibility;
   vector<pair<double, double> > objVariability;
   vector<int> x_var;
   vector< vector<double> > globalPlan;
   vector< vector< vector<double> > > nutriment_per_day;
   vector< vector< vector < int > > > time_id_day_table;
   vector< vector< int > > time_diff, time_diff_cat;
   vector<vector<int>> uniq_per_day;
};
class MPP:public MPP_Problem{
   public:
    MPP(){
    }
    ~MPP(){}
     //These functions are defined in "Model.cpp"
    void evaluate();
    void evaluate(vector<int> &sol, vector<double> &objFeas, vector<pair<double, double> > &objVar);
    double modelation(vector<double> &objFeas, vector<pair<double, double> > &objVar);

    //These functions are dfined in "recomb.cpp"
    void dependentMutation(double pm);
    void dependentCrossover(MPP &i2);
    void uniformCrossover(MPP &i2);
    void uniform2Crossover(MPP &i2);
    void pairBasedCrossover(MPP &i2);

    //These functions are defined in "localsearches.cpp"
    double localSearch(double finalTime, bool isMeasuringTime=false);
    void  First_Improvement_Hill_Climbing(vector<Neighbor> &neighbors, vector<int> &best_sol, vector<double> &bestFeas, vector<pair<double, double> > &bestVar);

    void First_Improvement_Hill_Climbing_swap(vector<Neighbor_swap> &neighbors, vector<int> &best_sol, vector<double> &bestFeas, vector<pair<double, double> > &bestVar);
    void First_Improvement_Hill_Climbing_swap2(vector<Neighbor_swap> &neighbors, vector<int> &best_sol, vector<double> &bestFeas, vector<pair<double, double> > &bestVar);


    void  First_Improvement_Hill_Climbing_Day(vector<Neighbor> &oneDayneighbors, vector<int> &best_solution, double &best_feasibility);
    void oneDaylocalSearch(vector<int> &solution, int day);

    //These functinos are defined in "MPP.cpp"
    void restart();
    void init (); 
    void feasibility_day(vector<int> &best_solution, double &best_feasibility);
		
    //functions that initialize the incremental evaluation of the local search in each day.
    void init_inc_eval_day(vector<int> &current_solution, vector<vector<double>> &nut_info);
    double inc_eval_day(Neighbor &neighbor, vector<vector<double>> &nut_info, vector<int> &best_solution, double best_feasibility);
    void update_inc_day(vector<vector<double>> &nut_info, Neighbor &new_neighbor, vector<int> &best_solution);
    int getDistance(MPP &ind2); 
    virtual void print(ostream &os) const;
    void exportcsv(){
       MPP_problem->exportcsv(x_var);
    }
    vector<int> x_var;
    double fitness;
    static MPP_Problem *MPP_problem;

  private:
    void calculateFeasibilityDegree(vector<int> &sol, double &feas_day, double &feas_global);
    void init_incremental_evaluation(struct Solution_LS &current);
    void inc_eval(struct Solution_LS &current, Neighbor &new_neighbor, vector<double> &new_objs, vector<pair<double, double> > &newVariability);
    void update_inc(struct Solution_LS &current, Neighbor &neighbor, vector<double> &newFeas, vector<pair<double, double> > &newVar);
    void swap_days(vector<int> &data, int day1, int day2);
    void perturb_opt(vector<int> &data, int day, int which);
    void perturb_day(vector<int> &data, int day);
    void update_dcn_pair(int diff, pair<int, int> &p_dcn);
    void calculateVariability(vector<int> &sol, vector<pair<double, double> > &objVar);
    bool comp_objs(vector<double> &objFeas1, vector<pair<double, double> > &objVar1, vector<double> &objFeas2, vector<pair<double, double> > &objVar2);
    //calculate neighbourhoods that are considered in each local search..
    void init_neighbours();

    int heaviestNut, heaviestType;
    double valorFac, variabilidadObj;//factibility and variability of the current solution..
    vector<double> objFeasibility;//{feasiblity, varibility x times}
    vector<pair<double, double> > objVariability;
    set<int> badDaysFeas, badDaysVar;
    //These vectors save neighbourhood which are considered into the local searches
    vector<Neighbor> neighbors, oneDayneighbors;
    vector<Neighbor_swap> neighbors_swap, neighbors_swap2;

   void Greedy_Improvement(vector<int> &x_var, vector<double> &objFeasibility, vector<pair<double, double> > &objVariability);

   void Improvement(vector<int> &best_sol, vector<double> &bestFeas, vector<pair<double, double> > &bestVar);
};
#endif
