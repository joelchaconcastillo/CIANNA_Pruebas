#include <bits/stdc++.h>
//v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner, v_both_snack;
#define EPSILON 1e-10

#define N_CATEGORIES 3
#define CATEGORY_1 1
#define CATEGORY_2 2
#define CATEGORY_BOTH 0
#define GLOBAL 1
#define DIARIA 2

//encoded times by each day of the individual
#define N_TIMES 6
#define N_OPT_DAY 8
#define BREAKFAST 0
#define MORNING_SNACK 1
#define STARTER_1 2
#define STARTER_2 3
#define MAIN_COURSE_1 4
#define MAIN_COURSE_2 5
#define EVENING_SNACK 6
#define DINNER 7



//#define BOTH_SNACK 8
////////crossover type......
#define PAIR_BASED_CROSSOVER 1
#define UNIFORM_CROSSOVER 2
#define UNIFORM2_CROSSOVER 3
#define WEIGHT_DAY 1.0e6
#define DAYS_FAVORITE 7*3
#define DAYS_NO_FAVORITE 7*4
#define ITERATIONS_LS 100

#define W_VAR_GLOBAL 100
#define W_VAR_GLOBAL_CAT 1
extern int crossoverType;
extern int nDias;

using namespace std;
struct infoDishes {
        int description;	
	string time_day; //time related with the dish
	vector<double> v_nutrient_value;    //nutriments meta-data...those indexes should be in the same order than the v_contraints vector...
	int category; //category, at this point is 1 or 2
	bool favorite; //true if this is a favorite dish...
};
struct constraint_nutrient
{
   double min, max;
   string name;
   int type;
};
class MPP_Problem{
	public:
		MPP_Problem();
		~MPP_Problem(){
		}
 		void load_data(int argc, char **argv);
		void load_constraints(char *Plates_file);
		void load_dishes(char *Constraints_file);
		void exportcsv(vector<int> &x_var);

		int random_dish(int time_dish);

		vector<vector<infoDishes> > v_opt_dishes;  // the same as ->  vector<int> v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner;
		vector<constraint_nutrient> v_constraints;
		unordered_map<string, int> dic_nut_id;
		vector<int> v_constraint_global, v_constraint_day;
		string out_filename;
		vector<vector<int>> conf_day, opt_conf, unique_opt_time;
		vector<int> inv_unique_opt_time;
		vector<double> weights;
		int max_description_id;
};
