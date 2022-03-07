#include <bits/stdc++.h>
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

extern double px;
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
		 int max_description_id;
		 int nDias;
};

