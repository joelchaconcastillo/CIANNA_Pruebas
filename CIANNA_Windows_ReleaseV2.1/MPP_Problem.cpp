#include "MPP_Problem.h"

const string WHITESPACE = " \n\r\t\f\v";
string ltrim(const string& s)
{
	size_t start = s.find_first_not_of(WHITESPACE);
	return (start == std::string::npos) ? "" : s.substr(start);
}

string rtrim(const string& s)
{
	size_t end = s.find_last_not_of(WHITESPACE);
	return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

string trim(const string& s)
{
	return rtrim(ltrim(s));
}
MPP_Problem::MPP_Problem(){
}
int MPP_Problem::random_dish(int time_dish){
   return rand()%((int)v_opt_dishes[time_dish].size());
}
void MPP_Problem::load_data(int argc, char **argv)
{
    if(argc < 6)
    {
	 cout << "N\'umero de argumentos inv\'alidos" <<endl;
	 exit(EXIT_FAILURE);
    }
    out_filename = string(argv[5]);
    nDias = atoi(argv[3]);

    ////reading the information.......
    v_opt_dishes.assign(N_OPT_DAY, vector<infoDishes>()); //N_times_dishes...
    load_constraints(argv[2]);
    load_dishes(argv[1]);
}
void MPP_Problem::load_dishes(char *c_filename)
{
        max_description_id = 0;
   	ifstream ifs;
	struct infoDishes str_dish;
	ifs.open(c_filename, ifstream::in);
	if (ifs.is_open())
	{
	   string line, word;
	   getline(ifs, line);
	   stringstream line_commas(line);
	   vector<string> column_names;
	//Here is assumed that the first two columns are of the "DESCRIPCION" and "TIEMPO" respectively, thereafter are the meta-data and the last two columns are "CATEGORIA" and "FAVORITO" respectively.
	   while(getline(line_commas ,word, ',')) //first getting the tag-information of each column....
	   column_names.push_back(trim(word));
	   while (ifs.good())
	   {
	       string cell;

	       str_dish.v_nutrient_value.assign(v_constraints.size(), 0.0);
               for(int i = 0; i < column_names.size(); i++)
               {
		  if(i==column_names.size()-1)
	             getline(ifs, cell, '\n');
		  else
	             getline(ifs, cell, ',');
		  cell = trim(cell);
	       if(cell.empty()) break; //The file has an extra empty line
		  if( column_names[i] == "DESCRIPCION")
               		//str_dish.description = stoi(cell);
                    str_dish.description = atoi(cell.c_str());
		  else if( column_names[i] == "TIEMPO")
               		str_dish.time_day = trim(cell);
		  else if( column_names[i] == "CATEGORIA")
	       		str_dish.category = atoi(cell.c_str());
                //str_dish.category = stoi(cell);
		  else if( column_names[i] == "FAVORITO")
               		//str_dish.favorite = (bool)stoi(cell);
               		str_dish.favorite = (bool)atoi(cell.c_str());
		  else{
			if( dic_nut_id.find(column_names[i]) != dic_nut_id.end())
		 	{
			   //double value = stod(cell);
			   double value = atof(cell.c_str());
//			   value = floor(value * 10 + 0.5)/10;
		           str_dish.v_nutrient_value[dic_nut_id[column_names[i]]] = value; //the nutrient values need to be stored in the same order that the contraints..
			}
			//else { cout << column_names[i]<<endl;}// cout << "error interno"<<endl; exit(EXIT_FAILURE);}
			}
               }
	       max_description_id = max(max_description_id, str_dish.description);
	       if(cell.empty()) break; //The file has an extra empty line

	       if(str_dish.time_day == "DESAYUNO") v_opt_dishes[BREAKFAST].push_back(str_dish);
	       else if(str_dish.time_day == "COLACION_MATUTINA") v_opt_dishes[MORNING_SNACK].push_back(str_dish);
	   //    else if(str_dish.time_day == "COMIDA_ENTRADA") v_opt_dishes[STARTER].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_1) v_opt_dishes[STARTER_1].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_2) v_opt_dishes[STARTER_2].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_BOTH) //this dish blong to both categories...
		else if(str_dish.time_day == "COMIDA_ENTRADA" )
		{
		 v_opt_dishes[STARTER_1].push_back(str_dish);
		 v_opt_dishes[STARTER_2].push_back(str_dish);
		}
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL") v_opt_dishes[MAIN_COURSE].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_1) v_opt_dishes[MAIN_COURSE_1].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_2) v_opt_dishes[MAIN_COURSE_2].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_BOTH) //both categories..
		else if(str_dish.time_day == "COMIDA_PRINCIPAL")
		{
		  v_opt_dishes[MAIN_COURSE_1].push_back(str_dish);
		  v_opt_dishes[MAIN_COURSE_2].push_back(str_dish);
		}
	       else if(str_dish.time_day == "COLACION_VESPERTINA") v_opt_dishes[EVENING_SNACK].push_back(str_dish);
	       else if(str_dish.time_day == "CENA") v_opt_dishes[DINNER].push_back(str_dish);
	       else if(str_dish.time_day == "COLACION_AMBAS")
	       {
 		 v_opt_dishes[MORNING_SNACK].push_back(str_dish);
 		 v_opt_dishes[EVENING_SNACK].push_back(str_dish);
	       }
	       else
	       {
		   cout <<"Tiempo del dia desconocido"<<endl;
		   exit(EXIT_FAILURE);
	       }
//	      if( str_dish.category == CATEGORY_BOTH)
//	      {
//		str_dish.category = CATEGORY_1;
//	       v_dishes.push_back(str_dish);
//		str_dish.category = CATEGORY_2;
//	       v_dishes.push_back(str_dish);
//	      }
//	      else
//	       v_dishes.push_back(str_dish);
	    }
		ifs.close();
	} else {
		cout << "\n\nError. No se ha podido leer el archivo de platos."<<endl;
		exit(EXIT_FAILURE);
	}
//	cout << "platillos... " <<endl;
//       for(int i = 0; i < N_OPT_DAY; i++) cout << v_opt_dishes[i].size() << " ";
//       cout << endl;
}
void MPP_Problem::load_constraints(char *c_filename)
{
   	ifstream ifs;
	struct constraint_nutrient str_constraint_nutrient;
	ifs.open(c_filename, ifstream::in);
	if (ifs.is_open())
	{
	   string line, word;
	  //headers.....
	   getline(ifs, line);
	   stringstream line_commas(line);
	   vector<string> column_names;
	   while(getline(line_commas ,word, ',')) //first getting the tag-information of each column....
	   column_names.push_back(trim(word));
	   while (ifs.good())
	   {
	       string cell;
	       getline(ifs, cell, ',');
	       if(trim(cell).empty()) break; //The file has an extra empty line
		if( trim(cell) == "GLOBAL") str_constraint_nutrient.type=GLOBAL;
		else if( trim(cell) == "DIARIA") str_constraint_nutrient.type=DIARIA;
               //str_constraint_nutrient.type = trim(cell);

	       getline(ifs, cell, ',');
	       str_constraint_nutrient.name = trim(cell);
	       getline(ifs, cell, ',');
	       str_constraint_nutrient.min = atof(cell.c_str());
           //str_constraint_nutrient.min = stod(cell);
	       getline(ifs, cell, '\n'); //last word...
           //str_constraint_nutrient.max = stod(cell);
           str_constraint_nutrient.max = atof(cell.c_str());
	       dic_nut_id[str_constraint_nutrient.name] = (int) v_constraints.size();
	       if( str_constraint_nutrient.type == GLOBAL)
		{
		 str_constraint_nutrient.min *=nDias;
		 str_constraint_nutrient.max *=nDias;
		 v_constraint_global.push_back((int)v_constraints.size());
		}
	       else if( str_constraint_nutrient.type == DIARIA) v_constraint_day.push_back((int)v_constraints.size());
	       else{
		   cout << "Se desconoce un tipo de restricci\'on, \'unicamente puede ser por d\'ia o global"<<endl;
		   exit(EXIT_FAILURE);
		}
	       v_constraints.push_back(str_constraint_nutrient);
	    }
		ifs.close();
	} else {
		cout << "\n\nError. No se ha podido leer el archivo de restricciones."<<endl;
		exit(EXIT_FAILURE);
	}
}
void MPP_Problem::exportcsv(vector<int> &x_var)
{
   ofstream ofs;
   ofs.open(out_filename.c_str());
   ofs << "DIA ";
   set<int> times_selected_per_day;
   for(int i = 0; i < conf_day.size(); i++)
     for(auto t =conf_day[i].begin(); t != conf_day[i].end(); t++) times_selected_per_day.insert(*t);
   for(auto a = times_selected_per_day.begin(); a != times_selected_per_day.end(); a++)
   {
	   if(*a == BREAKFAST) ofs << " , DESAYUNO ";
	   if(*a == MORNING_SNACK) ofs << " , COLACION_MATUTINA ";
	   if(*a == STARTER_1) ofs << " , COMIDA_ENTRADA ";
	   if(*a == STARTER_2) ofs << " , COMIDA_ENTRADA ";
	   if(*a == MAIN_COURSE_1) ofs << " , COMIDA_PRINCIPAL ";
	   if(*a == MAIN_COURSE_2) ofs << " , COMIDA_PRINCIPAL ";
	   if(*a == EVENING_SNACK) ofs << " , COLACION_VESPERTINA ";
	   if(*a == DINNER) ofs << " , CENA ";
   }
//   ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
   //ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
    for(auto i = dic_nut_id.begin(); i !=  dic_nut_id.end(); i++) ofs <<" , "<<i->first ;
	ofs<< "\n";

   for(int i = 0; i < nDias; i++)
   {
	ofs << i+1;
   	for(auto a = times_selected_per_day.begin(); a != times_selected_per_day.end(); a++)
	  ofs<< " , "<<v_opt_dishes[(*a)][x_var[i*N_OPT_DAY + (*a)]].description;
    	   for(auto ij = dic_nut_id.begin(); ij !=  dic_nut_id.end(); ij++)
	   {
		ofs <<" , \" (";
		for(int c = 0; c < conf_day.size(); c++)
		{
		double sum_nut = 0.0;
     		  for(auto t =conf_day[c].begin(); t != conf_day[c].end(); t++)
	              sum_nut +=v_opt_dishes[*t][x_var[i*N_OPT_DAY + (*t)]].v_nutrient_value[ij->second];
			if( c>0) ofs<<",";
			ofs <<sum_nut ;
		}
 		 ofs<<") ["<<v_constraints[ij->second].min<<","<<v_constraints[ij->second].max<<"] \"" ;
	   }
	ofs<<"\n";
   }
   ofs.close();
}

