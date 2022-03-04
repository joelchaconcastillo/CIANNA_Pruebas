/*
 * Developers: Joel Chacón Castillo, Carlos Segura
 * Last modification: 03/05/22
 * This file defines the model that is taken into account as the fitness of each individuals given a combination of dishes
 * */

#include "MPP.h"
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
	max1 = max(max1,fabs(v1-ref)/MPP_problem->weights[a-2]);
  }
  fitness = 1e4*(obj_values[0]+obj_values[1]) + max1;
}
