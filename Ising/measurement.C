#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<string>
#include<vector>
#include<time.h>
#include "ising.h"

using namespace std;

//--------------------------------------------------------------------------
//extern void stat_estimate(vector<double>&, vector<double>&);
extern void calc_jackknife_nosuberr(vector<double>& dat_jk, vector<double>& dat);
extern void calc_jackknife_final(vector<double>& res, vector<double>& dat_jk);

//--------------------------------------------------------------------------

double calc_mag_mc(vector<int>& spin)
{

  // edit

  return 0.0; // dummy
}

//--------------------------------------------------------------------------

void calc_mag(vector<double>& mag, vector<double>& mag_mc)
{

  // edit
  
}

//--------------------------------------------------------------------------

void calc_mag_chi(vector<double>& mag_chi, vector<double>& mag_mc)
{

  // edit
  
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

double calc_E_mc_1d(vector<int>& spin)
{
  // edit

  return 0.0; // dummy
}

double calc_E_mc_2d(vector<int>& spin)
{

  // edit

  return 0.0; // dummy
}

double calc_E_mc(vector<int>& spin)
{
  double E;
  
  if (n_dim == 1) {
    E = calc_E_mc_1d(spin);
  }
  else if (n_dim == 2) {
    E = calc_E_mc_2d(spin);
  }
  else {
    fprintf(stderr, "n_dim = %d is not supported.\n", n_dim);
    exit(1);
  }

  return E;
}

//--------------------------------------------------------------------------

void calc_E(vector<double>& E, vector<double>& E_mc)
{

  // edit

}

//--------------------------------------------------------------------------

void calc_C_heat(vector<double>& C_heat, vector<double>& E_mc)
{

  // edit
  
}

//--------------------------------------------------------------------------
