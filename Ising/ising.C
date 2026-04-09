#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<string>
#include<vector>
#include "ising.h"

using namespace std;

//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
extern void mc_init  (vector<int>&, int);
extern void mc_update(vector<int>&);

extern double calc_mag_mc(vector<int>&);
extern double calc_E_mc  (vector<int>&);

extern void calc_mag    (vector<double>&, vector<double>&);
extern void calc_abs_mag(vector<double>&, vector<double>&);
extern void calc_mag_chi(vector<double>&, vector<double>&);
extern void calc_E      (vector<double>&, vector<double>&);
extern void calc_C_heat (vector<double>&, vector<double>&);

//extern void stat_estimate(vector<double>&, vector<double>&);
//--------------------------------------------------------------------------

void ising(vector<double>& mag, vector<double>& abs_mag, vector<double>& mag_chi,
           vector<double>& E,   vector<double>& C_heat,
           vector<double>& therm_mag,
           int flg_init, int n_therm, int n_sweep, int n_mc)
{
  vector<int> spin(V);
  vector<double> mag_mc(n_mc);
  vector<double> E_mc  (n_mc);
  therm_mag.clear();
  therm_mag.reserve(n_therm);

  mc_init(spin, flg_init);

  for(int i=0; i<n_therm; i++) {
    mc_update(spin);
    therm_mag.push_back(calc_mag_mc(spin));
  }

  for(int i_mc=0; i_mc<n_mc; i_mc++) {
    for(int i=0; i<n_sweep; i++) {
      mc_update(spin);
    }
    mag_mc[i_mc] = calc_mag_mc(spin);
    E_mc  [i_mc] = calc_E_mc  (spin);
  }

  calc_mag    (mag,     mag_mc);
  calc_abs_mag(abs_mag, mag_mc);
  calc_mag_chi(mag_chi, mag_mc);
  calc_E      (E,       E_mc);
  calc_C_heat (C_heat,  E_mc);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
