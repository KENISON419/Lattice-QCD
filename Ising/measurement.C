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
  double sum_s = 0.0;
  for (int i = 0; i < V; i++) {
    sum_s += spin[i];
  }

  return mu0 * sum_s / V;
}

//--------------------------------------------------------------------------

void calc_mag(vector<double>& mag, vector<double>& mag_mc)
{
  vector<double> mag_jk;
  calc_jackknife_nosuberr(mag_jk, mag_mc);
  calc_jackknife_final(mag, mag_jk);
  
}

void calc_abs_mag(vector<double>& abs_mag, vector<double>& mag_mc)
{
  const int n = (int)mag_mc.size();
  vector<double> abs_mag_mc(n, 0.0);
  for (int i = 0; i < n; i++) {
    abs_mag_mc[i] = fabs(mag_mc[i]);
  }

  vector<double> abs_mag_jk;
  calc_jackknife_nosuberr(abs_mag_jk, abs_mag_mc);
  calc_jackknife_final(abs_mag, abs_mag_jk);
}

//--------------------------------------------------------------------------

void calc_mag_chi(vector<double>& mag_chi, vector<double>& mag_mc)
{
  const int n = (int)mag_mc.size();
  vector<double> chi_jk(n, 0.0);

  for (int i = 0; i < n; i++) {
    double sum_m = 0.0;
    double sum_m2 = 0.0;
    for (int j = 0; j < n; j++) {
      if (j == i) continue;
      sum_m += mag_mc[j];
      sum_m2 += mag_mc[j] * mag_mc[j];
    }
    const double den = (double)(n - 1);
    const double ave_m = sum_m / den;
    const double ave_m2 = sum_m2 / den;
    chi_jk[i] = beta_inv * V * (ave_m2 - ave_m * ave_m);
  }

  calc_jackknife_final(mag_chi, chi_jk);
  
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

double calc_E_mc_1d(vector<int>& spin)
{
  double E_bond = 0.0;
  double sum_s = 0.0;

  for (int i = 0; i < L; i++) {
    int ip = (i + 1) % L;
    E_bond += -J * spin[i] * spin[ip];
    sum_s += spin[i];
  }

  double E_field = -mu0 * H * sum_s;
  double E_tot = E_bond + E_field;

  return E_tot / V;
}

double calc_E_mc_2d(vector<int>& spin)
{
  double E_bond = 0.0;
  double sum_s = 0.0;

  for (int y = 0; y < L; y++) {
    for (int x = 0; x < L; x++) {
      int idx = y * L + x;
      int xp = (x + 1) % L;
      int yp = (y + 1) % L;

      E_bond += -J * spin[idx] * spin[y * L + xp];
      E_bond += -J * spin[idx] * spin[yp * L + x];
      sum_s += spin[idx];
    }
  }

  double E_field = -mu0 * H * sum_s;
  double E_tot = E_bond + E_field;

  return E_tot / V;
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
  vector<double> E_jk;
  calc_jackknife_nosuberr(E_jk, E_mc);
  calc_jackknife_final(E, E_jk);

}

//--------------------------------------------------------------------------

void calc_C_heat(vector<double>& C_heat, vector<double>& E_mc)
{
  const int n = (int)E_mc.size();
  vector<double> C_jk(n, 0.0);

  for (int i = 0; i < n; i++) {
    double sum_E = 0.0;
    double sum_E2 = 0.0;
    for (int j = 0; j < n; j++) {
      if (j == i) continue;
      sum_E += E_mc[j];
      sum_E2 += E_mc[j] * E_mc[j];
    }
    const double den = (double)(n - 1);
    const double ave_E = sum_E / den;
    const double ave_E2 = sum_E2 / den;
    C_jk[i] = beta_inv * beta_inv * V * (ave_E2 - ave_E * ave_E);
  }

  calc_jackknife_final(C_heat, C_jk);
  
}

//--------------------------------------------------------------------------
