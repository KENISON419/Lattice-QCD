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
  double m = 0.0;
  for (int i = 0; i < V; i++) {
    m += spin[i];
  }
  return mu0 * m / V;
}

//--------------------------------------------------------------------------

void calc_mag(vector<double>& mag, vector<double>& mag_mc)
{
  const int n_mc = mag_mc.size();
  vector<double> abs_m(n_mc);
  for (int i = 0; i < n_mc; i++) {
    abs_m[i] = fabs(mag_mc[i]);
  }

  vector<double> abs_m_jk(n_mc);
  calc_jackknife_nosuberr(abs_m_jk, abs_m);
  calc_jackknife_final(mag, abs_m_jk);
}

//--------------------------------------------------------------------------

void calc_mag_chi(vector<double>& mag_chi, vector<double>& mag_mc)
{
  const int n_mc = mag_mc.size();
  vector<double> m_jk(n_mc);
  vector<double> m2(n_mc), m2_jk(n_mc);
  vector<double> chi_jk(n_mc);

  for (int i = 0; i < n_mc; i++) {
    m2[i] = mag_mc[i] * mag_mc[i];
  }

  calc_jackknife_nosuberr(m_jk, mag_mc);
  calc_jackknife_nosuberr(m2_jk, m2);

  for (int i = 0; i < n_mc; i++) {
    chi_jk[i] = beta_inv * V * (m2_jk[i] - m_jk[i] * m_jk[i]);
  }

  calc_jackknife_final(mag_chi, chi_jk);
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

double calc_E_mc_1d(vector<int>& spin)
{
  double E = 0.0;
  double M = 0.0;

  for (int x = 0; x < L; x++) {
    const int xp = (x + 1) % L;
    E += -J * spin[x] * spin[xp];
    M += spin[x];
  }
  E += -mu0 * H * M;

  return E / V;
}

double calc_E_mc_2d(vector<int>& spin)
{
  double E = 0.0;
  double M = 0.0;

  for (int y = 0; y < L; y++) {
    for (int x = 0; x < L; x++) {
      const int i = x + L * y;
      const int xp = (x + 1) % L;
      const int yp = (y + 1) % L;
      const int ixp = xp + L * y;
      const int iyp = x + L * yp;

      E += -J * spin[i] * spin[ixp];
      E += -J * spin[i] * spin[iyp];
      M += spin[i];
    }
  }
  E += -mu0 * H * M;

  return E / V;
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
  const int n_mc = E_mc.size();
  vector<double> E_jk(n_mc);
  calc_jackknife_nosuberr(E_jk, E_mc);
  calc_jackknife_final(E, E_jk);
}

//--------------------------------------------------------------------------

void calc_C_heat(vector<double>& C_heat, vector<double>& E_mc)
{
  const int n_mc = E_mc.size();
  vector<double> e_jk(n_mc);
  vector<double> e2(n_mc), e2_jk(n_mc);
  vector<double> c_jk(n_mc);

  for (int i = 0; i < n_mc; i++) {
    e2[i] = E_mc[i] * E_mc[i];
  }

  calc_jackknife_nosuberr(e_jk, E_mc);
  calc_jackknife_nosuberr(e2_jk, e2);

  for (int i = 0; i < n_mc; i++) {
    c_jk[i] = beta_inv * beta_inv * V * (e2_jk[i] - e_jk[i] * e_jk[i]);
  }

  calc_jackknife_final(C_heat, c_jk);
}

//--------------------------------------------------------------------------
