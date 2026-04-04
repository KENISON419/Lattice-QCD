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
extern void init_genrand(unsigned long s);
extern double genrand_real2(void);

//--------------------------------------------------------------------------

//--------------------------------------------------------------------------

void mc_init(vector<int>& spin, int flg_init)
{
  init_genrand((unsigned long)time(NULL));

  for (int i = 0; i < V; i++) {
    if (flg_init == 0) {
      spin[i] = 1;
    }
    else {
      spin[i] = (genrand_real2() < 0.5) ? -1 : 1;
    }
  }
}

//--------------------------------------------------------------------------

void mc_update_1d(vector<int>& spin)
{
  for (int n = 0; n < V; n++) {
    const int x = (int)(genrand_real2() * L);
    const int xm = (x - 1 + L) % L;
    const int xp = (x + 1) % L;

    const int s = spin[x];
    const int nn_sum = spin[xm] + spin[xp];
    const double dE = 2.0 * s * (J * nn_sum + mu0 * H);

    if (dE <= 0.0 || genrand_real2() < exp(-beta_inv * dE)) {
      spin[x] = -s;
    }
  }
}

//--------------------------------------------------------------------------

void mc_update_2d(vector<int>& spin)
{
  for (int n = 0; n < V; n++) {
    const int x = (int)(genrand_real2() * L);
    const int y = (int)(genrand_real2() * L);

    const int xm = (x - 1 + L) % L;
    const int xp = (x + 1) % L;
    const int ym = (y - 1 + L) % L;
    const int yp = (y + 1) % L;

    const int i = x + L * y;
    const int ixm = xm + L * y;
    const int ixp = xp + L * y;
    const int iym = x + L * ym;
    const int iyp = x + L * yp;

    const int s = spin[i];
    const int nn_sum = spin[ixm] + spin[ixp] + spin[iym] + spin[iyp];
    const double dE = 2.0 * s * (J * nn_sum + mu0 * H);

    if (dE <= 0.0 || genrand_real2() < exp(-beta_inv * dE)) {
      spin[i] = -s;
    }
  }
}

//--------------------------------------------------------------------------

void mc_update(vector<int>& spin)
{
  if (n_dim == 1) {
    mc_update_1d(spin);
  }
  else if (n_dim == 2) {
    mc_update_2d(spin);
  }
  else {
    fprintf(stderr, "n_dim = %d is not supported.\n", n_dim);
    exit(1);
  }
}

//--------------------------------------------------------------------------
