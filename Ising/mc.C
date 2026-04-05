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
  static bool is_seeded = false;
  if (!is_seeded) {
    init_genrand((unsigned long)time(NULL));
    is_seeded = true;
  }

  for (int i = 0; i < V; i++) {
    if (flg_init == 0) {
      spin[i] = 1; // cold start
    }
    else {
      spin[i] = (genrand_real2() < 0.5) ? -1 : 1; // hot start
    }
  }

}

//--------------------------------------------------------------------------

void mc_update_1d(vector<int>& spin)
{
  for (int it = 0; it < V; it++) {
    int i = (int)(genrand_real2() * V);
    if (i >= V) i = V - 1;

    int ip = (i + 1) % L;
    int im = (i - 1 + L) % L;

    int s = spin[i];
    int neigh = spin[ip] + spin[im];
    double dE = 2.0 * s * (J * neigh + mu0 * H);

    if (dE <= 0.0 || genrand_real2() < exp(-beta_inv * dE)) {
      spin[i] = -s;
    }
  }
  
}

//--------------------------------------------------------------------------

void mc_update_2d(vector<int>& spin)
{
  for (int it = 0; it < V; it++) {
    int idx = (int)(genrand_real2() * V);
    if (idx >= V) idx = V - 1;

    int x = idx % L;
    int y = idx / L;

    int xp = (x + 1) % L;
    int xm = (x - 1 + L) % L;
    int yp = (y + 1) % L;
    int ym = (y - 1 + L) % L;

    int s = spin[idx];
    int neigh =
      spin[y * L + xp] + spin[y * L + xm] +
      spin[yp * L + x] + spin[ym * L + x];

    double dE = 2.0 * s * (J * neigh + mu0 * H);

    if (dE <= 0.0 || genrand_real2() < exp(-beta_inv * dE)) {
      spin[idx] = -s;
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
