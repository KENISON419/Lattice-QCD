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

  // edit here

}

//--------------------------------------------------------------------------

void mc_update_1d(vector<int>& spin)
{

  // edit
  
}

//--------------------------------------------------------------------------

void mc_update_2d(vector<int>& spin)
{

  // edit

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
