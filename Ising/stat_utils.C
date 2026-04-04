#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<string>
#include<vector>

using namespace std;

//--------------------------------------------------------------------------

void stat_estimate(vector<double>& res, vector<double>& dat)
{
  const int n = dat.size();
  if (n <= 0) {
    res[0] = 0.0;
    res[1] = 0.0;
    return;
  }

  double ave = 0.0;
  for (int i = 0; i < n; i++) {
    ave += dat[i];
  }
  ave /= n;

  double var = 0.0;
  for (int i = 0; i < n; i++) {
    const double d = dat[i] - ave;
    var += d * d;
  }
  var /= n;

  res[0] = ave;
  res[1] = (n > 1) ? sqrt(var / (n - 1.0)) : 0.0;
}

//--------------------------------------------------------------------------

void calc_jackknife_nosuberr(vector<double>& dat_jk, vector<double>& dat)
{
  const int n = dat.size();
  if (n <= 1) {
    dat_jk[0] = (n == 1) ? dat[0] : 0.0;
    return;
  }

  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += dat[i];
  }

  for (int i = 0; i < n; i++) {
    dat_jk[i] = (sum - dat[i]) / (n - 1.0);
  }
}

void calc_jackknife_final(vector<double>& res, vector<double>& dat_jk)
{
  const int n = dat_jk.size();
  if (n <= 0) {
    res[0] = 0.0;
    res[1] = 0.0;
    return;
  }

  double ave = 0.0;
  for (int i = 0; i < n; i++) {
    ave += dat_jk[i];
  }
  ave /= n;

  double var = 0.0;
  for (int i = 0; i < n; i++) {
    const double d = dat_jk[i] - ave;
    var += d * d;
  }
  var *= (n - 1.0) / n;

  res[0] = ave;
  res[1] = sqrt(var);
}
//--------------------------------------------------------------------------
