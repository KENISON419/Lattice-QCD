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
  const int n = (int)dat.size();
  if (n <= 0) {
    res.assign(2, 0.0);
    return;
  }

  double ave = 0.0;
  for (int i = 0; i < n; i++) ave += dat[i];
  ave /= n;

  double var = 0.0;
  for (int i = 0; i < n; i++) {
    double d = dat[i] - ave;
    var += d * d;
  }
  var /= n;

  double err = sqrt(var / n);

  res.resize(2);
  res[0] = ave;
  res[1] = err;

}

//--------------------------------------------------------------------------

void calc_jackknife_nosuberr(vector<double>& dat_jk, vector<double>& dat)
{
  const int n = (int)dat.size();
  dat_jk.assign(n, 0.0);
  if (n <= 1) return;

  double sum = 0.0;
  for (int i = 0; i < n; i++) sum += dat[i];

  for (int i = 0; i < n; i++) {
    dat_jk[i] = (sum - dat[i]) / (n - 1);
  }

}

void calc_jackknife_final(vector<double>& res, vector<double>& dat_jk)
{
  const int n = (int)dat_jk.size();
  if (n <= 0) {
    res.assign(2, 0.0);
    return;
  }

  double ave = 0.0;
  for (int i = 0; i < n; i++) ave += dat_jk[i];
  ave /= n;

  double sq = 0.0;
  for (int i = 0; i < n; i++) {
    double d = dat_jk[i] - ave;
    sq += d * d;
  }
  double err = (n > 1) ? sqrt((n - 1.0) / n * sq) : 0.0;

  res.resize(2);
  res[0] = ave;
  res[1] = err;

}
//--------------------------------------------------------------------------
