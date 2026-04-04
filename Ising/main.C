#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<string>
#include<vector>


using namespace std;

//--------------------------------------------------------------------------
int n_dim = 0;
int L = 0;
int V = 0;
double J = 0.0;
double mu0 = 0.0;
double T = 0.0;
double H = 0.0;
double beta = 0.0; // 1/T
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
static void set_args(int argc, char* argv[],
                     int& n_dim, int& L, double& J, double& mu0,
                     int& flg_init, int& n_therm, int& n_sweep, int& n_mc,
                     double& T_i, double& T_f, int& n_T,
                     double& H_i, double& H_f, int& n_H);
static void usage(void);
//--------------------------------------------------------------------------
void ising(vector<double>& mag, vector<double>& mag_chi,
           vector<double>& E,   vector<double>& C_heat,
           int flg_init, int n_therm, int n_seep, int n_mc);
//--------------------------------------------------------------------------

int
main(int argc, char* argv[])
{
  int flg_init;
  int n_therm, n_sweep, n_mc;
  double T_i, T_f;
  double H_i, H_f;
  int n_T, n_H;

  set_args(argc, argv, n_dim, L, J, mu0,
           flg_init, n_therm, n_sweep, n_mc,
           T_i, T_f, n_T, H_i, H_f, n_H);

  V = 1; // volume
  for(int i_dim = 0; i_dim < n_dim; i_dim++) {
    V *= L;
  }

  printf("# ndim= %d L= %d V= %d J= %f mu0= %f\n", n_dim, L, V, J, mu0);
  printf("# flg_init=%d n_therm=%d n_sweep=%d n_mc= %d\n", flg_init, n_therm, n_sweep, n_mc);
  printf("# T= %f %f %d H= %f %f %d\n", T_i, T_f, n_T, H_i, H_f, n_H);

  double T_sep = 0;
  if ( n_T != 0 ) {
    T_sep = (T_f - T_i) / n_T;
  }

  double H_sep = 0;
  if ( n_H != 0 ) {
    H_sep = (H_f - H_i) / n_H;
  }


  for(  int iT = 0; iT <= n_T; iT++) {
    for(int iH = 0; iH <= n_H; iH++) {
      T = T_i + iT * T_sep;
      H = H_i + iH * H_sep;
      beta = 1.0/T;

      vector<double> mag(2); // cen, err
      vector<double> mag_chi(2);
      vector<double> E(2);
      vector<double> C_heat(2);
      ising(mag, mag_chi, E, C_heat, flg_init, n_therm, n_sweep, n_mc);

      printf("T= %f H= %f "
             "mag= %f mag_err= %f "
             "mag_chi= %f mag_chi_err= %f "
             "E= %f E_err= %f "
             "C_heat= %f C_heat_err= %f\n",
             T, H,
             mag[0], mag[1],
             mag_chi[0], mag_chi[1],
             E[0], E[1],
             C_heat[0], C_heat[1]);
    }
  }

  return 0;
}


//--------------------------------------------------------------------------

static void set_args(int argc, char* argv[],
                     int& n_dim, int& L, double& J, double& mu0,
                     int& flg_init, int& n_therm, int& n_sweep, int& n_mc,
                     double& T_i, double& T_f, int& n_T,
                     double& H_i, double& H_f, int& n_H)
{
  if (argc == 1) usage();
  argc--; argv++;

  //--------------------------------------------------------------------------
  // analysing args
  //--------------------------------------------------------------------------

  while(argc > 0){
    if(strcmp(argv[0],"-n_dim")==0){
      if( argc < 2 ) { usage(); }
      n_dim = atoi(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-L")==0){
      if( argc < 2 ) { usage(); }
      L = atoi(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-flg_init")==0){
      if( argc < 2 ) { usage(); }
      flg_init = atoi(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-n_therm")==0){
      if( argc < 2 ) { usage(); }
      n_therm = atoi(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-n_sweep")==0){
      if( argc < 2 ) { usage(); }
      n_sweep = atoi(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-n_mc")==0){
      if( argc < 2 ) { usage(); }
      n_mc = atoi(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-J")==0){
      if( argc < 2 ) { usage(); }
      J = atof(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-mu0")==0){
      if( argc < 2 ) { usage(); }
      mu0 = atof(argv[1]);
      argc-=2; argv+=2;
      continue;
    }
    if(strcmp(argv[0],"-T")==0){
      if( argc < 4 ) { usage(); }
      T_i = atof(argv[1]);
      T_f = atof(argv[2]);
      n_T = atoi(argv[3]);
      argc-=4; argv+=4;
      continue;
    }
    if(strcmp(argv[0],"-H")==0){
      if( argc < 4 ) { usage(); }
      H_i = atof(argv[1]);
      H_f = atof(argv[2]);
      n_H = atoi(argv[3]);
      argc-=4; argv+=4;
      continue;
    }
    // unknown option
    argc-=1; argv+=1;
  }

}

//--------------------------------------------------------------------------
static void usage(void)
{
  fprintf(stderr,
          "usage: [OPTIONS]\n"
          "\t -n_dim <n_dim>\n"
          "\t -L <L>\n"
          "\t -J <J>\n"
          "\t -mu0 <mu0>\n"
          "\t -flg_init <0: cold, 1: hot>\n"
          "\t -n_therm <n_therm>\n"
          "\t -n_sweep <n_sweep>\n"
          "\t -n_mc <n_mc>\n"
          "\t -T <T_i T_f n_T>\n"
          "\t -H <H_i H_f n_H>\n"
          );

  exit(1);
}
//--------------------------------------------------------------------------

