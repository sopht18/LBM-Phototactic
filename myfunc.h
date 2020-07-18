#include <iostream>
#include "parameters.h"
// // inline~ だと早くなるらしい(?)
inline void set_w(double w[9]) { // create w of D2Q9
  w[0] = 4.0/9.0;
  for (int i = 1; i <= 8; i++) {
    if (i <= 4) {
      w[i] = 1.0/9.0;
    }
    else{
      w[i] = 1.0/36.0;
    }
  }
}
void set_w(); // この宣言抜いても動いたが，ならこれは何故必要？  inline だとあってもなくても良いのか？

inline void set_e(double e[2][10]) { // // create vec{e}
  e[0][0] = 0.0;
  e[1][0] = 0.0;
  e[0][1] = (dx/dt);
  e[1][1] = 0.0;
  e[0][2] = 0.0;
  e[1][2] = (dx/dt);
  e[0][5] = (dx/dt);
  e[1][5] = (dx/dt);
  e[0][6] = -(dx/dt);
  e[1][6] = (dx/dt);

  for (int k = 0; k <= 1; k++) {
    e[k][3] = -1.0*e[k][1];
    e[k][4] = -1.0*e[k][2];
    e[k][7] = -1.0*e[k][5];
    e[k][8] = -1.0*e[k][6];
  }
}
void set_e();

inline void sbound(double (&A)[LX][LY]) { // boundary condition for scalar
  for(int i = 1; i <= LX-2; i++){
    // // // neumann
    A[i][0]      = A[i][1];
    A[i][LY-1]   = A[i][LY-2];
    // // phototactic neumann
    // double al;
    // al = v0ph/b;
    // A[i][0]      = (2.0 -al)*A[i][1]/(2.0 +al);
    // A[i][LY-1]   = (2.0 +al)*A[i][LY-2]/(2.0 -al);
  }
  for(int i = 1; i <= LY-2; i++){
    // neumann
    // A[0][i]      = A[1][i];
    // A[LX-1][i]   = A[LX-2][i];
    // periodic
    A[0][i]      = A[LX-2][i];
    A[LX-1][i]   = A[1][i];
  }
}
void sbound();

void vbound(double (&V)[2][LX][LY]) {
  for(int k = 0; k < 2; k++){
    for(int i = 1; i <= LX-2; i++){
      // // neumann
      V[k][i][0]     = -V[k][i][1];
      V[k][i][LY-1]  = -V[k][i][LY-2];
    }
    for(int i = 1; i <= LY-2; i++){
      // neumann
      // V[k][0][i]     = -V[k][1][i];
      // V[k][LX-1][i]  = -V[k][LX-2][i];
      // // x-periodic
      V[k][0][i]     = V[k][LX-2][i];
      V[k][LX-1][i]  = V[k][1][i];
    }
  }
}
void vbound();

inline void slap(double (&p)[LX][LY], double (&lp)[LX][LY]) { //laplacian for scalar
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      lp[i][j] = (p[i+1][j] +p[i-1][j] +p[i][j+1] +p[i][j-1] -4.0*p[i][j])/(dx*dx);
    }
  }
}
void slap();

inline void grad(double (&p)[LX][LY], double (&grp)[2][LX][LY]){ // central diff
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      grp[0][i][j] = (p[i+1][j] -p[i-1][j])/(2.0*dx);
      grp[1][i][j] = (p[i][j+1] -p[i][j-1])/(2.0*dx);
    }
  }
}
void grad();

inline void div(double (&A)[2][LX][LY], double (&B)[LX][LY]){ // central diff
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      B[i][j] = (A[0][i+1][j] -A[0][i-1][j] +A[1][i][j+1] -A[1][i][j-1])/(2.0*dx);
    }
  }
}
void div();

inline void calc_mu(double (&p)[LX][LY], double (&mu)[LX][LY]) {
  // // rule : (1) operand , (2) product , (3) parameters
  double nn[LX][LY];
  slap(p, nn);
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      mu[i][j] = a*p[i][j] +b*p[i][j]*p[i][j]*p[i][j] -kappa*nn[i][j];
    }
  }
}
void calc_mu();

inline void init_zero(double (&A)[LX][LY]) {
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      A[i][j] = 0.0;
    }
  }
}
void init_zero();

inline void init_one(double (&A)[LX][LY]) {
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      A[i][j] = 1.0;
    }
  }
}
void init_one();

inline void calc_n(double (&f)[9][LX][LY], double (&n)[LX][LY]) {
  // // calc sum of local distribution density f
  double sum_f;
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      sum_f = 0.0;
      for (int k = 0; k < 9; k++) {
        sum_f += f[k][i][j];
      }
      n[i][j] = sum_f;
    }
  }
}
void calc_n();

inline void calc_velocity(double (&f)[9][LX][LY], double (&e)[2][10],
                          double (&nrho)[LX][LY], double (&u)[2][LX][LY]){
  double nu[2];
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      nu[0] = 0.0;
      nu[1] = 0.0;
      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < 9; l++) {
          nu[k] += f[l][i][j]*e[k][l];
        }
        u[k][i][j] = nu[k]/nrho[i][j];
      }
    }
  }
}
void calc_velocity();

inline void calc_forcing(double (&w)[9], double (&e)[2][10], double (&us)[2][LX][LY],
                            double (&force)[2][LX][LY] ,double (&forcing)[9][LX][LY]) {
  double cs2 = dx*dx/(dt*dt*3.0);
  double cs4 = dx*dx*dx*dx/(dt*dt*dt*dt*9.0);
  double v2[2], v4[2];
  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      for (int k = 0; k < 9; k++) {
        v2[0] = (e[0][k] -us[0][i][j])/cs2;
        v2[1] = (e[1][k] -us[1][i][j])/cs2;
        v4[0] = (e[0][k]*us[0][i][j] +e[1][k]*us[1][i][j])*e[0][k]/cs4;
        v4[1] = (e[0][k]*us[0][i][j] +e[1][k]*us[1][i][j])*e[1][k]/cs4;
        forcing[k][i][j] = (1.0 -dt*0.5/tau)*w[k]
                            *((v2[0] +v4[0])*force[0][i][j] +(v2[1] +v4[1])*force[1][i][j]);
      }
    }
  }
}
void calc_forcing();
