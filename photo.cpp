#include<stdio.h>
#include<stdlib.h>
#include<random>
#include<iostream>
#include<time.h>
#include<cmath>
#include<fstream>
#include <chrono>
#include <omp.h>
#include "myfunc.h"

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> pinit(0.0, 1.0); // 0 ~ 1　一様乱数
double Q[dat_num*LX*LY];// output phi
void output_Q(const char *file, int timestep)
{

   FILE *fp;
   char filename[100];

   // open file
   sprintf(filename,"%s%d",file,timestep);
   if(NULL==(fp=fopen(filename, "w"))){
     fprintf(stderr, "Cannot open file.\n");
     exit(EXIT_FAILURE);
   }

   // output Q
   for(int k = 0; k < dat_num*(LX-2)*(LY-2); k++){
     fwrite(&Q[k],sizeof(double),1,fp);
   }

   // close file
   fclose(fp);

}
double V[dat_num*LX*LY*2]; // for velocity field
void output_V(const char *file, int timestep)
{

   FILE *fp;
   char filename[100];

   // open file
   sprintf(filename,"%s%d",file,timestep);
   if(NULL==(fp=fopen(filename, "w"))){
     fprintf(stderr, "Cannot open file.\n");
     exit(EXIT_FAILURE);
   }

   // output Q
   for(int k = 0; k < dat_num*(LX-2)*(LY-2)*2; k++){
     fwrite(&V[k],sizeof(double),1,fp);
   }

   // close file
   fclose(fp);

}

void detect_nan(double a, int b) {
  if (std::isnan(a) == true) {
    std::cout << "At line " << b << "nan has occured" << '\n';
    exit(1);
  }
}

double f[9][LX][LY]; //f_k at (i,j) local distribution function
double ftemp[9][LX][LY]; // f_k to calc streaming process
double nrho[LX][LY]; // local density
double feq; // local equilibrium distribution
double p[LX][LY]; // order parameter
double p3[LX][LY]; // (phi)^3
double grdp[2][LX][LY]; // grad(phi)
double u[2][LX][LY]; // velocity field
double divu[LX][LY];
double w[9]; // D2Q9
double e[2][10]; // D2Q9 vector // これと1つ上の　w　もstaticにするべきだがやり方がわからない
double mu[LX][LY];
double force[2][LX][LY]; // force term of NS eq
double forcing[9][LX][LY]; // forcing term of Boltzmann eq

double lapp3[LX][LY], lapp[LX][LY], lap2p[LX][LY];

double uph[2][LX][LY]; // phototactic velocity
double divuph[LX][LY];
double j_p[2][LX][LY];

int l;
double nv, mean_phi;

int main(/*int argc, char const *argv[]*/) {
  set_w(w);
  set_e(e);

  init_one(nrho);

  for(int i = 1; i <= LX-2; i++){
    for(int j = 1; j <= LY-2; j++){
      // if (j < LY/2-1 ) {
      //   // こっちが下
      //   p[i][j] = 0.99 + (pinit(mt)-0.5)*0.02;
      // }
      // else{
      //   p[i][j] = -0.99 + (pinit(mt)-0.5)*0.02;
      // }
      p[i][j] = phi_in +(pinit(mt)-0.5)*0.2;

      for (int k = 0; k < 9; k++) {
        f[k][i][j] = nrho[i][j]*w[k];
        // f[k][i][j] = nrho[i][j]*w[k]*phi_in;
      }
    }
  }
  calc_n(f, nrho);

  // //************ main loop　*************** // //
  int i, j, k;
  int in, ip, jn, jp;
  for (int t = 0; t < T; t++) {
    // // std::cout << "step number is " << t << '\n';
    sbound(p); // phi(t=n)
    // calc_mu(p, mu);
    // sbound(mu);
    // // for(i = 1; i <= LX-2; i++){
    // //   for(j = 1; j <= LY-2; j++){
    // //     std::cout << i << "," << j << " = " << mu[i][j] << '\n';
    // //     detect_nan(mu[i][j], 78);
    // //   }
    // // }
    // for(i = 1; i <= LX-2; i++){
    //   for(j = 1; j <= LY-2; j++){
    //     // std::cout << p[i][j] << '\n';
    //     p3[i][j] = p[i][j]*p[i][j]*p[i][j];
    //     // detect_nan(p3[i][j], 84);
    //   }
    // }
    // sbound(p3);
    slap(p, lapp);
    grad(p, grdp);
    vbound(grdp);
    #pragma omp parallel for
    for(i = 1; i <= LX-2; i++){
      for(j = 1; j <= LY-2; j++){
        // force[0][i][j] = -p[i][j]*(mu[i+1][j] -mu[i-1][j])*0.5/dx; // -phi nabla mu
        // force[1][i][j] = -p[i][j]*(mu[i][j+1] -mu[i][j-1])*0.5/dy -grav*p[i][j];
        // force[0][i][j] = -zeta*lapp[i][j]*grdp[0][i][j];
        // force[1][i][j] = -zeta*lapp[i][j]*grdp[1][i][j] -grav*p[i][j];
        force[0][i][j] = 0.0;
        force[1][i][j] = -grav*p[i][j];
        // force[1][i][j] = -grav*(1.0 +p[i][j]);
        // force[1][i][j] = -(0.001 +grav*p[i][j]);
      }
    }

    calc_n(f, nrho);
    calc_velocity(f, e, nrho, u);
    #pragma omp parallel for
    for(i = 1; i <= LX-2; i++){
      for(j = 1; j <= LY-2; j++){
        for (k = 0; k < 2; k++) {
          u[k][i][j] += 0.5*dt*force[k][i][j]/nrho[i][j];
        }
      }
    }

    // //************* update order parameter ******************  //
    // grad(p,grdp);
    vbound(u);
    div(u, divu);

    // phototactic term
    #pragma omp parallel for
    for(i = 1; i <= LX-2; i++){
      for(j = 1; j <= LY-2; j++){
        uph[0][i][j] = 0.0;
        // uph[1][i][j] = v0ph*(((double) j)*std::exp(-((double) j)/lmd)/lmd );
        // uph[1][i][j] = v0ph*std::exp(-((double) j)/lmd)*lmd;
        uph[1][i][j] = v0ph;
      }
    }
    vbound(uph);
    div(uph, divuph);
    // #pragma omp parallel for
    // for(i = 1; i <= LX-2; i++){
    //   for(j = 1; j <= LY-2; j++){
    //     j_p[0][i][j] = u[0][i][j]*p[i][j] -b*grdp[0][i][j];
    //     j_p[1][i][j] = u[1][i][j]*p[i][j] +uph[1][i][j]*p[i][j] -b*grdp[1][i][j];
    //   }
    // }
    // vbound(j_p);
    #pragma omp parallel for
    for(i = 1; i <= LX-2; i++){
      for(j = 1; j <= LY-2; j++){
        p[i][j] += -dt*(p[i][j]*divu[i][j] +u[0][i][j]*grdp[0][i][j] +u[1][i][j]*grdp[1][i][j]
                     +p[i][j]*divuph[i][j] +uph[0][i][j]*grdp[0][i][j] +uph[1][i][j]*grdp[1][i][j] );
        // p[i][j] += -dt*(p[i][j]*divu[i][j] +u[0][i][j]*grdp[0][i][j] +u[1][i][j]*grdp[1][i][j]
        //              +p[i][j]*divuph[i][j] +uph[0][i][j]*grdp[0][i][j] +uph[1][i][j]*grdp[1][i][j] )
        //             +dt*b*lapp[i][j];
        // p[i][j] += -dt*0.5*(j_p[0][i+1][j] -j_p[0][i-1][j] +j_p[1][i][j+1] -j_p[1][i][j-1]);
      }
    }

    // lap(p^3), lap(p(t=n+1/2), lap(lap(p(t=n+1/2)))), respectively
    // slap(p3, lapp3);
    sbound(p);
    slap(p, lapp);
    // sbound(lapp);
    // slap(lapp, lap2p);
    for(i = 1; i <= LX-2; i++){
      for(j = 1; j <= LY-2; j++){
        // p[i][j] += dt*Gamma*(a*lapp[i][j] +b*lapp3[i][j] -kappa*lap2p[i][j]);
        p[i][j] += dt*b*lapp[i][j];
      }
    }

    // //*********** update distribution function *************// //

    // //********* collision and force
    // for (i = 0; i < 9; i++) {
    //   init_zero(forcing[i]);
    // }
    calc_forcing(w, e, u, force, forcing);

    //************ Collision --> Streaming ********// //
    // Colliding
    #pragma omp parallel for
    for(i = 1; i <= LX-2; i++){
      for(j = 1; j <= LY-2; j++){
        for(k = 0; k < 9; k++) {
          // calc_MBd(u[0][i][j], u[1][i][j], nrho[i][j], e[0][k], e[1][k], w[k], feq);
          feq = nrho[i][j]*w[k]*( 1.0 +3.0*(e[0][k]*u[0][i][j] +e[1][k]*u[1][i][j])/(csq)
                +4.5*(e[0][k]*u[0][i][j] +e[1][k]*u[1][i][j])*(e[0][k]*u[0][i][j] +e[1][k]*u[1][i][j])/(csq*csq)
                -1.5*(u[0][i][j]*u[0][i][j] +u[1][i][j]*u[1][i][j])/(csq));
          // detect_nan(feq, 177);
          // if (std::isnan(feq) == true) {
          //   std::cout << "At line 223, " << i << " , " << j << " ; " << k << " nan has occured" << '\n';
          //   std::cout << "n = " << nrho[i][j] << " , w[k] = " << w[k] << " , e = " << e[0][k] << " , "<< e[1][k]<< '\n';
          //   std::cout << "u = " << u[0][i][j] << " , " << u[1][i][j] << '\n';
          //   exit(1);
          // }
          ftemp[k][i][j] = f[k][i][j] -(1.0/tau)*(f[k][i][j] -feq) +dt*forcing[k][i][j];
          // ftemp[k][i][j] = f[k][i][j] -(1.0/tau)*(f[k][i][j] -feq);
          detect_nan(f[k][i][j], 265);
          // std::cout << f[k][i][j] << '\n';
        }
      }
    }
    // // // Streaming for Periodic BC
    // for(i = 1; i <= LX-2; i++){
    //   in = (i > 1)?(i -1):(LX -2); //三項演算子
    //   ip = (i < LX-2)?(i +1):(1);
    //   for(j = 1; j <= LY-2; j++){
    //     jn = (j > 1)?(j -1):(LY -2);
    //     jp = (j < LY-2)?(j +1):(1);
    //     f[0][i][j]   = ftemp[0][i][j];
    //     f[1][ip][j]  = ftemp[1][i][j];
    //     f[2][i][jp]  = ftemp[2][i][j];
    //     f[3][in][j]  = ftemp[3][i][j];
    //     f[4][i][jn]  = ftemp[4][i][j];
    //     f[5][ip][jp] = ftemp[5][i][j];
    //     f[6][in][jp] = ftemp[6][i][j];
    //     f[7][in][jn] = ftemp[7][i][j];
    //     f[8][ip][jn] = ftemp[8][i][j];
    //   }
    // }

    // Streaming for Bounceback
    // bulk
    // int in, ip, jn, jp;
    #pragma omp parallel for
    for(i = 1; i <= LX-2; i++){
      in = i -1;
      ip = i +1;
      for(j = 1; j <= LY-2; j++){
        jn = j -1;
        jp = j +1;
        f[0][i][j]   = ftemp[0][i][j];
        f[1][ip][j]  = ftemp[1][i][j];
        f[2][i][jp]  = ftemp[2][i][j];
        f[3][in][j]  = ftemp[3][i][j];
        f[4][i][jn]  = ftemp[4][i][j];
        f[5][ip][jp] = ftemp[5][i][j];
        f[6][in][jp] = ftemp[6][i][j];
        f[7][in][jn] = ftemp[7][i][j];
        f[8][ip][jn] = ftemp[8][i][j];
        // f[0][i][j]     = ftemp[0][i][j];
        // f[1][i+1][j]   = ftemp[1][i][j];
        // f[2][i][j+1]   = ftemp[2][i][j];
        // f[3][i-1][j]   = ftemp[3][i][j];
        // f[4][i][j-1]   = ftemp[4][i][j];
        // f[5][i+1][j+1] = ftemp[5][i][j];
        // f[6][i-1][j+1] = ftemp[6][i][j];
        // f[7][i-1][j-1] = ftemp[7][i][j];
        // f[8][i+1][j-1] = ftemp[8][i][j];
      }
    }
    /// // Streaming near the boundary
    // // Bounceback BC
    for (i = 1; i <= LX-2; i++) {
      // // north
      f[4][i][1]    = ftemp[2][i][1];
      f[7][i][1]    = ftemp[5][i][1];
      f[8][i][1]    = ftemp[6][i][1];
      // // south
      f[2][i][LY-2] = ftemp[4][i][LY-2];
      f[5][i][LY-2] = ftemp[7][i][LY-2];
      f[6][i][LY-2] = ftemp[8][i][LY-2];
    }
    for (j = 1; j <= LY-2; j++) {
      // // // east
      // f[3][LX-2][j] = ftemp[1][LX-2][j];
      // f[6][LX-2][j] = ftemp[8][LX-2][j];
      // f[7][LX-2][j] = ftemp[5][LX-2][j];
      // // // west
      // f[1][1][j]    = ftemp[3][1][j];
      // f[5][1][j]    = ftemp[7][1][j];
      // f[8][1][j]    = ftemp[6][1][j];

      // x-direction : Periodic
      // // east
      f[1][1][j] = ftemp[1][LX-2][j];
      f[8][1][j] = ftemp[8][LX-2][j];
      f[5][1][j] = ftemp[5][LX-2][j];
      // // west
      f[3][LX-2][j] = ftemp[3][1][j];
      f[7][LX-2][j] = ftemp[7][1][j];
      f[6][LX-2][j] = ftemp[6][1][j];
    }

    // //************ Streaming --> Collision ********// //
    // // // Periodic BC
    // for(i = 1; i <= LX-2; i++){
    //   in = (i > 1)?(i -1):(LX -2);
    //   ip = (i < LX-2)?(i +1):(1);
    //   for(j = 1; j <= LY-2; j++){
    //     jn = (j > 1)?(j -1):(LY -2);
    //     jp = (j < LY-2)?(j +1):(1);
    //     ftemp[0][i][j]   = f[0][i][j];
    //     ftemp[1][ip][j]  = f[1][i][j];
    //     ftemp[2][i][jp]  = f[2][i][j];
    //     ftemp[3][in][j]  = f[3][i][j];
    //     ftemp[4][i][jn]  = f[4][i][j];
    //     ftemp[5][ip][jp] = f[5][i][j];
    //     ftemp[6][in][jp] = f[6][i][j];
    //     ftemp[7][in][jn] = f[7][i][j];
    //     ftemp[8][ip][jn] = f[8][i][j];
    //   }
    // }
    //
    // calc_n(ftemp, nrho);
    // // // Colliding
    // for(i = 1; i <= LX-2; i++){
    //   for(j = 1; j <= LY-2; j++){
    //     for(k = 0; k < 9; k++) {
    //       // calc_MBd(u[0][i][j], u[1][i][j], nrho[i][j], e[0][k], e[1][k], w[k], feq);
    //       feq = nrho[i][j]*w[k]*( 1.0 +3.0*(e[0][k]*u[0][i][j] +e[1][k]*u[1][i][j])/(csq)
    //             +4.5*(e[0][k]*u[0][i][j] +e[1][k]*u[1][i][j])*(e[0][k]*u[0][i][j] +e[1][k]*u[1][i][j])/(csq*csq)
    //             -1.5*(u[0][i][j]*u[0][i][j] +u[1][i][j]*u[0][i][j])/(csq));
    //       // detect_nan(feq, 177);
    //       // if (std::isnan(feq) == true) {
    //       //   std::cout << "At line 223, " << i << " , " << j << " ; " << k << " nan has occured" << '\n';
    //       //   std::cout << "n = " << nrho[i][j] << " , w[k] = " << w[k] << " , e = " << e[0][k] << " , "<< e[1][k]<< '\n';
    //       //   std::cout << "u = " << u[0][i][j] << " , " << u[1][i][j] << '\n';
    //       //   exit(1);
    //       // }
    //       f[k][i][j] = ftemp[k][i][j] -(1.0/tau)*(ftemp[k][i][j] -feq) +dt*forcing[k][i][j];
    //       detect_nan(f[k][i][j], 1243);
    //       // std::cout << f[k][i][j] << '\n';
    //     }
    //   }
    // }


    // for(i = 1; i <= LX-2; i++){
    //   for(j = 1; j <= LY-2; j++){
    //     for(k = 0; k < 2; k++) {
    //       detect_nan(f[k][i][j], 250);
    //     }
    //   }
    // }

    // // output
    // if (t%1000==0) {
    //   for (i = 1; i <= LX-2; i++) {
    //     for (j = 1; j <= LY-2; j++) { // i,jの順番に注意
    //       outputfile << p[i][j] << "\t";
    //     }
    //     outputfile << '\n';
    //   }
    // }
    if (t%duration==0) {
      nv = 0.0;
      mean_phi = 0.0;
      l = t/duration;
      for (j = 1; j <= LY-2; j++) {
        for (i = 1; i <= LX-2; i++) {
          Q[l*(LY-2)*(LX-2) +(j-1)*(LX-2) +i-1] = p[i][j];
          // // check velocity
          nv += sqrt(u[0][i][j]*u[0][i][j] +u[1][i][j]*u[1][i][j])/((double) (LX-2)*(LY-2));
          mean_phi += p[i][j]/((double) (LX-2)*(LY-2));
        }
      }
      for (k = 0; k <= 1; k++){
        for (j = 1; j <= LY-2; j++) {
          for (i = 1; i <= LX-2; i++) {
            V[(2*l+k)*(LY-2)*(LX-2) +(j-1)*(LX-2) +i-1] = u[k][i][j];
          }
        }
      }
      std::cout << "progress:" << 100*((double) t)/((double) T)
                << "%" << '\n';
      std::cout << "numerical velocity : " << nv << '\n';
      std::cout << "mean phi : " << mean_phi << '\n';
      // std::cout << "numerical / typical: "
      //           << nv*((tau -0.5)/3)/(sqrt((8.0*a*a*(-a)*kappa)/(9.0*b*b))) << '\n';
      //           std::cout << "mean phi 2 = " << mean_phi << '\n';
    }
   // std::cout << t << '\n';
  }
  // std::cout << "main calc is over !" << '\n';
  // std::ofstream outputfile("rti_128.txt");
  // for(int l = 0; l < dat_num; l++){
  //   for (int i = 1; i <= LX-2; i++) {
  //     for (int j = 1; j <= LY-2; j++) { // i,jの順番に注意
  //       outputfile << outp[l][i][j] << "\t";
  //     }
  //     outputfile << '\n';
  //   }
  // }
  //
  // outputfile.close();
  output_Q("pho_zm05_128_", dat_num);
  output_V("pV_zm05_128_", dat_num);
  std::cout << "successfully computed!" << '\n';
  return 0;
}
