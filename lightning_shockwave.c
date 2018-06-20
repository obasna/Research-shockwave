/*Riemann*/
/*Two step Lax-Wendroff*/
#include <math.h>
#include <stdio.h>
#define xMAX 1000
#define tMAX 100
int main(void) {
  double km;
  double R,Tin,Tout,K,nu,gamma,deltax,time,rho2new,rho1new,u1new,u2new,p1new,p2new;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1],Tnew[xMAX+1],T[xMAX+1];
  R = 8.31451e7; /*気体定数*/
  K = 1;
  gamma = 1.4;
  km = 1e5;
  deltax = 20 * km / (double)xMAX;
  printf("%lf\n",deltax);
  /*Initial Condition HEAD*/
  nu = 5e-7; /*Courant*/
  Tin = 1700; /*Initial Temperature*/
  Tout = 100;
  for (int j = 0; j <= xMAX; ++j) {
    if (j<=xMAX/2-1) {
      rho[j] = 10e-5;
      u[j] = 0;
      p[j] = rho[j] * R * Tin;
    }else{
      rho[j] = 10e-10;
      u[j] = 0;
      p[j] = rho[j] * R * Tout;
      }
  }
  /*Initial Condition END*/
  /*Term of Density HEAD*/
  FILE *file;
	file = fopen("lightning-shock-rho.dat","w");
  fprintf(file, "%lf ",0.0);
  printf("縦軸Timestep・横軸各Grid ");
  for (int j = 0; j <= xMAX; j++) {
    fprintf(file,"%lf ",deltax * j / km);
    printf("%lf ",deltax * j / km);
  }
  fprintf(file, "\n");
  printf("\n");
  fprintf(file,"%d ",0);
  printf("%d ",0);
  for (int j = 0; j <= xMAX; j++) {
    fprintf(file, "%8.5e ",rho[j]);
    printf("%8.5e ",rho[j]);
  }
  fprintf(file, "\n");
  printf("\n");
  /*Time Evolution HEAD*/
  for (int n = 1; n <= tMAX; n++) {
    time = n * nu * deltax;
    for (int j = 1; j < xMAX; j++) {
      rho1new = (rho[j] + rho[j-1]) / 2 - nu * (rho[j] * u[j] - rho[j-1] * u[j-1]) / 2;
      rho2new = (rho[j+1] + rho[j]) / 2 - nu * (rho[j+1] * u[j+1] - rho[j] * u[j]) / 2;
      u1new = (rho[j] * u[j] + rho[j-1] * u[j-1]  - nu * (rho[j] * u[j] * u[j] + p[j] - rho[j-1] * u[j-1] * u[j-1] - p[j-1])) / 2 / rho1new;
      u2new = (rho[j+1] * u[j+1] + rho[j] * u[j]  - nu * (rho[j+1] * u[j+1] * u[j+1] + p[j+1] - rho[j] * u[j] * u[j] - p[j])) / 2 / rho2new;
      p1new = (gamma - 1) / 2 * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 + p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2 - nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1]) - rho1new * u1new * u1new);
      p2new = (gamma - 1) / 2 * (p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2 + p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - nu * ((gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1] - (gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j]) - rho2new * u2new * u2new);
      rhonew[j] = rho[j] - nu * (rho2new * u2new - rho1new * u1new);
      unew[j] = (rho[j] * u[j] - nu * (rho2new * u2new * u2new + p2new - rho1new * u1new * u1new - p1new)) / rhonew[j];
      pnew[j] = (gamma - 1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - nu * ((gamma / (gamma - 1) * p2new + rho2new * u2new * u2new / 2) * u2new - (gamma / (gamma - 1) * p1new + rho1new * u1new * u1new / 2) * u1new) - rhonew[j] * unew[j] * unew[j] / 2);
      Tnew[j] = pnew[j] / rhonew[j] / R;
    }
    /*Boundary Condition HEAD*/
    rhonew[0] = rhonew[1];
    unew[0] = 0;
    pnew[0] = pnew[1];
    rhonew[xMAX] = rhonew[xMAX-1];
    unew[xMAX] = 0;
    pnew[xMAX] = pnew[xMAX-1];
    /*Boundary Condition END*/
    /*Update HEAD*/
    for (int j = 0; j <= xMAX; j++) {
      rho[j] = rhonew[j];
      u[j] = unew[j];
      p[j] = pnew[j];
      T[j] = Tnew[j];
    }
    /*Update END*/
    /*Out Put HEAD*/
    fprintf(file, "%8.5e ",time);
    printf("time=%8.5e ", time);
    for (int j = 0; j <= xMAX; j++) {
      fprintf(file, "%8.5e ",rho[j]);
      printf("%8.5e ",rho[j]);
    }
    fprintf(file,"\n");
    printf("\n");
    /*Out Put END*/
  }
  /*Time Evolution END*/
  fclose(file);
}
