/*Riemann*/
#include <math.h>
#include <stdio.h>
#define xMAX 20
#define tMAX 10
int main(void) {
  double K,nu,gamma,rho1new,rho2new,u1new,u2new;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1];
  nu = 0.4; /*Courant*/
  K = 1;
  gamma = 1.4;
  /*Initial Condition HEAD*/
  for (int j = 0; j <= xMAX; ++j) {
    if (j<=xMAX/2-1) {
      rho[j] = 1;
      u[j] = 0;
      p[j] = 1;
    }else{
      rho[j] = 0.125;
      u[j] = 0;
      p[j] = 0.1;
    }
  }
  /*Initial Condition END*/
  FILE *file;
  file = fopen("shock-test1-pressure.dat","w");
  fprintf(file, "%lf ",0.0);
  printf("縦軸Timestep・横軸各Grid ");
  for (int j = 0; j <= xMAX; j++) {
    fprintf(file,"%lf ",deltax * j);
    printf("%lf ",deltax * j);
  }
  fprintf(file, "\n");
  printf("\n");
  fprintf(file,"%d ",0);
  printf("%d ",0);
  for (int j = 0; j <= xMAX; j++) {
    fprintf(file, "%lf ",p[j]);
    printf("%lf ",p[j]);
  }
  fprintf(file, "\n");
  printf("\n");
  /*Time Evolution HEAD*/
  for (int n = 1; n <= tMAX; n++) {
    time = n * nu * deltax;
    for (int j = 1; j < xMAX; j++) {
      rho1new = (rho[j] + rho[j-1]) / 2 - nu * (rho[j] * u[j] - rho[j-1] * u[j-1]) / 2;
      rho2new = (rho[j+1] + rho[j]) / 2 - nu * (rho[j+1] * u[j+1] - rho[j] * u[j]) / 2;
      p[j] = K * pow(rho[j],gamma);
      u1new = ((rho[j] * u[j] + rho[j-1] * u[j-1]) / 2 - nu * (rho[j] * u[j] * u[j] + p[j] - rho[j-1] * u[j-1] * u[j-1] - p[j-1]) / 2) / rho1new;
      u2new = ((rho[j+1] * u[j+1] + rho[j] * u[j]) / 2 - nu * (rho[j+1] * u[j+1] * u[j+1] + p[j+1] - rho[j] * u[j] * u[j] - p[j]) / 2) / rho2new;
      rhonew[j] = rho[j] - nu * (rho2new * u2new - rho1new * u1new);
      p1new = K * pow(rho1new,gamma);
      p2new = K * pow(rho2new,gamma);
      unew[j] = (rho[j] * u[j] - nu * (rho2new * u2new * u2new + p2new - rho1new * u1new * u1new - p1new)) / rhonew[j];
    }
    /*Boundary Condition HEAD*/
    rhonew[0] = rhonew[1];
    unew[0] = 0;
    rhonew[xMAX] = rhonew[xMAX-1];
    unew[xMAX] = 0;
    /*Boundary Condition END*/
    /*Update HEAD*/
    for (int j = 0; j <= xMAX; j++) {
      rho[j] = rhonew[j];
      u[j] = unew[j];
    }
    /*Update END*/
    /*Out Put HEAD*/
    fprintf(file, "%lf ",time);
    printf("time=%lf ", time);
    for (int j = 0; j <= xMAX; j++) {
      fprintf(file, "%lf ",p[j]);
      printf("%lf ",p[j]);
    }
    fprintf(file,"\n");
    printf("\n");
    /*Out Put END*/
  }
  /*Time Evolution END*/
  fclose(file);
}
