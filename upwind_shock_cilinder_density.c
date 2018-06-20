/*Riemann*/
#include <math.h>
#include <stdio.h>
#define xMAX 1000
#define tMAX 500
int main(void) {
  double K,nu,gamma,deltax,time;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1];
  nu = 0.4; /*Courant*/
  K = 1;
  gamma = 1.4;
  deltax = 1 / (double)xMAX;
  printf("%lf\n",deltax);
  /*Initial Condition HEAD*/
  for (int j = 0; j <= xMAX; ++j) {
    if (j <= xMAX / 2 - 1) {
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
  /*Term of Density HEAD*/
  FILE *file;
	file = fopen("upwind-shock-test1-rho.dat","w");
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
    fprintf(file, "%lf ",rho[j]);
    printf("%lf ",rho[j]);
  }
  fprintf(file, "\n");
  printf("\n");
  /*Time Evolution HEAD*/
  for (int n = 1; n <= tMAX; n++) {
    time = n * nu * deltax;
    for (int j = 1; j < xMAX; j++) {
      if (u[j]>0) {
        rhonew[j] = rho[j] - nu * (rho[j] * u[j] - rho[j-1] * u[j-1]);
        unew[j] = (rho[j] * u[j] - nu * (rho[j] * u[j] * u[j] - rho[j-1] * u[j-1] * u[j-1]) - nu / 2 * (p[j+1] - p[j-1])) / rhonew[j];
        pnew[j] = (gamma-1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 - nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1]));
      }else{
        rhonew[j] = rho[j] + nu * (rho[j] * u[j] - rho[j+1] * u[j+1]);
        unew[j] = (rho[j] * u[j] + nu * (rho[j] * u[j] * u[j] - rho[j+1] * u[j+1] * u[j+1]) - nu / 2 * (p[j+1] - p[j-1])) / rhonew[j];
        pnew[j] = (gamma-1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 + nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1]));
      }
      if (j==xMAX/2-2||j==xMAX/2-1||j==xMAX/2) {
          printf("後ここに注目rho[j]    %lf\n",rho[j]);
          printf("後ここに注目rho[j+1]  %lf\n",rho[j+1]);
          printf("後ここに注目u[j]      %lf\n",u[j]);
          printf("後ここに注目u[j+1]    %lf\n",u[j+1]);
          printf("後ここに注目rhonew[j] %lf\n",rhonew[j]);
          printf("後ここに注目p[j]      %lf\n",p[j]);
          printf("後ここに注目p[j+1]    %lf\n",p[j+1]);
          printf("後ここに注目unew[j]   %lf\n",unew[j]);
          printf("後ここに注目pnew[j]   %lf\n",pnew[j]);
          printf("後ここに注目u[0]      %lf\n",u[0]);
          printf("後ここに注目u[xMAX]   %lf\n",u[xMAX]);
        }
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
    }
    /*Update END*/
    /*Out Put HEAD*/
    fprintf(file, "%lf ",time);
    printf("time=%lf ", time);
    for (int j = 0; j <= xMAX; j++) {
      fprintf(file, "%lf ",rho[j]);
      printf("%lf ",rho[j]);
    }
    fprintf(file,"\n");
    printf("\n");
    /*Out Put END*/
  }
  /*Time Evolution END*/
  fclose(file);
}
