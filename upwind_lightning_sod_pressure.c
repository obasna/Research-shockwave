/*Riemann*/
/*Include Energy Equation*/
#include <math.h>
#include <stdio.h>
#define xMAX 1000
#define tMAX 10000000
int main(void) {
  double km;
  double c,Tout,Tin,K,R,nu,gamma,deltax,deltat,time,rho2new,rho1new,u1new,u2new,p1new,p2new;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1],T[xMAX+1];
  R = 8.31451e7; /*気体定数*/
  c = 3.4e4;/*音速*/
  gamma = 1.4;
  km = 1e5;
  deltat = 1e-7;/*CFL条件では*/
  deltax = 20 * km / (double)xMAX;
  printf("%lf\n",c*deltat/deltax);
  /*Initial Condition HEAD*/
  nu = /*c * */deltat / deltax; /*Courant*/
  Tin = 2000; /*Initial Temperature*/
  Tout = 1600;
  for (int j = 0; j <= xMAX; ++j) {
    if (j<=xMAX/2-1) {
      rho[j] = 1e-5;
      u[j] = 0;
      p[j] = rho[j] * R * Tin;
      T[j] = Tin;
    }else{
      rho[j] = 1.25e-6;
      u[j] = 0;
      p[j] = rho[j] * R * Tout;
      T[j] = Tout;
      }
  }
  /*Initial Condition END*/
  /*Term of Density HEAD*/
  FILE *file;
	file = fopen("upwind-lightning-sod-pressure-grid1000.dat","w");
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
    fprintf(file, "%8.5e ",p[j]);
    printf("%8.5e ",p[j]);
  }
  fprintf(file, "\n");
  printf("\n");
  /*Time Evolution HEAD*/
  for (int n = 1; n <= tMAX; n++) {
    time = deltat * n;
    for (int j = 1; j < xMAX; j++) {
      /*if (j==xMAX/2-2||j==xMAX/2-1||j==xMAX/2) {
        printf("前ここに注目rho[j] %8.5e\n",rho[j]);
        printf("前ここに注目u[j]   %8.5e\n",u[j]);
        printf("前ここに注目p[j]   %8.5e\n",p[j]);
        printf("\n");
      }*/
      if (u[j]>0) {
        rhonew[j] = rho[j] - nu * (rho[j] * u[j] - rho[j-1] * u[j-1]);
        unew[j] = (rho[j] * u[j] - nu * (rho[j] * u[j] * u[j] - rho[j-1] * u[j-1] * u[j-1]) - nu / 2 * (p[j+1] - p[j-1])) / rhonew[j];
        pnew[j] = (gamma-1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 - nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1]));
      }else{
        rhonew[j] = rho[j] + nu * (rho[j] * u[j] - rho[j+1] * u[j+1]);
        unew[j] = (rho[j] * u[j] + nu * (rho[j] * u[j] * u[j] - rho[j+1] * u[j+1] * u[j+1]) - nu / 2 * (p[j+1] - p[j-1])) / rhonew[j];
        pnew[j] = (gamma-1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 + nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1]));
      }
      /*if (j==xMAX/2-2||j==xMAX/2-1||j==xMAX/2) {
        printf("後ここに注目rho[j]    %8.5e\n",rho[j]);
        printf("後ここに注目rho[j+1]  %8.5e\n",rho[j+1]);
        printf("後ここに注目rho1new   %8.5e\n",rho1new);
        printf("後ここに注目rho2new   %8.5e\n",rho2new);
        printf("後ここに注目u[j]      %8.5e\n",u[j]);
        printf("後ここに注目u[j+1]    %8.5e\n",u[j+1]);
        printf("後ここに注目u1new     %8.5e\n",u1new);
        printf("後ここに注目u2new     %8.5e\n",u2new);
        printf("後ここに注目rhonew[j] %8.5e\n",rhonew[j]);
        printf("後ここに注目p[j]      %8.5e\n",p[j]);
        printf("後ここに注目p[j+1]    %8.5e\n",p[j+1]);
        printf("後ここに注目p1new     %8.5e\n",p1new);
        printf("後ここに注目p2new     %8.5e\n",p2new);
        printf("後ここに注目unew[j]   %8.5e\n",unew[j]);
        printf("後ここに注目pnew[j]   %8.5e\n",pnew[j]);
        printf("後ここに注目u[0]      %8.5e\n",u[0]);
        printf("後ここに注目u[xMAX]   %8.5e\n",u[xMAX]);
        printf("\n");
      }*/
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
      T[j] = p[j] / rho[j] / R;
    }
    /*Update END*/
    /*Out Put HEAD*/
    if (n%10000==0) {
      fprintf(file, "%8.5e ",time);
      printf("time=%8.5e ", time);
      for (int j = 0; j <= xMAX; j++) {
        fprintf(file, "%8.5e ",p[j]);
        printf("%8.5e ",p[j]);
      }
      fprintf(file,"\n");
      printf("\n");
    }
    /*Out Put END*/
  }
  /*Time Evolution END*/
  fclose(file);
}
