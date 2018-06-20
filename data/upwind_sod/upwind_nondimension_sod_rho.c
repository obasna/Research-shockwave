/*Riemann*/
/*Include Energy Equation*/
#include <math.h>
#include <stdio.h>
/*sgnの定義 HEAD*/
double sgn(double A){
    return (A>0)-(A<0);
}
#define xMAX 1000
#define tMAX 2000000
int main(void) {
  int f;
  double c,Tout,Tin,K,R,nu,gamma,deltax,deltat,time,rho2new,rho1new,u1new,u2new,p1new,p2new;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1],T[xMAX+1];
  R = 1; /*気体定数*/
  c = 3.4e4;/*音速*/
  gamma = 1.4;
  nu = 1e-5;
  deltax = 1 / (double)xMAX;
  /*Initial Condition HEAD*/
  for (int j = 0; j <= xMAX; ++j) {
    if (j <= xMAX / 2 - 1) {
      rho[j] = 1;
      u[j] = 0;
      p[j] = 1;
      T[j] = p[j] / rho[j] / R;
    }else{
      rho[j] = 0.125;
      u[j] = 0;
      p[j] = 0.1;
      T[j] = p[j] / rho[j] / R;
      }
  }
  /*Initial Condition END*/
  /*Term of Density HEAD*/
  FILE *file;
	file = fopen("upwind-sod-nondimension-T-grid1000.dat","w");
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
    fprintf(file, "%8.5e ",T[j]);
    printf("%8.5e ",T[j]);
  }
  fprintf(file, "\n");
  printf("\n");
  /*Time Evolution HEAD*/
  for (int n = 1; n <= tMAX; n++) {
    time = n * nu * deltax;
    for (int j = 1; j < xMAX; j++) {
      /*if (j==xMAX/2-2||j==xMAX/2-1||j==xMAX/2) {
        printf("前ここに注目rho[j] %8.5e\n",rho[j]);
        printf("前ここに注目u[j]   %8.5e\n",u[j]);
        printf("前ここに注目p[j]   %8.5e\n",p[j]);
        printf("\n");
      }*/
      /*if (u[j]>0) {
        rhonew[j] = rho[j] - nu * (rho[j] * u[j] - rho[j-1] * u[j-1]);
        unew[j] = (rho[j] * u[j] - nu * (rho[j] * u[j] * u[j] - rho[j-1] * u[j-1] * u[j-1]) - nu / 2 * (p[j+1] - p[j-1])) / rhonew[j];
        pnew[j] = (gamma - 1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 - nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1]));
      }else{
        rhonew[j] = rho[j] + nu * (rho[j] * u[j] - rho[j+1] * u[j+1]);
        unew[j] = (rho[j] * u[j] + nu * (rho[j] * u[j] * u[j] - rho[j+1] * u[j+1] * u[j+1]) - nu / 2 * (p[j+1] - p[j-1])) / rhonew[j];
        pnew[j] = (gamma - 1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 + nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1]));
      }*/
        f = j - sgn(u[j]);
        rhonew[j] = rho[j] - sgn(u[j]) * nu * (rho[j] * u[j] - rho[f] * u[f]);
        unew[j] = (rho[j] * u[j] - sgn(u[j]) * nu * (rho[j] * u[j] * u[j] - rho[f] * u[f] * u[f]) - nu / 2 * (p[j+1] - p[j-1])) / rhonew[j];
        pnew[j] = (gamma - 1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 - sgn(u[j]) * nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[f] + rho[f] * u[f] * u[f] / 2) * u[f]));
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
    if (n%1000==0) {
      fprintf(file, "%8.5e ",time);
      printf("time=%8.5e ", time);
      for (int j = 0; j <= xMAX; j++) {
        fprintf(file, "%8.5e ",T[j]);
        printf("%8.5e ",T[j]);
      }
      fprintf(file,"\n");
      printf("\n");
    }
    /*Out Put END*/
  }
  /*Time Evolution END*/
  fclose(file);
}
