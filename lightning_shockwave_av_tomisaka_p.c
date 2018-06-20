/*Riemann*/
/*Include Energy Equation*/
#include <math.h>
#include <stdio.h>
#define xMAX 1000
#define tMAX 10000000
int main(void) {
  double km;
  double k,q0rho,q1rho,q2rho,q3rho,q4rho,q0u,q1u,q2u,q3u,q4u,q0p,q1p,q2p,q3p,q4p,rho0new,rho1new,rho2new,rho3new,u0new,u1new,u2new,u3new,p0new,p1new,p2new,p3new,q1rhonew,q2rhonew,q1unew,q2unew,q1pnew,q2pnew;
  double c,Tout,Tin,K,R,nu,gamma,deltax,deltat,time;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1],T[xMAX+1];
  k = 0.1;
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
	file = fopen("lightning-sod-p-grid1000-tomisaka.dat","w");
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
      q0rho = k / 8 * (rho[j-1] - 2 * rho[j-2] + rho[j-3]);
      q1rho = k / 8 * (rho[j] - 2 * rho[j-1] + rho[j-2]);
      q2rho = k / 8 * (rho[j+1] - 2 * rho[j] + rho[j-1]);
      q3rho = k / 8 * (rho[j+2] - 2 * rho[j+1] + rho[j]);
      q4rho = k / 8 * (rho[j+3] - 2 * rho[j+2] + rho[j+1]);
      q0u = k / 8 * (rho[j-1] * u[j-1] - 2 * rho[j-2] * u[j-2] + rho[j-3] * u[j-3]);
      q1u = k / 8 * (rho[j] * u[j] - 2 * rho[j-1] * u[j-1] + rho[j-2] * u[j-2]);
      q2u = k / 8 * (rho[j+1] * u[j+1] - 2 * rho[j] * u[j] + rho[j-1] * u[j-1]);
      q3u = k / 8 * (rho[j+2] * u[j+2] - 2 * rho[j+1] * u[j+1] + rho[j] * u[j]);
      q4u = k / 8 * (rho[j+3] * u[j+3] - 2 * rho[j+2] * u[j+2] + rho[j+1] * u[j+1]);
      q0p = k / 8 * (p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2 - 2 * (p[j-2] / (gamma - 1) + rho[j-2] * u[j-2] * u[j-2] / 2) + p[j-3] / (gamma - 1) + rho[j-3] * u[j-3] * u[j-3] / 2);
      q2p = k / 8 * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - 2 * (p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2) + p[j-2] / (gamma - 1) + rho[j-2] * u[j-2] * u[j-2] / 2);
      q2p = k / 8 * (p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2 - 2 * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2) + p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2);
      q3p = k / 8 * (p[j+2] / (gamma - 1) + rho[j+2] * u[j+2] * u[j+2] / 2 - 2 * (p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2) + p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2);
      q4p = k / 8 * (p[j+3] / (gamma - 1) + rho[j+3] * u[j+3] * u[j+3] / 2 - 2 * (p[j+2] / (gamma - 1) + rho[j+2] * u[j+2] * u[j+2] / 2) + p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2);
      rho0new = (rho[j-2] + rho[j-1]) / 2 - nu / 2 * (rho[j-1] * u[j-1] - rho[j-2] * u[j-2] - q1rho + q0rho);
      rho1new = (rho[j-1] + rho[j]) / 2 - nu / 2 * (rho[j] * u[j] - rho[j-1] * u[j-1] - q2rho + q1rho);
      rho2new = (rho[j] + rho[j+1]) / 2 - nu / 2 * (rho[j+1] * u[j+1] - rho[j] * u[j] - q3rho + q2rho);
      rho3new = (rho[j+1] + rho[j+2]) / 2 - nu / 2 * (rho[j+2] * u[j+2] - rho[j+1] * u[j+1] - q4rho + q3rho);
      u0new = ((rho[j-2] * u[j-2] + rho[j-1] * u[j-1]) / 2 - nu / 2 * (rho[j-1] * u[j-1] * u[j-1] + p[j-1] - rho[j-2] * u[j-2] * u[j-2] - p[j-2] - q1u + q0u)) / rho0new;
      u1new = ((rho[j-1] * u[j-1] + rho[j] * u[j]) / 2 - nu / 2 * (rho[j] * u[j] * u[j] + p[j] - rho[j-1] * u[j-1] * u[j-1] - p[j-1] - q2u + q1u)) / rho1new;
      u2new = ((rho[j] * u[j] + rho[j+1] * u[j+1]) / 2 - nu / 2 * (rho[j+1] * u[j+1] * u[j+1] + p[j+1] - rho[j] * u[j] * u[j] - p[j] - q3u + q2u)) / rho2new;
      u3new = ((rho[j+1] * u[j+1] + rho[j+2] * u[j+2]) / 2 - nu / 2 * (rho[j+2] * u[j+2] * u[j+2] + p[j+2] - rho[j+1] * u[j+1] * u[j+1] - p[j+1] - q4u + q3u)) / rho3new;
      p0new = (gamma - 1) * ((p[j-2] / (gamma - 1) + rho[j-2] * u[j-2] * u[j-2] / 2 + p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2) / 2 - rho0new * u0new * u0new / 2 - nu / 2 * ((gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1] - (gamma / (gamma - 1) * p[j-2] + rho[j-2] * u[j-2] * u[j-2] / 2) * u[j-2] - q1p + q0p));
      p1new = (gamma - 1) * ((p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2 + p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2) / 2 - rho1new * u1new * u1new / 2 - nu / 2 * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1] - q2p + q1p));
      p2new = (gamma - 1) * ((p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 + p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2) / 2 - rho2new * u2new * u2new / 2 - nu / 2 * ((gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1] - (gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - q3p + q2p));
      p3new = (gamma - 1) * ((p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2 + p[j+2] / (gamma - 1) + rho[j+2] * u[j+2] * u[j+2] / 2) / 2 - rho3new * u3new * u3new / 2 - nu / 2 * ((gamma / (gamma - 1) * p[j+2] + rho[j+2] * u[j+2] * u[j+2] / 2) * u[j+2] - (gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1] - q4p + q3p));
      q1rhonew = k / 8 * (rho2new - 2 * rho1new + rho0new);
      q2rhonew = k / 8 * (rho3new - 2 * rho2new + rho1new);
      q1unew = k / 8 * (rho2new * u2new * u2new + p2new - 2 * (rho1new * u1new * u1new + p1new) + rho0new * u0new * u0new + p0new);
      q2unew = k / 8 * (rho3new * u3new * u3new + p3new - 2 * (rho2new * u2new * u2new + p2new) + rho1new * u1new * u1new + p1new);
      q1pnew = k / 8 * (p2new / (gamma - 1) + rho2new * u2new * u2new / 2 - 2 * (p1new / (gamma - 1) + rho1new * u1new * u1new / 2) + p0new / (gamma - 1) + rho0new * u0new * u0new / 2);
      q2pnew = k / 8 * (p3new / (gamma - 1) + rho3new * u3new * u3new / 2 - 2 * (p2new / (gamma - 1) + rho2new * u2new * u2new / 2) + p1new / (gamma - 1) + rho1new * u1new * u1new / 2);
      rhonew[j] = rho[j] - nu * (rho2new * u2new - rho1new * u1new - q2rhonew + q1rhonew);
      unew[j] = (rho[j] * u[j] - nu * (rho2new * u2new * u2new + p2new - rho1new * u1new * u1new - p1new - q2unew + q1unew)) / rhonew[j];
      pnew[j] = (gamma - 1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 - nu * ((gamma / (gamma - 1) * p2new + rho2new * u2new * u2new / 2) * u2new - (gamma / (gamma - 1) * p1new + rho1new * u1new * u1new / 2) * u1new - q2pnew + q1pnew));
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
