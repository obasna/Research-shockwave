/*Riemann*/
/*Include Energy Equation*/
#include <math.h>
#include <stdio.h>
#define xMAX 1000
#define tMAX 1
int main(void) {
  double eta,f1new,f2new,delta0,delta1,delta2,delta3,deltah1,deltah2,u3new[xMAX+1],p3new[xMAX+1],uhnew[xMAX+1],phnew[xMAX+1];
  double c,Tout,Tin,K,R,nu,gamma,deltax,deltat,time,rho2new,rho1new,u1new,u2new,p1new,p2new;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1],T[xMAX+1];
  eta = 0.5;
  nu = 0.4; /*Courant*/
  K = 1;
  R = 8.31451e7; /*気体定数*/
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
	file = fopen("shock-test-rho-all.dat","w");
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
      /*if (j==xMAX/2-1) {
        printf("前ここに注目rho1new   %lf\n",rho1new);
        printf("前ここに注目rho2new   %lf\n",rho2new);
        printf("前ここに注目u1new     %lf\n",u1new);
        printf("前ここに注目u2new     %lf\n",u2new);
        printf("前ここに注目rhonew[j] %lf\n",rhonew[j]);
        printf("前ここに注目p1new     %lf\n",p1new);
        printf("前ここに注目p2new     %lf\n",p2new);
        printf("前ここに注目unew[j]   %lf\n",unew[j]);
        printf("前ここに注目pnew[j]   %lf\n",pnew[j]);
        printf("前ここに注目p[j]      %lf\n",p[j]);
      }*/
      rho1new = (rho[j] + rho[j-1]) / 2 - nu * (rho[j] * u[j] - rho[j-1] * u[j-1]) / 2;
      rho2new = (rho[j+1] + rho[j]) / 2 - nu * (rho[j+1] * u[j+1] - rho[j] * u[j]) / 2;
      u1new = (rho[j] * u[j] + rho[j-1] * u[j-1]  - nu * (rho[j] * u[j] * u[j] + p[j] - rho[j-1] * u[j-1] * u[j-1] - p[j-1])) / 2 / rho1new;
      u2new = (rho[j+1] * u[j+1] + rho[j] * u[j]  - nu * (rho[j+1] * u[j+1] * u[j+1] + p[j+1] - rho[j] * u[j] * u[j] - p[j])) / 2 / rho2new;
      p1new = (gamma - 1) / 2 * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 + p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2 - nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1]) - rho1new * u1new * u1new);
      p2new = (gamma - 1) / 2 * (p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2 + p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - nu * ((gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1] - (gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j]) - rho2new * u2new * u2new);
      rhonew[j] = rho[j] - nu * (rho2new * u2new - rho1new * u1new);
      u3new[j] = (rho[j] * u[j] - nu * (rho2new * u2new * u2new + p2new - rho1new * u1new * u1new - p1new)) / rhonew[j];
      p3new[j] = (gamma - 1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - nu * ((gamma / (gamma - 1) * p2new + rho2new * u2new * u2new / 2) * u2new - (gamma / (gamma - 1) * p1new + rho1new * u1new * u1new / 2) * u1new) - rhonew[j] * u3new[j] * u3new[j] / 2);
      uhnew[j] = u3new[j] + eta * (u[j+1] - 2 * u[j] + u[j-1]);
      phnew[j] = p3new[j] + eta * (p[j+1] - 2 * p[j] + p[j-1]);
    if (j==xMAX/2-2||j==xMAX/2-1||j==xMAX/2) {
        printf("後ここに注目rho[j]    %lf\n",rho[j]);
        printf("後ここに注目rho[j+1]  %lf\n",rho[j+1]);
        printf("後ここに注目rho1new   %lf\n",rho1new);
        printf("後ここに注目rho2new   %lf\n",rho2new);
        printf("後ここに注目u[j]      %lf\n",u[j]);
        printf("後ここに注目u[j+1]    %lf\n",u[j+1]);
        printf("後ここに注目u1new     %lf\n",u1new);
        printf("後ここに注目u2new     %lf\n",u2new);
        printf("後ここに注目rhonew[j] %lf\n",rhonew[j]);
        printf("後ここに注目p[j]      %lf\n",p[j]);
        printf("後ここに注目p[j+1]    %lf\n",p[j+1]);
        printf("後ここに注目p1new     %lf\n",p1new);
        printf("後ここに注目p2new     %lf\n",p2new);
        printf("後ここに注目u3new[j]   %lf\n",u3new[j]);
        printf("後ここに注目p3new[j]   %lf\n",p3new[j]);
        printf("後ここに注目uhnew[j]   %lf\n",uhnew[j]);
        printf("後ここに注目phnew[j]   %lf\n",phnew[j]);
        printf("\n");
      }
    }
    /*Last Step for solving velocity HEAD*/
    for (int j = 1; j < xMAX; j++) {
      delta0 = uhnew[j-1] - uhnew[j-2];
      delta1 = uhnew[j] - uhnew[j-1];
      delta2 = uhnew[j+1] - uhnew[j];
      delta3 = uhnew[j+2] - uhnew[j+1];
      deltah1 = eta * (u3new[j] - u3new[j-1]);
      deltah2 = eta * (u3new[j+1] - u3new[j]);
      if (deltah1 > 0) {
        f1new = fmax(0,fmin(fmin(delta0,fabs(delta1)),fmin(fabs(delta1),delta2)));
      }else if (delta1 == 0) {
        f1new = 0;
      }else{
        f1new = -fmax(0,fmin(fmin((-1)*delta0,fabs(delta1)),fmin(fabs(delta1),(-1)*delta2)));
      }
      if (deltah2 > 0) {
        f2new = fmax(0,fmin(fmin(delta1,fabs(delta2)),fmin(fabs(delta2),delta3)));
      }else if (delta2 == 0) {
        f2new = 0;
      }else{
        f2new = -fmax(0,fmin(fmin((-1)*delta1,fabs(delta2)),fmin(fabs(delta2),(-1)*delta3)));
      }
      unew[j] = uhnew[j] - (f2new - f1new);
    }
    /*Last Step for solving velocity END*/
    /*Last Step for solving velocity HEAD*/
    for (int j = 1; j < xMAX; j++) {
      delta0 = phnew[j-1] - phnew[j-2];
      delta1 = phnew[j] - phnew[j-1];
      delta2 = phnew[j+1] - phnew[j];
      delta3 = phnew[j+2] - phnew[j+1];
      deltah1 = eta * (p3new[j] - p3new[j-1]);
      deltah2 = eta * (p3new[j+1] - p3new[j]);
      if (deltah1 > 0) {
        f1new = fmax(0,fmin(fmin(delta0,fabs(delta1)),fmin(fabs(delta1),delta2)));
      }else if (delta1 == 0) {
        f1new = 0;
      }else{
        f1new = -fmax(0,fmin(fmin((-1)*delta0,fabs(delta1)),fmin(fabs(delta1),(-1)*delta2)));
      }
      if (deltah2 > 0) {
        f2new = fmax(0,fmin(fmin(delta1,fabs(delta2)),fmin(fabs(delta2),delta3)));
      }else if (delta2 == 0) {
        f2new = 0;
      }else{
        f2new = -fmax(0,fmin(fmin((-1)*delta1,fabs(delta2)),fmin(fabs(delta2),(-1)*delta3)));
      }
      pnew[j] = phnew[j] - (f2new - f1new);
    }
    /*Last Step for solving velocity END*/
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
