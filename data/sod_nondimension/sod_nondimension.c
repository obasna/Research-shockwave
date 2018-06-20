/*Riemann*/
/*Two step Lax-Wendroff including the artificial viscosity*/
/*Non-dimension*/
#include <math.h>
#include <stdio.h>
/*sgnの定義 HEAD*/
double sgn(double A){
    return (A>0)-(A<0);
}
/*sgnの定義　END*/
#define xMAX 1000
#define tMAX 20000000
int main(void) {
  double test;
  double rhof1new,rhof2new,rhodelta0,rhodelta1,rhodelta2,rhodelta3,rhodeltah1,rhodeltah2,uf1new,uf2new,udelta0,udelta1,udelta2,udelta3,udeltah1,udeltah2,pf1new,pf2new,pdelta0,pdelta1,pdelta2,pdelta3,pdeltah1,pdeltah2;
  double eta,rho3new[xMAX+1],u3new[xMAX+1],p3new[xMAX+1],rhohnew[xMAX+1],uhnew[xMAX+1],phnew[xMAX+1];
  double c,Tout,Tin,K,R,nu,gamma,deltax,deltat,time,rho2new,rho1new,u1new,u2new,p1new,p2new;
  double rhonew[xMAX+1],rho[xMAX+1],unew[xMAX+1],u[xMAX+1],pnew[xMAX+1],p[xMAX+1],T[xMAX+1];
  eta = 0.1;
  R = 1;
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
	file = fopen("sod-nondimension-grid1000-T-long.dat","w");
  fprintf(file, "%lf ",0.0);
  printf("縦軸Timestep・横軸各Grid ");
  for (int j = 0; j <= xMAX; j++) {
    fprintf(file,"%lf ", deltax * j);
    printf("%lf ", deltax * j);
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
      rho1new = (rho[j] + rho[j-1]) / 2 - nu * (rho[j] * u[j] - rho[j-1] * u[j-1]) / 2;
      rho2new = (rho[j+1] + rho[j]) / 2 - nu * (rho[j+1] * u[j+1] - rho[j] * u[j]) / 2;
      u1new = (rho[j] * u[j] + rho[j-1] * u[j-1]  - nu * (rho[j] * u[j] * u[j] + p[j] - rho[j-1] * u[j-1] * u[j-1] - p[j-1])) / 2 / rho1new;
      u2new = (rho[j+1] * u[j+1] + rho[j] * u[j]  - nu * (rho[j+1] * u[j+1] * u[j+1] + p[j+1] - rho[j] * u[j] * u[j] - p[j])) / 2 / rho2new;
      p1new = (gamma - 1) / 2 * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 + p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2 - nu * ((gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j] - (gamma / (gamma - 1) * p[j-1] + rho[j-1] * u[j-1] * u[j-1] / 2) * u[j-1]) - rho1new * u1new * u1new);
      p2new = (gamma - 1) / 2 * (p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2 + p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - nu * ((gamma / (gamma - 1) * p[j+1] + rho[j+1] * u[j+1] * u[j+1] / 2) * u[j+1] - (gamma / (gamma - 1) * p[j] + rho[j] * u[j] * u[j] / 2) * u[j]) - rho2new * u2new * u2new);
      rho3new[j] = rho[j] - nu * (rho2new * u2new - rho1new * u1new);
      u3new[j] = (rho[j] * u[j] - nu * (rho2new * u2new * u2new + p2new - rho1new * u1new * u1new - p1new)) / rho3new[j];
      p3new[j] = (gamma - 1) * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2 - nu * ((gamma / (gamma - 1) * p2new + rho2new * u2new * u2new / 2) * u2new - (gamma / (gamma - 1) * p1new + rho1new * u1new * u1new / 2) * u1new) - rho3new[j] * u3new[j] * u3new[j] / 2);
      rhohnew[j] = rho3new[j] + eta * (rho[j+1] - 2 * rho[j] + rho[j-1]);
      uhnew[j] = (rho3new[j] * u3new[j] + eta * (rho[j+1] * u[j+1] - 2 * rho[j] * u[j] + rho[j-1] * u[j-1])) / rhohnew[j];
      phnew[j] = (gamma - 1) * (p3new[j] / (gamma - 1) + rho3new[j] * u3new[j] * u3new[j] / 2 - rhohnew[j] * uhnew[j] * uhnew[j] / 2 + eta * (p[j+1] / (gamma - 1) + rho[j+1] * u[j+1] * u[j+1] / 2 - 2 * (p[j] / (gamma - 1) + rho[j] * u[j] * u[j] / 2) + p[j-1] / (gamma - 1) + rho[j-1] * u[j-1] * u[j-1] / 2));
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
        printf("後ここに注目rho3new[j]%8.5e\n",rho3new[j]);
        printf("後ここに注目u3new[j]  %8.5e\n",u3new[j]);
        printf("後ここに注目p3new[j]  %8.5e\n",p3new[j]);
        printf("後ここに注目rhohnew[j]%8.5e\n",rhohnew[j]);
        printf("後ここに注目uhnew[j]  %8.5e\n",uhnew[j]);
        printf("後ここに注目phnew[j]  %8.5e\n",phnew[j]);
        printf("\n");
      }*/
    }
    /*調整 HEAD*/
    rhohnew[0] = rhohnew[1];
    rhohnew[xMAX] = rhohnew[xMAX-1];
    uhnew[0] = uhnew[1];
    uhnew[xMAX] = uhnew[xMAX-1];
    phnew[0] = phnew[1];
    phnew[xMAX] = phnew[xMAX-1];
    /*調整　END*/
    /*Last Step for solving density HEAD*/
    for (int j = 1; j < xMAX-1; j++) {
      rhodelta0 = rhohnew[j-1] - rhohnew[j-2];
      rhodelta1 = rhohnew[j] - rhohnew[j-1];
      rhodelta2 = rhohnew[j+1] - rhohnew[j];
      rhodelta3 = rhohnew[j+2] - rhohnew[j+1];
      rhodeltah1 = eta * (rho3new[j] - rho3new[j-1]);
      rhodeltah2 = eta * (rho3new[j+1] - rho3new[j]);
      rhof1new = sgn(rhodeltah1) * fmax(0,fmin(fmin(sgn(rhodeltah1)*rhodelta0,fabs(rhodeltah1)),fmin(fabs(rhodeltah1),sgn(rhodeltah1) * rhodelta2)));
      rhof2new = sgn(rhodeltah2) * fmax(0,fmin(fmin(sgn(rhodeltah2)*rhodelta1,fabs(rhodeltah2)),fmin(fabs(rhodeltah2),sgn(rhodeltah2) * rhodelta3)));
      rhonew[j] = rhohnew[j] - (rhof2new - rhof1new);
      udelta0 = rhohnew[j-1] * uhnew[j-1] - rhohnew[j-2] * uhnew[j-2];
      udelta1 = rhohnew[j] * uhnew[j] - rhohnew[j-1] * uhnew[j-1];
      udelta2 = rhohnew[j+1] * uhnew[j+1] - rhohnew[j] * uhnew[j];
      udelta3 = rhohnew[j+2] * uhnew[j+2] - rhohnew[j+1] * uhnew[j+1];
      udeltah1 = eta * (rho3new[j] * u3new[j] - rho3new[j-1] * u3new[j-1]);
      udeltah2 = eta * (rho3new[j+1] * u3new[j+1] - rho3new[j] * u3new[j]);
      uf1new = sgn(udeltah1) * fmax(0,fmin(fmin(sgn(udeltah1) * udelta0,fabs(udeltah1)),fmin(fabs(udeltah1),sgn(udeltah1) * udelta2)));
      uf2new = sgn(udeltah2) * fmax(0,fmin(fmin(sgn(udeltah2) * udelta1,fabs(udeltah2)),fmin(fabs(udeltah2),sgn(udeltah2) * udelta3)));
      unew[j] = (rhohnew[j] * uhnew[j] - (uf2new - uf1new)) / rhonew[j];
      pdelta0 = phnew[j-1] / (gamma - 1) + rhohnew[j-1] * uhnew[j-1] * uhnew[j-1] / 2 - phnew[j-2] / (gamma - 1) - rhohnew[j-2] * uhnew[j-2] * uhnew[j-2] / 2;
      pdelta1 = phnew[j] / (gamma - 1) + rhohnew[j] * uhnew[j] * uhnew[j] / 2 - phnew[j-1] / (gamma - 1) - rhohnew[j-1] * uhnew[j-1] * uhnew[j-1] / 2;
      pdelta2 = phnew[j+1] / (gamma - 1) + rhohnew[j+1] * uhnew[j+1] * uhnew[j+1] / 2 - phnew[j] / (gamma - 1) - rhohnew[j] * uhnew[j] * uhnew[j] / 2;
      pdelta3 = phnew[j+2] / (gamma - 1) + rhohnew[j+2] * uhnew[j+2] * uhnew[j+2] / 2 - phnew[j+1] / (gamma - 1) - rhohnew[j+1] * uhnew[j+1] * uhnew[j+1] / 2;
      pdeltah1 = eta * (p3new[j] / (gamma - 1) + rho3new[j] * u3new[j] * u3new[j] / 2 - p3new[j-1] / (gamma - 1) - rho3new[j-1] * u3new[j-1] * u3new[j-1] / 2);
      pdeltah2 = eta * (p3new[j+1] / (gamma - 1) + rho3new[j+1] * u3new[j+1] * u3new[j+1] / 2 - p3new[j] / (gamma - 1) - rho3new[j] * u3new[j] * u3new[j] / 2);
      pf1new = sgn(pdeltah1) * fmax(0,fmin(fmin(sgn(pdeltah1) * pdelta0,fabs(pdeltah1)),fmin(fabs(pdeltah1),sgn(pdeltah1) * pdelta2)));
      pf2new = sgn(pdeltah2) * fmax(0,fmin(fmin(sgn(pdeltah2) * pdelta1,fabs(pdeltah2)),fmin(fabs(pdeltah2),sgn(pdeltah2) * pdelta3)));
      pnew[j] = (gamma - 1) * (phnew[j] / (gamma - 1) + rhohnew[j] * uhnew[j] * uhnew[j] / 2 - rhonew[j] * unew[j] * unew[j] / 2 - (pf2new - pf1new));
    }
    /*Last Step for solving velocity END*/
    /*Boundary Condition HEAD*/
    rhonew[0] = rhonew[1];
    unew[0] = 0;
    pnew[0] = pnew[1];
    rhonew[xMAX-1] = rhonew[xMAX-2];
    unew[xMAX-1] = 0;
    pnew[xMAX-1] = pnew[xMAX-2];
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
