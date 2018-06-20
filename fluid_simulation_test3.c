#include <math.h>
#include <stdio.h>
#define xMAX 15
#define tMAX 10
int main(void) {
  FILE *file;
	file = fopen("num-test3.dat","w");
  int x0 = 6;/*boundary value*/
  double x,y0,nu;
  double ynew[xMAX+1],y[xMAX+1],yold[xMAX+1];
  nu = 0.5;
  /*First step HEAD*/
  fprintf(file, "%lf ",0.0);
  printf("縦軸Timestep・横軸各Grid ");
  for (int n = 0; n < xMAX; n++) {
    if (n <= x0) {
      yold[n] = 1;
      y[n] = 1;
    }else{
      yold[n] = 0;
      y[n] = 0;
    }
    fprintf(file,"%d ",n);
    printf("%d ",n);
  }
  /*First step END*/
  fprintf(file, "\n");
  printf("\n");
  fprintf(file,"%d ",0);
  printf("%d ",0);
  for (int n = 0; n < xMAX; n++) {
    fprintf(file, "%lf ",y[n]);
    printf("%lf ",y[n]);
  }
  fprintf(file, "\n");
  printf("\n");
  y[xMAX]=0;
  /*Time Evolution HEAD*/
  for (int j = 1; j <= tMAX; j++) {
    fprintf(file, "%d ",j);
    printf("%d ", j);
    for (int n = 1; n < xMAX; n++) {
      ynew[n] = yold[n] - nu * (y[n+1] - y[n-1]);
    }
    /*Update HEAD*/
    for (int n = 1; n < xMAX; n++) {
      if (ynew[n]<0) {
        ynew[n] = 0;
      }
      yold[n] = y[n];
      y[n] = ynew[n];
    }
    /*Update END*/
    for (int n = 0; n < xMAX; n++) {
      fprintf(file, "%lf ",y[n]);
      printf("%lf ",y[n]);
    }
    fprintf(file,"\n");
    printf("\n");
  }
  /*Time Evolution END*/
  fclose(file);
}
