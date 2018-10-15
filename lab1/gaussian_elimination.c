
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define num double
#define eps 1e-12
#define true 1
#define false 0
int gaussian_elimination_1(int n, num**AB,num*X){
  int i,j,k;
  num m,s;


  for(i = 0; i < n - 1; i++)
  {
    for(j = i + 1; j < n; j++)
    {
      if(fabs(AB[i][i]) < eps) return false;
      m = -AB[j][i] / AB[i][i];
      for(k = i + 1; k <= n; k++)
        AB[j][k] += m * AB[i][k];
    }
  }


  for(i = n - 1; i >= 0; i--)
  {
    s = AB[i][n];
    for(j = n - 1; j >= i + 1; j--)
      s -= AB[i][j] * X[j];
    if(fabs(AB[i][i]) < eps) return false;
    X[i] = s / AB[i][i];
  }
  return true;
}
int main(int argc, char**argv){
    if(argc<2){
        fprintf(_IO_stderr,"Not enough arguments supplied to program, give me the size of matrix.\n");
        exit(1);
    }
    int size= atoi(argv[1]);
    
    int i,j;
    num**A=(num**) malloc(size*sizeof(num*));
    num*B=malloc(size*sizeof(num));
    for(i=0;i<size;i++)
        A[i] = malloc(size* sizeof(num));
    printf("Give A matrix parameters:\n");
    for(i=0;i<size;i++)
        for(j=0;j<size;j++){
            scanf("%lf",&A[i][j]);
        }
    printf("Give B matrix parameters: \n");
    for(i=0;i<size;i++)
        scanf("%lf",&B[i]);

    num**AB=(num**)malloc(size*sizeof(num*));
    for(i=0;i<size;i++)
        AB[i]=(num*)malloc((size+1)*sizeof(num));
    
    for(i=0;i<size;i++){
        for(j=0;j<size+1;j++){
            if(j!=size)
                AB[i][j]=A[i][j];
            else
                AB[i][j]=B[i];
        }
    }
    num*result= (num*)malloc(size*sizeof(num));

    for(i=0;i<size;i++){
        for(j=0;j<size+1;j++){
            printf("%lf ",AB[i][j]);
        }
        printf("\n");
    }
    if(gaussian_elimination_1(size,AB,result)==false){
        printf("Unsolvable matrixes given\n");
        exit(1);
    }
    printf("Results:\n");
    for(i=0;i<size;i++)
        printf("%lf  ",result[i]);

    return 0;
    
}