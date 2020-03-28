#include <stdio.h>
#include <stdlib.h>

void print_number(int **number, int n) {
  if (n==0) {
    printf("1");
  } else {
    for (int i=0; number[n][i]>=0; i++) {
      printf("(");
      print_number(number, number[n][i]);
      printf(")");
    }
  }
}

int main(int argc, char **argv) {
  int j;
  int n;
  int **number;
  int *prime;
  int nprime;
  int *factor;
  int nfactor;

  if (argc != 2) {
    printf("\n");
    printf("syntax: count n\n");
    printf("\n");
    printf("Counts to n using a special numbering system.\n");
    printf("\n");
    return 0;
  }

  n=atoi(argv[1]);

  number=new int *[n];
  prime=new int[n];

  number[0]=new int[1];
  number[0][0]=-1;
  number[1]=new int[2];
  number[1][0]=0;
  number[1][1]=-1;
  prime[0]=2;
  nprime=1;

  factor=new int[n];

  for (int i=2; i<n; i++) {
    int n=i+1;
    nfactor=0;
    for (j=0; j<nprime && prime[j]<=n; j++) {
      //printf("%d %d %d\n", i+1, j, n);
      if (n % prime[j] == 0) {
        factor[nfactor]=j;
	nfactor++;
	n=n/prime[j];
	j=-1;
      }
      //printf("%d %d %d\n", i+1, j, n);
    }
    if (nfactor==0) {
      prime[nprime]=i+1;
      nfactor=1;
      factor[0]=nprime;
      nprime++;
    }
    number[i]=new int[nfactor+1];
    for (int j=0; j<nfactor; j++) number[i][j]=factor[j];
    number[i][nfactor]=-1;
  }

  /*
  for (int i=0; i<nprime; i++) {
    printf("%4d %4d\n", i+1, prime[i]);
  }
  printf("\n");
  for (int i=0; i<n; i++) {
    printf("%4d:", i+1);
    for (int j=0; number[i][j]>=0; j++) printf(" %4d", number[i][j]+1);
    printf("\n");
  }
  */

  for (int i=0; i<nprime; i++) {
    printf("%4d %5d ", i+1, prime[i]);
    print_number(number, prime[i]-1);
    printf("\n");
  }
  printf("\n");

  for (int i=0; i<n; i++) {
    printf("%5d ", i+1);
    print_number(number, i);
    printf("\n");
  }

}

