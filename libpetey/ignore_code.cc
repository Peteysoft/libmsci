#include <stdio.h>
#include <stdlib.h>

//simple utility to ignore certain exit codes and return 0 instead...

int main (int argc, char **argv) {

  int ncodes=argc-2;
  int code[ncodes];
  int err;
  char *command=argv[argc-1];

  //to prevent the process from making innumerable copies of itself:
  if (argc < 2) {
    fprintf(stderr, "Usage: ignore code [code1 [code2 [code3 ...]]] [\"]command[\"]\n");
    fprintf(stderr, "Purpose: ignores listed exit codes from command and returns\n");
    fprintf(stderr, "0 instead.\n");
    exit(1);
  } 

  for (int i=0; i< ncodes; i++) sscanf(argv[i+1], "%d", code+i);

  err=system(command);
  //printf("%s returned %d exit status\n", command, err);

  for (int i=0; i<ncodes; i++) {
    if (code[i] == err) {
      fprintf(stderr, "Warning: %s returned %d exit status\n", command, err);
      exit(0);
    }
  }

  return err;
} 
  
