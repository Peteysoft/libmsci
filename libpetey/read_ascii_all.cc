#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "linked.h"

#include "read_ascii_all.h"

namespace libpetey {

  char * fget_line(FILE *fs, int flag) {
    linked_list<char> line;
    int c1;
    char c2;
    char *result;
    long n;		//not used

    if (flag) {
      //remove line break:
      c1=fgetc(fs);
      c2=(char) c1;
      while (feof(fs)==0 && c1!=EOF && c2 != '\n') {
        line.add(c2);
        c1=fgetc(fs);
        c2=(char) c1;
      }
    } else {
      //keep line break:
      do {
        c1=fgetc(fs);
        if (c1==EOF || feof(fs)) break;
        c2=(char) c1;
        line.add(c2);
      } while (c2!='\n');
    }
    //EOF and only EOF returns NULL pointer:
    if (feof(fs)==0) {
      line.add('\0');
    }
    result=line.make_array(n);

    return result;
  }

  //reads an entire file in and returns it in an array of strings:
  //uses a linked list...
  char ** read_ascii_all(FILE *fs, long *n, int flag) {
    char *line;
    char **all;

    linked_list<char *> data;

    //we keep the newline if flag is 1, delete it if 0...
    while ((line=fget_line(fs, flag)) != NULL) {
      //printf("%s\n", line);
      data.add(line);
    }

    all=data.make_array(*n);

    return all;
  }

  //gee, not much to it.  Hardly seems worth writing at all...
  //(even less too it now that we've writtenn fget_line)

  char ** read_ascii_all(const char *filename, long *n, int maxll) {
    char **result;
    FILE *fs;
    fs=fopen(filename, "r");
    if (fs==NULL) return NULL;
    result=read_ascii_all(fs, n, maxll);
    fclose(fs);
    return result;
  }

  //find the beginnings of all whitespace-separated sub-strings:
  int * split_string(const char *string, int &n) {
    linked_list<int> loc;	//locations
    long n1;
    int i, j;
    int finish=0;
    int *result;

    i=0; 
    do {
      for (; isspace(string[i]); i++) {
        if (string[i]=='\0') {
          finish=1;
          break;
        }
      }
      if (finish) break;
      loc.add(i);
      for (; isspace(string[i])==0; i++) {
        if (string[i]=='\0') {
          finish=1;
          break;
        }
      }
    } while (finish==0);

    result=loc.make_array(n1);
    n=n1;

    return result;
  }

  //find the beginnings of all whitespace-separated sub-strings:
  char ** split_string_destructive(char *string, int &n) {
    char **result;
    linked_list<char *> loc;	//locations
    long n1;
    int i, j;
    int finish=0;

    i=0; 
    while (string[i]!='\0') {
      for (; isspace(string[i]); i++) {
        if (string[i]=='\0') {
          finish=1;
          break;
        }
      }
      if (finish) break;
      loc.add(string+i);
      for (; isspace(string[i])==0; i++) {
        if (string[i]=='\0') {
          finish=1;
          break;
        }
      }
      if (finish) break;
      string[i]='\0';
      i++;
    }

    result=loc.make_array(n1);
    n=n1;
    return result;
  }

} //end namespace libpetey

