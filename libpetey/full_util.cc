#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "error_codes.h"

#include "linked.h"
#include "full_util.h"
#include "read_ascii_all.h"

namespace libpetey {

template <class real, class integer>
real ** allocate_matrix(integer m, integer n) {
  real **mat;

  mat=new real *[m];
  mat[0]=new real[m*n];
  for (integer i=1; i<m; i++) {
    mat[i]=mat[0]+n*i;
  }

  return mat;
}

template <class real, class integer>
void zero_matrix(real ** mat, integer m, integer n) {
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) mat[i][j]=0;
  }
}

template <class real, class integer>
real ** zero_matrix(integer m, integer n) {
  real **result;
  result=allocate_matrix<real, integer>(m, n);
  zero_matrix(result, m, n);
  return result;
}

template <class real, class integer>
void identity_matrix(real ** mat, integer m, integer n) {
  integer dmin;
  zero_matrix(mat, m, n);
  if (m < n) dmin=m; else dmin=n;
  for (integer i=0; i<dmin; i++) mat[i][i]=1;
}

template <class real, class integer>
real ** identity_matrix(integer m, integer n) {
  real **result;
  result=allocate_matrix<real, integer>(m, n);
  identity_matrix(result, m, n);
  return result;
}

template <class real>
void delete_matrix(real ** mat) {
  if (mat==NULL) return;
  delete[] mat[0];
  delete[] mat;
}

template <class real, class integer>
void copy_matrix(real **m1, real **m2, integer m, integer n) {
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      m2[i][j]=m1[i][j];
    }
  }
}

template <class real, class integer>
real ** copy_matrix(real **mat, integer m, integer n) {
  real **result;
  result=allocate_matrix<real, integer>(m, n);
  copy_matrix(mat, result, m, n);
  return result;
}

//should be more efficient:
template <class real, class integer>
void matrix_mult_t(real **plier, real **cand, real **result, integer m, integer p, integer n) {

  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      result[i][j]=0;
      for (integer k=0; k<p; k++) {
        result[i][j]+=plier[i][k]*cand[j][k];
      }
    }
  }

}

template <class real, class integer>
void matrix_mult(real **plier, real **cand, real **result, integer m, integer p, integer n) {

  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      result[i][j]=0;
      for (integer k=0; k<p; k++) {
        result[i][j]+=plier[i][k]*cand[k][j];
      }
    }
  }

}

template <class real, class integer>
real ** matrix_mult(real **plier, real **cand, integer m, integer p, integer n) {
  real **result;
  result=allocate_matrix<real, integer>(m, n);
  matrix_mult(plier, cand, result, m, p, n);
  return result;
}

template <class real, class integer>
void vector_mult(real **plier, real *cand, real *result, integer m, integer n) {
  for (integer i=0; i<m; i++) {
    result[i]=0;
    for (integer j=0; j<n; j++) result[i]+=plier[i][j]*cand[j];
  }
}

template <class real, class integer>
real * vector_mult(real **plier, real *cand, integer m, integer n) {
  real *result;
  result=new real[n];
  vector_mult(plier, cand, result, m, n);
  return result;
}

template <class real, class integer>
void left_vec_mult(real *plier, real **cand, real *result, integer m, integer n) {
  for (integer i=0; i<n; i++) result[i]=0;
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      result[j]+=plier[i]*cand[i][j];
    }
  }
}

template <class real, class integer>
real * left_vec_mult(real *plier, real **cand, integer m, integer n) {
  real *result;
  result=new real[n];
  left_vec_mult<real, integer>(plier, cand, result, m, n);
  return result;
}


template <class real, class integer>
void matrix_add(real **mat1, real ** mat2, integer m, integer n) {
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      mat1[i][j]=mat1[i][j]+mat2[i][j];
    }
  }
}

//for square matrices (inplace):
template <class real, class integer>
void matrix_transpose(real **mat, integer m) {
  real swap;
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<i; j++) {
      swap=mat[j][i];
      mat[j][i]=mat[i][j];
      mat[i][j]=swap;
    }
  }

}

template <class real, class integer>
real ** matrix_transpose(real **mat, integer m, integer n) {

  real **matnew;

  //fuuuuuuuuuuuu
  matnew=allocate_matrix<real, integer>(n, m);

  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      matnew[j][i]=mat[i][j];
    }
  }

  return matnew;
}

template <class real, class integer>
real ** scan_matrix(FILE *fptr, integer &m, integer &n, int headflag) {
  real **mat;
  int line_pos=0;
  int lc;
  int ncon;
  int err=0;
  char format1[4];
  char format2[6];

  if (sizeof(real) == 8) {
    strcpy(format1, "%lg");
    strcpy(format2, "%lg%n");
  } else {
    strcpy(format1, "%g");
    strcpy(format2, "%g%n");
  }

  //does the file include a header?
  if (headflag==0) {
    int64_t m1, n1;
    char *line;

    //"more sophisticated" method is slow on cygwin because memory allocation is slow
    //use a dirt simple read method instead:
    //(want to speed up statistical classification routines that use ASCII files...)
    ncon=fscanf(fptr, "%ld %ld", &m1, &n1);
    if (ncon!=2) {
      fprintf(stderr, "scan_matrix: error reading header\n");
      return NULL;
    }
    m=m1;
    n=n1;
    mat=allocate_matrix<real, integer>(m, n);
    for (integer i=0; i<m; i++) {
      for (integer j=0; j<n; j++) {
        ncon=fscanf(fptr, format1, mat[i]+j);
	if (ncon!=1) {
          fprintf(stderr, "scan_matrix: file ended before finished scanning; row=%d\n", (int32_t) i);
          m=i-1;				//truncate matrix to where failure occured
          err=FILE_READ_ERROR;
          break;
	}
      }
      if (err!=0) break;
    }
    return mat;

    //unreachable code:
    if ((line=fget_line(fptr))==NULL) {
      fprintf(stderr, "scan_matrix: no data found\n");
      return NULL;
    }
    ncon=sscanf(line, "%ld %ld\n", &m1, &n1);
    m=m1;
    n=n1;
    if (ncon!=2) {
      fprintf(stderr, "scan_matrix: error reading header\n");
      return NULL;
    }
    mat=allocate_matrix<real, integer>(m, n);
    delete [] line;

    for (integer i=0; i<m; i++) {
      if ((line=fget_line(fptr))==NULL) {
        fprintf(stderr, "scan_matrix: file ended before finished scanning; row=%d\n", (int32_t) i);
        m=i-1;				//truncate matrix to where failure occured
        err=FILE_READ_ERROR;
        break;
      }
      if (err!=0) break;
      line_pos=0;
      for (integer j=0; j<n; j++) {
        //ncon=sscanf(line+line_pos, "%f%n", mat[i]+j, &lc);
        ncon=sscanf(line+line_pos, format2, mat[i]+j, &lc);
        if (ncon!=1) {
          fprintf(stderr, "scan_matrix: error converting data or not enough on line\n");
          fprintf(stderr, "             row=%d; column=%d\n", (int32_t) i, (int32_t) j);
          m=i-1;			//truncate matrix to where failure occured
          err=FILE_READ_ERROR;
          break;
        }
        line_pos+=lc;
      }
      delete [] line;
      if (err!=0) break;
    }
  } else {
    char **line;
    real val;
    long m1;
    int n1;
    int *loc;

    line=read_ascii_all(fptr, &m1, 1);
    m=m1;
    if (m==0) {
      fprintf(stderr, "scan_matrix: no data found\n");
      return NULL;
    }
    loc=split_string(line[0], n1);
    n=n1;
    if (n==0) {
      fprintf(stderr, "scan_matrix: first line missing data elements\n");
      for (integer i=0; i<m; i++) delete [] line[i];
      delete [] line;
      return NULL;
    }
    mat=allocate_matrix<real, integer>(m, n);

    //for (integer i=0; i<n; i++) sscanf(line[0]+loc[i], "%g", mat[0]+i);
    for (integer i=0; i<n; i++) sscanf(line[0]+loc[i], format1, mat[0]+i);
    delete [] loc;
    delete [] line[0];
    
    for (integer i=1; i<m; i++) {
      line_pos=0;
      for (integer j=0; j<n; j++) {
        //ncon=sscanf(line[i]+line_pos, "%f%n", mat[i]+j, &lc);
        ncon=sscanf(line[i]+line_pos, format2, mat[i]+j, &lc);
        if (ncon!=1) {
          fprintf(stderr, "scan_matrix: error converting data or not enough on line\n");
          fprintf(stderr, "             row=%d; column=%d\n", (int32_t) i, (int32_t) j);
          m=i-1;			//truncate matrix to where failure occured
          err=FILE_READ_ERROR;
          break;
        }
        line_pos+=lc;
      }
      delete [] line[i];
      if (err!=0) break;
    }
    delete [] line;
  }
    
  return mat;
}

//can't get the bloody thing to link otherwise:
template <class real, class integer>
real ** scan_matrix(FILE *fptr, integer &m, integer &n) {
  return scan_matrix<real, integer>(fptr, m, n, 0);
}

template <class integer>
void print_matrix(FILE *fptr, int **mat, integer m, integer n) {
  fprintf(fptr, "%d %d\n", (int32_t) m, (int32_t) n);
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      fprintf(fptr, "%d ", mat[i][j]);
    }
    fprintf(fptr, "\n");
  }
}

template <class integer>
void print_matrix(FILE *fptr, float **mat, integer m, integer n) {
  fprintf(fptr, "%d %d\n", (int32_t) m, (int32_t) n);
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      fprintf(fptr, "%g ", mat[i][j]);
    }
    fprintf(fptr, "\n");
  }
}

template <class integer>
void print_matrix(FILE *fptr, double **mat, integer m, integer n) {
  fprintf(fptr, "%d %d\n", (int32_t) m, (int32_t) n);
  for (integer i=0; i<m; i++) {
    for (integer j=0; j<n; j++) {
      fprintf(fptr, "%20.12lg ", mat[i][j]);
    }
    fprintf(fptr, "\n");
  }
}

//getting too long to be inlined...
template <class real, class integer>
real ** read_matrix(FILE *fptr, integer &m, integer &n) {
  real **mat;
  long check;
  long ntot=0;
  n=0;
  fread(&n, 1, sizeof(n), fptr);

  check=ftell(fptr);
  //not a seekable stream--use a linked list to get the data:
  if (check < 0) {
    linked_list<real> data;
    real *data2;
    real val=0;

    while (fread(&val, sizeof(val), 1, fptr)==1) data.add(val);
    data2=data.make_array(ntot);
    //assert(ntot % nvar==0);
    m=ntot/n;
    mat=new real *[m];
    mat[0]=data2;
    for (integer i=0; i<m; i++) mat[i]=mat[0]+i*n;
    if (ntot % n != 0) {
      fprintf(stderr, "read_matrix: header (%d) and data elements (%d) inconsistent\n", (int32_t) n, (int32_t) ntot);
      m=-1; n=-1;
    }
  } else {
    fseek(fptr, 0, SEEK_END);
    ntot=(ftell(fptr)-sizeof(n));
    //check for consistency::
    if (ntot % (n*sizeof(real)) != 0) {
      fprintf(stderr, "read_matrix: header (%d) and data length (%d) inconsistent\n", (int32_t) n, (int32_t) ntot);
      m=-1; n=-1;
      return NULL;
    }
    m=ntot/sizeof(real)/n;
    fseek(fptr, sizeof(n), SEEK_SET);
    //printf("read_matrix: reading a [%dx%d] matrix\n", m, n);
    mat=allocate_matrix<real, integer>(m, n);
    fread(mat[0], sizeof(real), m*n, fptr);
  }

  return mat;
}

template <class real, class integer>
size_t write_matrix(FILE *fptr, real **mat, integer m, integer n) {
  size_t nwrit;
  nwrit+=fwrite(&n, 1, sizeof(n), fptr);
  nwrit+=fwrite(mat[0], sizeof(real), m*n, fptr);
  return nwrit;
}

template <class real, class integer>
real matrix_norm(real **mat, integer m, integer n) {
  real result=0;
  for (integer i=0; i<m*n; i++) result+=mat[0][i]*mat[0][i];
  return sqrt(result);
}

template float ** allocate_matrix<float, int16_t>(int16_t, int16_t);
template float ** allocate_matrix<float, int32_t>(int32_t, int32_t);
template float ** allocate_matrix<float, int64_t>(int64_t, int64_t);

template double ** allocate_matrix<double, int16_t>(int16_t, int16_t);
template double ** allocate_matrix<double, int32_t>(int32_t, int32_t);
template double ** allocate_matrix<double, int64_t>(int64_t, int64_t);

template int32_t ** allocate_matrix<int32_t, int32_t>(int32_t, int32_t);
template int64_t ** allocate_matrix<int64_t, int32_t>(int32_t, int32_t);

template void zero_matrix<float, int16_t>(float **, int16_t, int16_t);
template void zero_matrix<float, int32_t>(float **, int32_t, int32_t);
template void zero_matrix<float, int64_t>(float **, int64_t, int64_t);

template void zero_matrix<double, int16_t>(double **, int16_t, int16_t);
template void zero_matrix<double, int32_t>(double **, int32_t, int32_t);
template void zero_matrix<double, int64_t>(double **, int64_t, int64_t);

template float ** zero_matrix<float, int16_t>(int16_t, int16_t);
template float ** zero_matrix<float, int32_t>(int32_t, int32_t);
template float ** zero_matrix<float, int64_t>(int64_t, int64_t);

template double ** zero_matrix<double, int16_t>(int16_t, int16_t);
template double ** zero_matrix<double, int32_t>(int32_t, int32_t);
template double ** zero_matrix<double, int64_t>(int64_t, int64_t);

template void identity_matrix<float, int16_t>(float **, int16_t, int16_t);
template void identity_matrix<float, int32_t>(float **, int32_t, int32_t);
template void identity_matrix<float, int64_t>(float **, int64_t, int64_t);

template void identity_matrix<double, int16_t>(double **, int16_t, int16_t);
template void identity_matrix<double, int32_t>(double **, int32_t, int32_t);
template void identity_matrix<double, int64_t>(double **, int64_t, int64_t);

template float ** identity_matrix<float, int16_t>(int16_t, int16_t);
template float ** identity_matrix<float, int32_t>(int32_t, int32_t);
template float ** identity_matrix<float, int64_t>(int64_t, int64_t);

template double ** identity_matrix<double, int16_t>(int16_t, int16_t);
template double ** identity_matrix<double, int32_t>(int32_t, int32_t);
template double ** identity_matrix<double, int64_t>(int64_t, int64_t);

template void delete_matrix<float>(float **);
template void delete_matrix<double>(double **);

template void delete_matrix<int32_t>(int32_t **);
template void delete_matrix<int64_t>(int64_t **);

template int ** copy_matrix<int, int16_t>(int **, int16_t, int16_t);
template int ** copy_matrix<int, int32_t>(int **, int32_t, int32_t);
template int ** copy_matrix<int, int64_t>(int **, int64_t, int64_t);

template float ** copy_matrix<float, int16_t>(float **, int16_t, int16_t);
template float ** copy_matrix<float, int32_t>(float **, int32_t, int32_t);
template float ** copy_matrix<float, int64_t>(float **, int64_t, int64_t);

template double ** copy_matrix<double, int16_t>(double **, int16_t, int16_t);
template double ** copy_matrix<double, int32_t>(double **, int32_t, int32_t);
template double ** copy_matrix<double, int64_t>(double **, int64_t, int64_t);

template void matrix_mult<float, int16_t>(float **, float **, float **, int16_t, int16_t, int16_t);
template void matrix_mult<float, int32_t>(float **, float **, float **, int32_t, int32_t, int32_t);
template void matrix_mult<float, int64_t>(float **, float **, float **, int64_t, int64_t, int64_t);

template void matrix_mult<double, int16_t>(double **, double **, double **, int16_t, int16_t, int16_t);
template void matrix_mult<double, int32_t>(double **, double **, double **, int32_t, int32_t, int32_t);
template void matrix_mult<double, int64_t>(double **, double **, double **, int64_t, int64_t, int64_t);

template int32_t ** matrix_mult<int32_t, int32_t>(int32_t **, int32_t **, int32_t, int32_t, int32_t);

template float ** matrix_mult<float, int16_t>(float **, float **, int16_t, int16_t, int16_t);
template float ** matrix_mult<float, int32_t>(float **, float **, int32_t, int32_t, int32_t);
template float ** matrix_mult<float, int64_t>(float **, float **, int64_t, int64_t, int64_t);

template double ** matrix_mult<double, int16_t>(double **, double **, int16_t, int16_t, int16_t);
template double ** matrix_mult<double, int32_t>(double **, double **, int32_t, int32_t, int32_t);
template double ** matrix_mult<double, int64_t>(double **, double **, int64_t, int64_t, int64_t);

template void matrix_mult_t<float, int16_t>(float **, float **, float **, int16_t, int16_t, int16_t);
template void matrix_mult_t<float, int32_t>(float **, float **, float **, int32_t, int32_t, int32_t);
template void matrix_mult_t<float, int64_t>(float **, float **, float **, int64_t, int64_t, int64_t);

template void matrix_mult_t<double, int16_t>(double **, double **, double **, 
		int16_t, int16_t, int16_t);
template void matrix_mult_t<double, int32_t>(double **, double **, double **, 
		int32_t, int32_t, int32_t);
template void matrix_mult_t<double, int64_t>(double **, double **, double **, 
		int64_t, int64_t, int64_t);

template void vector_mult<float, int16_t>(float **, float *, float *, int16_t, int16_t);
template void vector_mult<float, int32_t>(float **, float *, float *, int32_t, int32_t);
template void vector_mult<float, int64_t>(float **, float *, float *, int64_t, int64_t);

template void vector_mult<double, int16_t>(double **, double *, double *, int16_t, int16_t);
template void vector_mult<double, int32_t>(double **, double *, double *, int32_t, int32_t);
template void vector_mult<double, int64_t>(double **, double *, double *, int64_t, int64_t);

template float * vector_mult<float, int16_t>(float **, float *, int16_t, int16_t);
template float * vector_mult<float, int32_t>(float **, float *, int32_t, int32_t);
template float * vector_mult<float, int64_t>(float **, float *, int64_t, int64_t);

template double * vector_mult<double, int16_t>(double **, double *, int16_t, int16_t);
template double * vector_mult<double, int32_t>(double **, double *, int32_t, int32_t);
template double * vector_mult<double, int64_t>(double **, double *, int64_t, int64_t);

template void left_vec_mult<float, int16_t>(float *, float **, float *, int16_t, int16_t);
template void left_vec_mult<float, int32_t>(float *, float **, float *, int32_t, int32_t);
template void left_vec_mult<float, int64_t>(float *, float **, float *, int64_t, int64_t);

template void left_vec_mult<double, int16_t>(double *, double **, double *, int16_t, int16_t);
template void left_vec_mult<double, int32_t>(double *, double **, double *, int32_t, int32_t);
template void left_vec_mult<double, int64_t>(double *, double **, double *, int64_t, int64_t);

template float * left_vec_mult<float, int16_t>(float *, float **, int16_t, int16_t);
template float * left_vec_mult<float, int32_t>(float *, float **, int32_t, int32_t);
template float * left_vec_mult<float, int64_t>(float *, float **, int64_t, int64_t);

template double * left_vec_mult<double, int16_t>(double *, double **, int16_t, int16_t);
template double * left_vec_mult<double, int32_t>(double *, double **, int32_t, int32_t);
template double * left_vec_mult<double, int64_t>(double *, double **, int64_t, int64_t);

template void matrix_add<float, int16_t>(float **, float **, int16_t, int16_t);
template void matrix_add<float, int32_t>(float **, float **, int32_t, int32_t);
template void matrix_add<float, int64_t>(float **, float **, int64_t, int64_t);

template void matrix_add<double, int16_t>(double **, double **, int16_t, int16_t);
template void matrix_add<double, int32_t>(double **, double **, int32_t, int32_t);
template void matrix_add<double, int64_t>(double **, double **, int64_t, int64_t);

//fuck, half of this file is just gonna be template instantiations...

template void matrix_transpose<float, int16_t>(float **, int16_t);
template void matrix_transpose<float, int32_t>(float **, int32_t);
template void matrix_transpose<float, int64_t>(float **, int64_t);

template void matrix_transpose<double, int16_t>(double **, int16_t);
template void matrix_transpose<double, int32_t>(double **, int32_t);
template void matrix_transpose<double, int64_t>(double **, int64_t);

template int ** matrix_transpose<int, int32_t>(int **, int32_t, int32_t);

template float ** matrix_transpose<float, int16_t>(float **, int16_t, int16_t);
template float ** matrix_transpose<float, int32_t>(float **, int32_t, int32_t);
template float ** matrix_transpose<float, int64_t>(float **, int64_t, int64_t);

template double ** matrix_transpose<double, int16_t>(double **, int16_t, int16_t);
template double ** matrix_transpose<double, int32_t>(double **, int32_t, int32_t);
template double ** matrix_transpose<double, int64_t>(double **, int64_t, int64_t);

template float ** scan_matrix<float, int16_t>(FILE *, int16_t &, int16_t &, int);
template float ** scan_matrix<float, int32_t>(FILE *, int32_t &, int32_t &, int);
template float ** scan_matrix<float, int64_t>(FILE *, int64_t &, int64_t &, int);

template double ** scan_matrix<double, int16_t>(FILE *, int16_t &, int16_t &, int);
template double ** scan_matrix<double, int32_t>(FILE *, int32_t &, int32_t &, int);
template double ** scan_matrix<double, int64_t>(FILE *, int64_t &, int64_t &, int);

template float ** read_matrix<float, int16_t>(FILE *, int16_t &, int16_t &);
template float ** read_matrix<float, int32_t>(FILE *, int32_t &, int32_t &);
template float ** read_matrix<float, int64_t>(FILE *, int64_t &, int64_t &);

template double ** read_matrix<double, int16_t>(FILE *, int16_t &, int16_t &);
template double ** read_matrix<double, int32_t>(FILE *, int32_t &, int32_t &);
template double ** read_matrix<double, int64_t>(FILE *, int64_t &, int64_t &);

template void print_matrix<int16_t>(FILE *, int **, int16_t, int16_t);
template void print_matrix<int32_t>(FILE *, int **, int32_t, int32_t);
template void print_matrix<int64_t>(FILE *, int **, int64_t, int64_t);

template void print_matrix<int16_t>(FILE *, float **, int16_t, int16_t);
template void print_matrix<int32_t>(FILE *, float **, int32_t, int32_t);
template void print_matrix<int64_t>(FILE *, float **, int64_t, int64_t);

template void print_matrix<int16_t>(FILE *, double **, int16_t, int16_t);
template void print_matrix<int32_t>(FILE *, double **, int32_t, int32_t);
template void print_matrix<int64_t>(FILE *, double **, int64_t, int64_t);

template size_t write_matrix<float, int16_t>(FILE *, float **, int16_t, int16_t);
template size_t write_matrix<float, int32_t>(FILE *, float **, int32_t, int32_t);
template size_t write_matrix<float, int64_t>(FILE *, float **, int64_t, int64_t);

template size_t write_matrix<double, int16_t>(FILE *, double **, int16_t, int16_t);
template size_t write_matrix<double, int32_t>(FILE *, double **, int32_t, int32_t);
template size_t write_matrix<double, int64_t>(FILE *, double **, int64_t, int64_t);

template float matrix_norm<float, int16_t>(float **, int16_t, int16_t);
template float matrix_norm<float, int32_t>(float **, int32_t, int32_t);
template float matrix_norm<float, int64_t>(float **, int64_t, int64_t);

template double matrix_norm<double, int16_t>(double **, int16_t, int16_t);
template double matrix_norm<double, int32_t>(double **, int32_t, int32_t);
template double matrix_norm<double, int64_t>(double **, int64_t, int64_t);

} //end namespace libpetey

