#include <stdio.h>

#include "ctraj_defaults.h"

namespace ctraj {

void ctraj_optargs(FILE *fs, const char *optargs, int flag) {
  int qflag=0;

  qflag=flag && (int) 1;

  for (int i=0; optargs[i]!='\0'; i++) {
    switch (optargs[i]) {
      case ('-'):
	//ugly, but the desired character is already in use...
        fprintf(fs, "  --   restrict to Southern hemisphere\n");
        break;
      case ('+'):
	//ugly, but the desired character is already used...
        fprintf(fs, "  -+   restrict to Northern hemisphere\n");
        break;
      case ('0'):
        fprintf(fs, "  -0   initial time index [0]\n");
        break;
      case ('1'):
        fprintf(fs, "  -1   index of record in first file [0]\n");
        break;
      case ('2'):
        fprintf(fs, "  -2   index of record in second file [0]\n");
        break;
      case ('a'):
        fprintf(fs, "  -a   missing/fill value\n");
        break;
      case ('A'):
        fprintf(fs, "  -A   number of Arnoldi vectors [%d]\n", NARNOLDI);
        break;
      case ('b'):
        fprintf(fs, "  -b   boundary condition [0.]\n");
        fprintf(fs, "         <number>=fixed\n");
        fprintf(fs, "         f=floating\n");
        break;
      case ('B'):
        fprintf(fs, "  -B   data file page size in bytes\n");
        break;
      case ('c'):
        fprintf(fs, "  -c   maximum degrees of arc between pairs of points [%g]\n", MAXARC);
        break;
      case ('C'):
        fprintf(fs, "  -C   count the number of measurements in PC proxy\n");
        break;
      case ('d'):
        fprintf(fs, "  -d   width of date field [%d]\n", TFIELD_WIDTH);
        break;
      case ('D'):
        fprintf(fs, "  -D   factal dimension type [u]\n");
        fprintf(fs, "         u=uncertainty exponent\n");
        fprintf(fs, "         0=box-counting dimension\n");
        fprintf(fs, "         l=box-counting dimension using line segments\n");
        fprintf(fs, "         (types can be combined, e.g. u0l)\n");
        break;
      case ('e'):
        fprintf(fs, "  -e   number of epsilon values [%d]\n", NEPS);
        break;
      case ('E'):
        fprintf(fs, "  -E   generate meridionally symmetric field\n");
        break;
      case ('f'):
        fprintf(fs, "  -f   final date [last available]\n");
        break;
      case ('F'):
        fprintf(fs, "  -F   final value of dependent variable (epsilon/z) [%f]\n", EPS_MAX);
        break;
      case ('g'):
        fprintf(fs, "  -g   use geometric progression\n");
        break;
      case ('h'):
        fprintf(fs, "  -h   non-dimensional 'coarse' time step [%g]\n", TSTEP_COARSE);
        break;
      case ('H'):
        //fprintf(fs, "  -H   for generating histograms (?)\n");
        break;
      case ('i'):
        fprintf(fs, "  -i   initial date [first available]\n");
        break;
      case ('I'):
        fprintf(fs, "  -I   initial value for dependent variable (epsilon/z) [%f]\n", EPS_MIN);
        break;
      case ('k'):
        fprintf(fs, "  -k   number of R-K 'fine' time steps per coarse step [%d]\n", TSTEP_NFINE);
        break;
      case ('K'):
        fprintf(fs, "  -K   fit constant term in PC proxy\n");
        break;
      case ('l'):
        fprintf(fs, "  -l   lead time for the measurement window [default is the same as for the tracer]\n");
        break;
      case ('L'):
        fprintf(fs, "  -L   number of dimensions [2]\n");
        break;
      case ('m'):
        fprintf(fs, "  -m   minimum change in path in km [%g]\n", DSMIN);
        break;
      case ('M'):
        fprintf(fs, "  -M   maximum number of points in advected contour\n", DSMIN);
        break;
      case ('n'):
        if (qflag) {
          fprintf(fs, "  -n   grid-points per side [%d]\n", NGRID_Q);
        } else {
          fprintf(fs, "  -n   grid-points per side [%d]\n", NGRID);
        }
        break;
      case ('N'):
        fprintf(fs, "  -N   number of time-grids [to end of file]\n");
        break;
      case ('o'):
        fprintf(fs, "  -o   wrap off\n");
        break;
      case ('O'):
        fprintf(fs, "  -O   write interval [%d]\n", WRITE_INT);
        break;
      case ('p'):
        fprintf(fs, "  -p   data path [./]\n");
        break;
      case ('P'):
        fprintf(fs, "  -P   calculate Pearson correlation coefficient\n");
        break;
      case ('q'):
        fprintf(fs, "  -q   sidelengths of tracer field in the form:\n");
        fprintf(fs, "       d1xd2xd3...\n");
        break;
      case ('Q'):
        fprintf(fs, "  -Q   query time grids (write to stdout)\n");
        break;
      case ('r'):
        if (qflag) {
          fprintf(fs, "  -r   sidelength/2 [%g]\n", SIDELENGTH_Q);
        } else {
          fprintf(fs, "  -r   sidelength/2 [%g]\n", SIDELENGTH);
        }
        break;
      case ('R'):
        fprintf(fs, "  -R   (return) record size\n");
        break;
      case ('s'):
        fprintf(fs, "  -s   maximum change in path in km [%g]\n", DSMAX);
        break;
      case ('S'):
        fprintf(fs, "  -S   sort measurement data by date\n");
        break;
      case ('t'):
        fprintf(fs, "  -t   number of constant parameters\n");
        break;
      case ('T'):
        fprintf(fs, "  -T   use theta (potential temperature) levels\n");
        break;
      case ('u'):
        fprintf(fs, "  -u   epsilon-uncertain: minimum number [%d]\n", NUNC_MIN);
        break;
      case ('U'):
        fprintf(fs, "  -U   epsilon-uncertain: maximum trials [%d]\n", UNC_MAXN);
        break;
      case ('V'):
        fprintf(fs, "  -V   type of velocity/tracer field [0]\n");
        fprintf(fs, "         0=2-D global, azimuthal equidistant projection\n");
        fprintf(fs, "         1=2-D Cartesian, date-based time grid\n");
        fprintf(fs, "         2=analytical/n-D Cart.\n");
        break;
      case ('v'):
        fprintf(fs, "  -v   number of singular vectors [%d]\n", NEIG);
        break;
      case ('w'):
        fprintf(fs, "  -w   convert vertical velocities\n");
        break;
      case ('W'):
        fprintf(fs, "  -W   output all (interpolated) tracer fields\n");
        break;
      case ('x'):
        if (qflag) {
          fprintf(fs, "  -x   number of x grids [%d]\n", NGRID_Q);
        } else {
          fprintf(fs, "  -x   number of longitude grids [%d]\n", NLON);
        }
        break;
      case ('X'):
        fprintf(fs, "  -X   length of tracer field in x-direction [%g]\n", SIDELENGTH_Q*2);
        break;
      case ('y'):
        if (qflag) {
          fprintf(fs, "  -y   number of y grids [%d]\n", NGRID_Q);
        } else {
          fprintf(fs, "  -y   number of latitude grids [%d]\n", NLAT);
        }
        break;
      case ('Y'):
        fprintf(fs, "  -Y   length of tracer field in y-direction [%g]\n", SIDELENGTH_Q*2);
        break;
      case ('z'):
        fprintf(fs, "  -z   number of contours [%d]\n", NZ);
        break;
      case ('?'):
        fprintf(fs, "  -?   print this help screen\n");
        break;
    }
  }

}

}

