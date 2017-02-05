#include <string.h>

#include "randomize.h"
#include "agf_lib.h"

using namespace libpetey;
using namespace libagf;

int main(int argc, char **argv) {
  FILE *docfs=stdout;
  int **coding_matrix=NULL;
  agf_command_opts opt_args;
  int n;
  int nrow;
  int strictflag;
  int err;

  opt_args.Qtype=0;
  err=agf_parse_command_opts(argc, argv, "Q:GnS:Y", &opt_args);
  if (err!=0) exit(err);

  if (argc<1) {
    fprintf(docfs, "syntax: print_control [-Q type] [-G] [-n] [-S nSV] {n | train} [nrow]\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  type   = 0 hierarchical (default)\n");
    fprintf(docfs, "           1 one against all\n");
    fprintf(docfs, "           2 partition adjacent classes\n");
    fprintf(docfs, "           3 random coding matrix (n<65)\n");
    fprintf(docfs, "           4 \"exhaustive\" coding matrix (n<34)\n");
    fprintf(docfs, "           5 orthogonal coding matrix (n%%4==0)\n");
    fprintf(docfs, "           6 one against one\n");
    fprintf(docfs, "           7 design hierarchical scheme based on training data\n");
    fprintf(docfs, "           8 orthogonal coding matrix (completely brute force: n<25; n%%4==0)\n");
    fprintf(docfs, "  n      = number of classes\n");
    fprintf(docfs, "  train  = base name of training data files (if applicable):\n");
    fprintf(docfs, "             .vec for vector data\n");
    fprintf(docfs, "             .cls for class data\n");
    fprintf(docfs, "  nrow = number of rows\n");
    fprintf(docfs, "  -G   = \"strict\" flag: all classes are included in each partition\n");
    fprintf(docfs, "  -Y   = top row of orthognal coding matrix is implicitly all ones\n");
    exit(0);
  }

  if (opt_args.Qtype!=7) {
    err=sscanf(argv[0], "%d", &n);
    if (err!=1) {
      fprintf(stderr, "print_control: error parsing command line first argument\n");
      exit(COMMAND_OPTION_PARSE_ERROR);
    }
  }

  ran_init();

  switch(opt_args.Qtype) {
    case(0):
      print_control_hier(stdout, n);
      break;
    case(1):
      coding_matrix=one_against_all(n);
      nrow=n;
      break;
    case(2):
      coding_matrix=partition_adjacent(n);
      nrow=n-1;
      break;
    case(3):
      nrow=atoi(argv[1]);
      coding_matrix=random_coding_matrix(n, nrow, opt_args.Gflag);
      break;
    case(4):
      coding_matrix=exhaustive_coding_matrix(n);
      nrow=pow(2, n-1)-1;
      break;
    case(5):
      coding_matrix=ortho_coding_matrix_nqbf(n, opt_args.Gflag, opt_args.Yflag);
      nrow=n;
      break;
    case(6):
      coding_matrix=one_against_one(n);
      nrow=(n-1)*n/2;
      break;
    case(7): {
        char *fname;		//filename
        dim_ta D;		//dimensionality of problem
        nel_ta ntrain;		//number of training samples
        real_a **x;		//features data
        cls_ta *cls;		//class data
        real_a *d;		//Hausdorff distance between all classes
	nel_ta n1;		//number of samples in class data
	int rmflag=0;		//remove normalization file?
        //dendrogram based on distance between classes:
        cluster_tree<real_a, cls_ta> dg;
        char options[3]="\"\"";

        if ((opt_args.svd>0 || opt_args.normflag) && opt_args.normfile == NULL) {
          opt_args.normfile=new char[L_tmpnam];
	  tmpnam(opt_args.normfile);
	  rmflag=1;
        }
	x=agf_get_features<real_a>(argv[0], &opt_args, D, ntrain, 1);
	if (x==NULL || ntrain<=0) {
          fprintf(stderr, "Error reading input file, %s.vec\n", argv[0]);
	  exit(FILE_READ_ERROR);
	}
	fname=new char[strlen(argv[0])+5];
	sprintf(fname, "%s.cls", argv[0]);
        cls=read_clsfile<cls_ta>(fname, n1);
	if (cls==NULL || n1<=0) {
          fprintf(stderr, "Error reading input file, %s\n", fname);
	  exit(FILE_READ_ERROR);
	}
	if (cls==NULL || n1<=0) {
          fprintf(stderr, "Sample count mismatch: %d in %s.vec; %d in %s\n", n, argv[0], n1, fname);
	  exit(SAMPLE_COUNT_MISMATCH);
	}
        d=class_triangle(x, cls, ntrain, D);
        n=0;
        for (nel_ta i=0; i<ntrain; i++) if (cls[i]>=n) n=cls[i]+1;
        dg.build_all(d, n, D);
        dg.print(stdout, options);
	delete [] d;
	delete [] cls;
	delete [] x[0];
	delete [] x;
	if (rmflag) remove(opt_args.normfile);
      }
      break;
    case(8):
      coding_matrix=ortho_coding_matrix_brute_force(n, opt_args.Yflag);
      nrow=n;
      break;
    default:
      print_control_hier(stdout, n);
  }

  if (opt_args.Qtype > 0 && opt_args.Qtype!=7) {
    if (opt_args.Qtype==5 || (opt_args.Qtype==8 && opt_args.Yflag)) {
      print_control_nonhier(stdout, coding_matrix+1, nrow-1, n);
    } else {
      print_control_nonhier(stdout, coding_matrix, nrow, n);
    }
    delete [] coding_matrix[0];
    delete [] coding_matrix;
  }

  ran_end();

  exit(0);

}
