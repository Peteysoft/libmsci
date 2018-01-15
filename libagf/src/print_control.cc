#include <string.h>

#include "full_util.h"
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
  char **name=NULL;		//list of binary classifiers
  cls_ta *label=NULL;		//list of class labels

  opt_args.Qtype=0;
  err=agf_parse_command_opts(argc, argv, "Q:GnS:Yy:", &opt_args);
  if (err!=0) exit(err);

  if (argc<1) {
    fprintf(docfs, "syntax: print_control [-Q type] [-G] [-n] [-S nSV] {n | file} [nrow]\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  type   = 0 hierarchical (default)\n");
    fprintf(docfs, "           1 one against all\n");
    fprintf(docfs, "           2 partition adjacent classes\n");
    fprintf(docfs, "           3 random coding matrix (n<65)\n");
    fprintf(docfs, "           4 \"exhaustive\" coding matrix (n<34)\n");
    fprintf(docfs, "           5 one against one\n");
    fprintf(docfs, "           6 design hierarchical scheme based on training data\n");
    fprintf(docfs, "           7 hierarchical non-hierarchical\n");
    fprintf(docfs, "           8 orthogonal coding matrix (completely brute force: n<25; n%%4==0)\n");
    fprintf(docfs, "           9 convert hierarchical to non-hierarchical\n");
    fprintf(docfs, "  n      = number of classes\n");
    fprintf(docfs, "  file   = base name of data files (if applicable):\n");
    fprintf(docfs, "             .vec for vector data\n");
    fprintf(docfs, "             .cls for class data\n");
    fprintf(docfs, "  nrow = number of rows\n");
    fprintf(docfs, "  -G   = \"strict\" flag: all classes are included in each partition\n");
    fprintf(docfs, "  -Y   = top row of orthognal coding matrix is implicitly all ones\n");
    exit(0);
  }

  if (opt_args.Qtype!=6 && opt_args.Qtype!=9) {
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
      coding_matrix=one_against_all<int>(n);
      nrow=n;
      break;
    case(2):
      coding_matrix=partition_adjacent<int>(n);
      nrow=n-1;
      break;
    case(3):
      nrow=atoi(argv[1]);
      coding_matrix=random_coding_matrix<int>(n, nrow, opt_args.Gflag);
      break;
    case(4):
      coding_matrix=exhaustive_coding_matrix<int>(n);
      nrow=pow(2, n-1)-1;
      break;
    case(5):
      coding_matrix=one_against_one<int>(n);
      nrow=(n-1)*n/2;
      break;
    case(6): {
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
	x=agf_get_features<real_a>(argv[0], &opt_args, ntrain, D, 0);
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
    case(7):
      coding_matrix=hierarchical_nonhierarchical<int>(n);
      nrow=n-1;
      break;
    case(8):
      int sum;
      //int **cm;
      //int **product;
      nrow=4*((n-1)/4+1);
      coding_matrix=NULL;
      do {
        if (coding_matrix!=NULL) delete_matrix(coding_matrix);
        coding_matrix=ortho_coding_matrix_greedy<int>(nrow, opt_args.Yflag);
	//pile brute force upon brute force:
	//(making sure there are classes on both sides of the fence...)
	//print_matrix(stdout, coding_matrix, nrow, n);
	for (int i=0; i<nrow; i++) {
	  sum=0;
          for (int j=0; j<n; j++) sum+=coding_matrix[i][j];
	  if (abs(sum)==n) break;
	}
      } while (abs(sum)==n);
      break;
    case(9):
      multiclass_hier<real_a, cls_ta> *dum;
      FILE *fs;
      cls_ta ncls;
      fs=fopen(argv[0], "r");
      dum=new multiclass_hier<real_a, cls_ta>(fs, 0, opt_args.path, "");
      fclose(fs);
      //second method:
      dum->get_code(coding_matrix, name, nrow, n);
      label=new cls_ta[n];
      dum->class_list(label);
      break;
    case(10):
      coding_matrix=ortho_coding_matrix_brute_force<int>(n);
      nrow=n;
      break;
    default:
      print_control_hier(stdout, n);
  }

  if (opt_args.Qtype > 0 && opt_args.Qtype!=6) {
    //if (opt_args.Qtype==5 || (opt_args.Qtype==8 && opt_args.Yflag)) {
    //print_matrix(stdout, coding_matrix, nrow, n);
    if (opt_args.Qtype==9 && opt_args.Yflag) {
      print_control_nonhier(stdout, coding_matrix+1, nrow-1, n, name, label);
    } else {
      print_control_nonhier(stdout, coding_matrix, nrow, n, name, label);
    }
    delete [] coding_matrix[0];
    delete [] coding_matrix;
  }

  if (name!=NULL) delete [] name;
  if (label!=NULL) delete [] label;
  ran_end();

  exit(0);

}
