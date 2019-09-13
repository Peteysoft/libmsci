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

  //for empirically designed multi-class models:
  nel_ta ntrain;		//number of training samples
  dim_ta D;			//dimensionality of problem
  real_a **x=NULL;		//features data
  cls_ta *cls=NULL;		//class data
  real_a *d=NULL;		//Hausdorff distance btw. all pairs of classes
				//(upper or lower triangular?)

  opt_args.Qtype=0;
  err=agf_parse_command_opts(argc, argv, "Q:GnS:Yy:P", &opt_args);
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
    fprintf(docfs, "  -P   = print out straight coding matrix\n");
    exit(0);
  }

  if (opt_args.Qtype!=6 && opt_args.Qtype!=9 && opt_args.Qtype!=12) {
    err=sscanf(argv[0], "%d", &n);
    if (err!=1) {
      fprintf(stderr, "print_control: error parsing command line first argument\n");
      exit(COMMAND_OPTION_PARSE_ERROR);
    }
  }

  //read in training data for designing multi-class model:
  if (opt_args.Qtype==6 || opt_args.Qtype==12) {
    char *fname;		//filename
    nel_ta n1;			//number of samples in class data
    int rmflag=0;		//remove normalization file?

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
    if (rmflag) remove(opt_args.normfile);
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
        //dendrogram based on distance between classes:
        cluster_tree<real_a, cls_ta> dg;
        char options[3]="\"\"";

        dg.build_all(d, n, D);
        dg.print(stdout, options);
      }
      break;
    case(7):
      coding_matrix=hierarchical_nonhierarchical<int>(n);
      nrow=n-1;
      break;
    case(8):
      int sum;
      int **codet;
      //int **cm;
      //int **product;
      nrow=4*((n-1)/4+1);
      codet=NULL;
      do {
        if (codet!=NULL) delete_matrix(codet);
        codet=ortho_coding_matrix_greedy<int>(n, nrow, opt_args.Yflag);
	//pile brute force upon brute force:
	//(making sure there are classes on both sides of the fence...)
	//print_matrix(stdout, coding_matrix, nrow, n);
	for (int i=0; i<nrow; i++) {
	  sum=0;
          for (int j=0; j<n; j++) sum+=codet[j][i];
	  if (abs(sum)==n) break;
	}
      } while (abs(sum)==n);
      coding_matrix=matrix_transpose(codet, n, nrow);
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
      //should move this mess into the function...
      int bflag;
      int np;
      if (argc>1) {
        nrow=atoi(argv[1]);
	if (argc>2) np=atoi(argv[2]); else np=-1;
      } else {
        nrow=4*((n-1)/4+1);
      }
      codet=NULL;
      //do {
      //  if (codet!=NULL) delete_matrix(codet);
        codet=orthogonal_coding_matrix<int>(n, nrow, np);
        if (opt_args.Pflag) {
          print_matrix(stdout, codet, nrow, n);
          printf("\n");
        }
        //pile brute force upon brute force:
        //(making sure there are classes on both sides of the fence...)
        //print_matrix(stdout, coding_matrix, nrow, n);
	//*** doesn't check for repeats ...
        for (int i=0; i<nrow; i++) {
          for (int j=0; j<i; j++) {
            bflag=1;
            for (int k=0; k<n; k++) {
              if (codet[i][k] != codet[j][k]) {
                bflag=0;
	        break;
	      }
              if (bflag) break;
            }
            if (bflag) break;
          }
          for (int j=0; j<i; j++) {
            bflag=1;
            for (int k=0; k<n; k++) {
              if (codet[i][k] !=  - codet[j][k]) {
                bflag=0;
	        break;
              }
              if (bflag) break;
            }
	    if (bflag) break;
	  }
	/* dumb idea: destroys orthogonality property...
	  if (bflag=0 || check_coding_row(codet[i], n)==0) {
            bflag=0;
	    break;
	  }
        }
      } while (bflag==0);
	*/
	//remove bad rows:
        if (bflag=0 || check_coding_row(codet[i], n)==0) {
          nrow--;
	  //delete [] codet[i];
	  for (int j=i; j<nrow; j++) codet[j]=codet[j+1];
	  i--;
	}
      }
      coding_matrix=codet;
      break;
    case(11):
      nrow=4*((n-1)/4+1);
      coding_matrix=orthogonal_coding_matrix<int>(n, nrow);
      nrow=n;
      break;
    case(12): {
        real_a *dtriangle;
        nel_ta k;
        nrow=atoi(argv[1]);
        dtriangle=new real_a[ntrain*(ntrain-1)/2+1];
        for (nel_ta i=0; i<ntrain; i++) {
          for (nel_ta j=0; j<i; j++) {
            dtriangle[k]=metric2(x[i], x[j], D);
            k++;
          }
        }
        coding_matrix=optimal_coding_matrix<int>(dtriangle, cls, ntrain, nrow);
	delete [] dtriangle;
      }
      break;
    case(16):
      nrow=4*((n-1)/4+1);
      coding_matrix=random_coding_matrix<int>(n, nrow, 1);
      break;
    default:
      print_control_hier(stdout, n);
  }

  if (opt_args.Qtype > 0 && opt_args.Qtype!=6) {
    //if (opt_args.Qtype==5 || (opt_args.Qtype==8 && opt_args.Yflag)) {
    if (opt_args.Pflag) {
      print_matrix(stdout, coding_matrix, nrow, n);
      printf("\n");
    }
    if (opt_args.Qtype==9 && opt_args.Yflag) {
      print_control_nonhier(stdout, coding_matrix+1, nrow-1, n, name, label);
    } else {
      print_control_nonhier(stdout, coding_matrix, nrow, n, name, label);
    }
    //delete [] coding_matrix[0];
    //delete [] coding_matrix;
  }

  if (name!=NULL) delete [] name;
  if (label!=NULL) delete [] label;
  delete [] d;
  ran_end();

  delete [] cls;
  if (x!=NULL) {
    delete [] x[0];
    delete [] x;
  }

  exit(0);

}
