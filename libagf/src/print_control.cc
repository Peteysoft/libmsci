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
  err=agf_parse_command_opts(argc, argv, "Q:G", &opt_args);
  if (err!=0) exit(err);

  if (argc<1) {
    fprintf(docfs, "syntax: print_control [-Q type] [-G] n [nrow]\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  type = 0 hierarchical (default)\n");
    fprintf(docfs, "         1 one against all\n");
    fprintf(docfs, "         2 partition adjacent classes\n");
    fprintf(docfs, "         3 random coding matrix (n<65)\n");
    fprintf(docfs, "         4 \"exhaustive\" coding matrix (n<34)\n");
    fprintf(docfs, "         5 orthogonal coding matrix (n<25; n%%4==0)\n");
    fprintf(docfs, "         6 one against one\n");
    fprintf(docfs, "  n    = number of classes\n");
    fprintf(docfs, "  nrow = number of rows\n");
    fprintf(docfs, "  -G   = \"strict\" flag: all classes are included in each partition\n");
    exit(0);
  }

  err=sscanf(argv[0], "%d", &n);
  if (err!=1) {
    fprintf(stderr, "print_control: error parsing command line first argument\n");
    exit(COMMAND_OPTION_PARSE_ERROR);
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
      coding_matrix=ortho_coding_matrix_nqbf(n, opt_args.Gflag);
      nrow=n;
      break;
    case(6):
      coding_matrix=one_against_one(n);
      nrow=(n-1)*n/2;
      break;
    case(7):
      coding_matrix=ortho_coding_matrix_brute_force(n);
      nrow=n;
      break;
    default:
      print_control_hier(stdout, n);
  }

  if (opt_args.Qtype > 0 && opt_args.Qtype < 8) {
    print_control_nonhier(stdout, coding_matrix, nrow, n);
    delete [] coding_matrix[0];
    delete [] coding_matrix;
  }

  ran_end();

  exit(0);

}
