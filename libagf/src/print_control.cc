#include "randomize.h"
#include "agf_lib.h"

using namespace libpetey;
using namespace libagf;

int main(int argc, char **argv) {
  FILE *docfs=stdout;
  agf_command_opts opt_args;
  int n;
  int nrow;
  int strictflag;
  int err;

  err=agf_parse_command_opts(argc, argv, "Q:G", &opt_args);
  if (err!=0) exit(err);

  if (argc<1) {
    fprintf(docfs, "syntax: print_control [-Q type] n [nrow]\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  type = 0 hierarchical (default)\n");
    fprintf(docfs, "         1 one against all\n");
    fprintf(docfs, "         2 partition adjacent classes\n");
    fprintf(docfs, "         3 random coding matrix\n");
    fprintf(docfs, "         4 \"exhaustive\" coding matrix\n");
    fprintf(docfs, "  n    = number of classes\n");
    fprintf(docfs, "  nrow = number of rows\n");
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
      print_control_1vsall(stdout, n);
      break;
    case(2):
      print_control_adj(stdout, n);
      break;
    case(3):
      nrow=atoi(argv[1]);
      print_control_random(stdout, n, nrow, opt_args.Gflag);
      break;
    case(4):
      print_control_exhaustive(stdout, n);
      break;
    default:
      print_control_hier(stdout, n);
  }

  ran_end();

  exit(0);

}
