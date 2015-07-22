#include "agf_lib.h"

using namespace libagf;

int main(int argc, char **argv) {
  FILE *docfs=stdout;
  agf_command_opts opt_args;
  int n;
  int err;

  err=agf_parse_command_opts(argc, argv, "Q:", &opt_args);
  if (err!=0) exit(err);

  if (argc<1) {
    fprintf(docfs, "syntax: print_control [-Q type] n\n");
    fprintf(docfs, "where:\n");
    fprintf(docfs, "  type = 0 hierarchical (default)\n");
    fprintf(docfs, "         1 one against all\n");
    fprintf(docfs, "         2 partition adjacent classes\n");
    fprintf(docfs, "  n    = number of classes\n");
    exit(0);
  }

  err=sscanf(argv[0], "%d", &n);
  if (err!=1) {
    fprintf(stderr, "print_control: error parsing command line first argument\n");
    exit(COMMAND_OPTION_PARSE_ERROR);
  }

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
    default:
      print_control_hier(stdout, n);
  }

  exit(0);

}
