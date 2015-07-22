#ifndef PARSE_COMMAND_OPTS_H_INCLUDED
#define PARSE_COMMAND_OPTS_H_INCLUDED


namespace libpetey {

  //returns number of options found
  //if there is a non-fatal error, returns -(# of found options)

  int parse_command_opts(int argc,	// number of command line args 
		char **argv, 		// arguments passed to command line
		const char *code, 	// code for each option
		const char *format, 	// format code for each option
		void **parm,		// returned parameters
                int *flag,		// found flags
		int opts=0);		// option flags (optional)

}

#endif

