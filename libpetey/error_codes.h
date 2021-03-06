//
// Error codes (should be updated to Linux standard error codes)
//

#ifndef AGF_ERR_H_INCLUDED
#define AGF_ERR_H_INCLUDED 1

#include <stdio.h>

//error codes:
//last digit: warning or fatal
//second last digit: type of error
//third last digit: class of error
#define INSUFFICIENT_COMMAND_ARGS 1
#define UNABLE_TO_OPEN_FILE_FOR_READING 101
#define FILE_READ_ERROR 111
#define FILE_READ_WARNING 256
#define ALLOCATION_FAILURE 303
#define UNABLE_TO_OPEN_FILE_FOR_WRITING 201
#define FILE_WRITE_ERROR 211
#define FILE_WRITE_WARNING 512
#define DIMENSION_MISMATCH 401
#define SAMPLE_COUNT_MISMATCH 411
#define NUMERICAL_ERROR 501
#define SYNTAX_ERROR 511
#define MAX_ITER_EXCEEDED 503
#define OUT_OF_DATA 523
#define PARAMETER_OUT_OF_RANGE 768
#define COMMAND_OPTION_PARSE_ERROR 1024
#define FATAL_COMMAND_OPTION_PARSE_ERROR 21
#define INTERNAL_ERROR 901
#define OTHER_ERROR 911
#define OTHER_WARNING 1280

#endif
