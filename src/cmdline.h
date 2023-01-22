#ifndef __CMDLINE
#define __CMDLINE

#include "common.h"

void print_help();
int parse_command_line(int argc, char** argv, parameters* params);
#endif
