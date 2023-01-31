#include <iostream>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "common.h"



/* The codes for the dictionary is taken from The C Programming language 
 * 2nd edition - Brian Kernighan and Dennis Ritchie */
unsigned hash(char *s)
{
    unsigned hashval;
    for (hashval = 0; *s != '\0'; s++)
      hashval = *s + 31 * hashval;
    return hashval % HASHSIZE;
}

void* getMem(size_t size)
{
	void* ret;

	ret = malloc(size);
	if(ret == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(0);
	}

	return ret;
}

/*void init_params(parameters** params)
{
	// initialize parameters
	*params = (parameters*) getMem(sizeof(parameters));
	(*params)->ref_graph = NULL;
	(*params)->gaf = NULL;
}*/

void set_str(char** target, char* source)
{
	if(*target != NULL)
	{
		free((*target));
	}

	if(source != NULL)
	{
		(*target) = (char*) getMem(sizeof(char) * (strlen(source) + 1));
		strncpy( (*target), source, (strlen(source) + 1));
	}
	else
	{
		(*target) = NULL;
	}
}

void print_error( char* msg)
{
	/* print error message and exit */
	printf("\n%s\n", msg);
	printf("Invoke parameter -h for help.\n");
	exit(EXIT_COMMON);
}


int decompose_cigars(std::string cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp)
{
	/*get the Cigar*/
	size_t cigar_offset = 0, str_offset = 0, cigar_cnt = 0;
	char* cigar_copy = (char*) cigar.c_str();	
	char *tmp_str = new char[6];
	while(cigar_offset < cigar.length())
	{
		if (isdigit(*(cigar_copy + cigar_offset)) == 0)
		{
			cigarOp.push_back (*(cigar_copy + cigar_offset));
			cigarLen.push_back (atoi(tmp_str));
			
			delete[] tmp_str; 	
			tmp_str = new char[6];
			str_offset = 0;
			cigar_cnt++;			
		}
		else 
		{
			*(tmp_str + str_offset) = *(cigar_copy + cigar_offset);
			str_offset++;
		}
		cigar_offset++;
	}
	delete[] tmp_str;
	
	return cigar_cnt;
}
