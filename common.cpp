#include <iostream>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "common.h"

// Track memory usage
long long memUsage = 0;

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

	memUsage = memUsage + size;
	return ret;
}

void init_params(parameters** params)
{
	/* initialize parameters */
	*params = (parameters*) getMem(sizeof(parameters));
	(*params)->ref_graph = NULL;
	(*params)->gaf = NULL;
}

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

/*void copy_string(const string &other)
{

	int len = strlen(other);
	other = new char [len + 1];
	strcpy(m_pString, other.m_pString);
 	return *this;
}*/

/*
char *substring(char *string, int position, int length)
{
   char *p;
   int c;

   p = (char*) malloc(length+1);

   if (p == NULL)
   {
      printf("Unable to allocate memory.\n");
      exit(1);
   }

   for (c = 0; c < length; c++)
   {
      *(p+c) = *(string+position-1);
      string++;
   }

   *(p+c) = '\0';

   return p;
}

char* substr(const char *src, int start_index, int end_index)
{
    int len = end_index - start_index;
    char *dest = (char*) malloc(sizeof(char) * (len + 1));
	
	strncpy(dest, (src + start_index), len);

    return dest;
}



double getMemUsage()
{
	return memUsage / 1048576.0;
}

FILE* safe_fopen( char* path, char* mode)
{
	FILE* file;
	char err[500];

	file = fopen(path, mode);
	if(!file)
	{
		sprintf(err, "[PSVPAN INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error(err);
	}
	return file;
}

gzFile safe_fopen_gz( char* path, char* mode)
{
	//Safe file open. Try to open a file; exit if file does not exist
	gzFile file;
	char err[500];

	file = gzopen( path, mode);
	if(!file)
	{
		sprintf(err, "[PSVPAN INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error(err);
	}
	return file;
}
*/
void print_error( char* msg)
{
	/* print error message and exit */
	printf("\n%s\n", msg);
	printf("Invoke parameter -h for help.\n");
	exit(EXIT_COMMON);
}
