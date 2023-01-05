#ifndef __COMMON
#define __COMMON

// Track memory usage
extern long long memUsage;

#define EXIT_SUCCESS 0
#define EXIT_COMMON 1

#define RETURN_SUCCESS 1
#define RETURN_ERROR 0

#define HASHSIZE 1001

typedef struct _parameters
{
	char* gaf;
	char* ref_graph;
} parameters;

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

/* String functions
char* substr(const char *src, int start_index, int end_index);


void* getMem(size_t size);
double getMemUsage();
FILE* safe_fopen( char* path, char* mode);		
//gzFile safe_fopen_gz( char* path, char* mode);
void print_error( char* msg);

*/
unsigned hash(char *s);
void init_params(parameters** params);
void set_str(char **target, char *source);


#endif
