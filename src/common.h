#ifndef __COMMON
#define __COMMON

// Track memory usage
//extern long long memUsage;

#define EXIT_SUCCESS 0
#define EXIT_COMMON 1

#define RETURN_SUCCESS 0
#define RETURN_ERROR -1

#define HASHSIZE 1001
#define MINSVSIZE 50

#define TEST_SAMPLE_SIZE 100000

typedef struct _parameters
{
	std::string gaf;
	std::string ref_graph;
	std::string fasta;
	
	_parameters() {
    }
} parameters;

//#define max(x, y) (((x) > (y)) ? (x) : (y))
//#define min(x, y) (((x) < (y)) ? (x) : (y))

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
