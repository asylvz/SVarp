#ifndef __COMMON
#define __COMMON

#include <fstream>

#define EXIT_SUCCESS 0
#define EXIT_COMMON 1

#define RETURN_SUCCESS 0
#define RETURN_ERROR -1

#define HASHSIZE 1001
#define MINSVSIZE 50

#define TEST_SAMPLE_SIZE 10000000
#define MIN_READ_SUPPORT 5
#define MIN_SV_DISTANCE 500

typedef struct _parameters
{
	std::string gaf;
	std::string ref_graph;
	std::string fasta;
	std::string phase_tags;

	_parameters() {
    }
} parameters;

extern std::ofstream logFile; //Defined in psvpan.c

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
