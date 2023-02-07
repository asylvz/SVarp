#ifndef __COMMON
#define __COMMON

#include <fstream>
#include <vector>

#define EXIT_SUCCESS 0
#define EXIT_COMMON 1

#define RETURN_SUCCESS 0
#define RETURN_ERROR -1

#define HASHSIZE 1001
#define MINSVSIZE 50

#define MINMAPQ 20
#define TEST_SAMPLE_SIZE 5000000
#define MIN_READ_SUPPORT 5
#define MIN_SV_DISTANCE 500

#define DELETION 'D'
#define INSERTION 'I'

const std::string REMAP_OUTPUT = "remap_output.gaf";
const std::string FASTA_OUTPUT = "merged_cns.fa";

typedef struct _parameters
{
	std::string gaf;
	std::string ref_graph;
	std::string fasta;
	std::string phase_tags;

	_parameters() {
    }
} parameters;

extern std::ofstream logFile; //Defined in svarp.c

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
int decompose_cigars(std::string cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp);

#endif
