// Track memory usage
extern long long memUsage;

#define EXIT_SUCCESS 0
#define EXIT_COMMON 1

#define RETURN_SUCCESS 1
#define RETURN_ERROR 0

#define HASHSIZE 101

typedef struct _parameters
{
	char* gaf;
	char* ref_graph;
} parameters;

/* String functions */
void set_str(char **target, char *source); /* Even safer than strncpy */
char* substr(const char *src, int start_index, int end_index);

/* memory */
void* getMem(size_t size);
double getMemUsage();
FILE* safe_fopen( char* path, char* mode);		
//gzFile safe_fopen_gz( char* path, char* mode);
void init_params(parameters** params);
void print_error( char* msg);


unsigned hash(char *s);