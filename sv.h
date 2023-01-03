
typedef struct _variant
{
	char* sv_type; 
	int sv_size; 
	int ref_size; 
	int ref_start; 
	int ref_end;
    char* contig;
	char* reads[101];
	char* node;
	char node_strand;
} variant;
