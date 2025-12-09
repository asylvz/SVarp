#ifndef __COMMON
#define __COMMON

#include <iosfwd>
#include "logfile.h"
#include <string>
#include <vector>
#include <map>

#define EXIT_SUCCESS 0
#define EXIT_COMMON 1

#define RETURN_SUCCESS 0
#define RETURN_ERROR -1

#define HASHSIZE 1001
#define MINSVSIZE 50

#define MINMAPQ 5
#define MINMAPQREMAP 5
#define TEST_SAMPLE_SIZE 250000000
#define MAX_FETCH_LEN 1000000
#define MIN_READ_START_END_WINDOW 200
#define MAX_CONTIG_DEPTH 100 //>100X coverage for a contig is unexpected (e.g., MT)

#define BUFLEN 262144 //2^18
#define DELETION 'D'
#define INSERTION 'I'
#define MISMATCH 'X'
#define INTER 'N'
#define INTRA 'R'
#define INTRA_ALIGNMENT 'Z'

typedef struct _parameters
{
	std::string assembler;
	std::string gaf;
	std::string ref_graph;
	bool asm_mode = false;
	bool debug = false;
	bool no_remap = false;
	bool skip_untagged = false;
	std::string fasta;
	std::string phase_tags;
	std::string output_path;
	std::string remap_gaf_path;
	std::string remap_gaf_path_realigned;
	std::string vcf_path;
	std::string log_path;
	std::string sample_name;
	LogFile fp_svtigs;
	LogFile fp_logs;
	int threads;
	int support;
	int dist_threshold;
	double min_map_ratio; //Used in filtering
	
	//Graphaligner parameters
	int min_alignment_score;
	double min_precise_clipping;

	_parameters() {
    }
} parameters;


class Gaf
{
private:
public:
	
	std::string query_name;
	int query_length;
	int query_start;
	int query_end;
	std::string strand;
	std::string path;
	int path_length;
	int path_start;
	int path_end;
	int residue_matches;
	int alignment_block_length;
	int mapping_quality;
	bool is_primary;
	std::string cigar;
	float aln_score;
	
	Gaf(std::string& _query_name, int& _query_length, int& _query_start, int& _query_end, std::string& _strand, std::string& _path, 
			int& _path_length, int& _path_start, int& _path_end, int& _residue_matches, int& _alignment_block_length, int& _mapping_quality, 
			bool& _is_primary, std::string& _cigar)
	{

		query_name = _query_name;
		query_length = _query_length;
		query_start = _query_start;
		query_end = _query_end;
		strand = _strand;
		path = _path;
		path_length = _path_length;
		path_start = _path_start;
		path_end = _path_end;
		residue_matches = _residue_matches;
		alignment_block_length = _alignment_block_length;
		mapping_quality = _mapping_quality;
		is_primary = _is_primary;
		cigar = _cigar;
	}
	Gaf()
	{

		query_name = "";
		query_length = -1;
		query_start = -1;
		query_end = -1;
		strand = "";
		path = "";
		path_length = -1;
		path_start = -1;
		path_end = -1;
		residue_matches = -1;
		alignment_block_length = -1;
		mapping_quality = -1;
		is_primary = false;
		cigar = "";
        aln_score = 0.0f;
	}

};
// External global logFile symbol is unused in the repo; remove to avoid exposing implementation details
// extern std::ofstream logFile; //Defined in svarp.c

unsigned hash(char *s);
void init_params(parameters** params);
void set_str(char **target, char *source);
int decompose_cigars(const std::string& cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp);
std::string exec(const std::string& command, bool return_out); 
void error(const char* const msg);
double overlap_ratio(int x_start, int x_end, int y_start, int y_end);
int parse_gaf_line(std::string& line, Gaf& gafline);
int run_and_log(const std::string& cmd, parameters& params, const std::string& label = "", int retries = 0, int backoff_seconds = 1, bool fatal = false);
std::string find_executable(const std::string &progname, const std::vector<std::string> &extra_dirs = {});

std::string& reverse_complement(std::string& seq);
const std::vector<std::string>& find_prev_next_nodes(std::map <std::string, std::vector<std::string>> inout_nodes, std::string node);
#endif
