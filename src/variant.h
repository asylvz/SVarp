#ifndef __SV
#define __SV

#include <map>
#include <unordered_map>
#include <vector>
#include "reference.h" 
#include "common.h"

class SVCluster
{
private:
public:

	std::string name;
	char sv_type; 
	std::string contig;
	std::string node;
	std::set <std::string> reads_h1;
	std::set <std::string> reads_h2;
	std::set <std::string> reads_untagged;
	int start_pos = 0;
	int end_pos = 0; //largest position or deletion end among the signals
	int ref_pos = 0;
	bool phased;
	bool filter = false;
	std::string	path;
	int rank = -1; //SR of the node, -1 without tag
	SVCluster()
	{
	}
	virtual ~SVCluster() {};
};


class SVtig
{
private:
public:

	std::string name; //Format e.g., H1-node_name_start_pos_in_node
	int pos; //This is the reference poisiton, not the node position
	//int coverage; //Number of supporting reads
	std::string contig;
	std::set <std::string> reads; //Read names that support this SVtig
	bool output = false; //Whether to output after remapping (filtered if false)
	std::string seq; //svtig genomic sequence
	std::string remap_path; //Graph path this svtig aligns to, from remapping
	double map_ratio = -1; //Fraction of the svtig covered by graph alignments; -1 if not remapped
	int max_gap = 0; //largest uncovered stretch between graph alignments (bp)
	int max_indel = 0; //largest indel inside a graph alignment (bp)
	//Non-reference nodes on remap_path, "first-last:contig" joined by ';'
	std::string alt_nodes;
};


class Variant
{
private:
public:
		
	char sv_type = 0; // For intra: "DELETION" or "INSERTION"
	char type = 0; 	// "INTER" or "INTRA"
	int sv_size = 0;
	int pos_in_node = 0;
	int pos_in_node_end = 0;
	int pos_in_ref = 0;
	int pos_in_ref_end = 0;
	int node_count = 0;
	std::string contig;
	std::string genotype;
	std::set <std::string> reads_h1;	
	std::set <std::string> reads_h2;
	std::set <std::string> reads_untagged;
	bool phased;
	std::string node;
	char node_strand;
	std::string path;
	bool duplicate = false;

	Variant()
	{
	}

	virtual ~Variant() {};
};


Variant* generate_sv_node(std::map<std::string, gfaNode*>& gfa, Gaf& line, const int base_pos, int var_len, char sv_type);
int merge_svs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::map<std::string, std::vector<SVCluster*>>& final_svtigs, EdgeMap& incoming, EdgeMap& outgoing);

int merge_neighbor_nodes(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, std::vector<SVCluster*>>& init_svtigs, EdgeMap& incoming, EdgeMap& outgoing);

int mapping_start_end(std::map<std::string, gfaNode*>& gfa, Gaf& line, std::map<std::string, Variant*>& variations_inter, int min_end = MIN_READ_START_END_WINDOW);
int find_deletions(parameters* params, std::map<std::string, std::vector<SVCluster*>> deletions);

#endif
