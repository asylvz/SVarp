#ifndef __SV
#define __SV

#include <map>
#include <unordered_map>
#include <vector>
#include "reference.h" 
#include "common.h"

class Svtig
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
	int start_pos;
	int end_pos;
	int ref_pos;
	bool phased;
	bool filter = false;
	std::string	path;
	Svtig()
	{
	}
	virtual ~Svtig() {};
};


class FinalSvtig
{
private:
public:

	std::string name; //Format e.g., H1-node_name:start_pos_in_node
	int pos; //This is the reference poisiton, not the node position
	//int coverage; //Number of supporting reads
	std::string contig;
	std::set <std::string> reads; //Read names that support this SVtig
	bool output = false; //Whether to output after remapping (filtered if false)
};


class Variant
{
private:
public:
		
	char sv_type; // For intra: "DELETION" or "INSERTION"
	char type; 	// "INTER" or "INTRA"
	int sv_size;
	int pos_in_node;
	int pos_in_node_end;
	int pos_in_ref;
	int pos_in_ref_end;
	int node_count;
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
int refine_svs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::map<std::string, std::vector<Svtig*>>& final_svtigs, std::map <std::string, std::vector<std::string>>& incoming, std::map <std::string, std::vector<std::string>>& outgoing);

int mapping_start_end(std::map<std::string, gfaNode*>& gfa, Gaf& line, std::map<std::string, Variant*>& variations_inter);
int find_deletions(parameters* params, std::map<std::string, std::vector<Svtig*>> deletions);

#endif
