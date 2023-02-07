#ifndef __SV
#define __SV

#include <map>
#include <vector>
#include "reference.h" 


class svtig
{
private:
public:

	char sv_type; 
	std::string contig;
	std::set <std::string> reads_h1;	
	std::set <std::string> reads_h2;
	std::vector <int> breakpoints; 		//The SV breakpoints that contribute to this svtig
	int sv_size;
	int start_pos;
	int end_pos;
	bool phased;
	
	svtig()
	{
	}
	virtual ~svtig() {};
};


class variant
{
private:
public:

	char sv_type; 
	int sv_size; 
	int ref_size; 
	int ref_start; 
	int ref_end;
	std::string contig;
	std::set <std::string> reads_h1;	
	std::set <std::string> reads_h2;
	bool phased;
	std::string node;
	char node_strand;

	variant()
	{
	}

	variant(char _sv_type, int _sv_size, int _ref_size, int _ref_start, int _ref_end, std::string _contig, std::string _node, char _node_strand)
	{
		sv_type = _sv_type;
		sv_size = _sv_size;
		ref_size = _ref_size;
		ref_start = _ref_start;
		ref_end = _ref_end;
		contig = _contig;
		node = _node;
		node_strand = _node_strand;
	}
	virtual ~variant() {};
};


variant* generate_sv_node(std::map<std::string, gfaNode*>& gfa, std::string path, int path_start, int path_end, int ref_pos, int cigar_len, char sv_type);

int refine_svs(std::map<std::string, variant*> initial_variations, std::map<std::string, std::vector<svtig*>>& final_ins, std::map<std::string, std::vector<svtig*>>& final_del);


int find_deletions(parameters* params, std::map<std::string, std::vector<svtig*>> deletions);

#endif
