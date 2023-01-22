#ifndef __SV
#define __SV

#include <map>
#include <vector>
#include "gfa.h" 


#define DELETION 'D'
#define INSERTION 'I'

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

variant* generate_sv_node(std::map<std::string, gfaNode*>, std::string, int, int, int, int, char);


#endif
