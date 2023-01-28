#ifndef __REFERENCE
#define __REFERENCE

#include <string>
#include <set>
#include "common.h"


typedef struct _contig
{
	//std::string contig_name;
	//long read_count;
	long mapped_bases;
	long contig_length;
	double coverage;

	_contig() {
		mapped_bases = 0;
		contig_length = 0;
		coverage = 0;
    }
} Contig;


class gfaNode
{
private:
public:

	std::string name;
	std::string sequence;
	int len;
	int offset;
	std::string contig;	
	
	gfaNode()
	{
	}
	gfaNode(std::string _name, std::string _sequence, int _len, std::string _contig, int _offset)
	{
		name = _name;
		sequence = _sequence;
		len = _len;
		contig = _contig;
		offset = _offset;
	}
	virtual ~gfaNode() {};
};


int read_gfa(parameters* params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa);
#endif
