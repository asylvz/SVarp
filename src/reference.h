#ifndef __REFERENCE
#define __REFERENCE

#include <string>
#include <set>
#include "common.h"

typedef struct _contig
{
public:
	long mapped_bases;
	long mapped_reads;
	long contig_length;
	double coverage;

	_contig() {
		mapped_bases = 0;
		mapped_reads = 0;
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
	gfaNode(const std::string& _name, const std::string& _sequence, int _len, const std::string& _contig, int _offset)
	{
		name = _name;
		sequence = _sequence;
		len = _len;
		contig = _contig;
		offset = _offset;
	}
	virtual ~gfaNode() {};
};


int contig_coverage(std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, Gaf& line);
int read_gfa(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map <std::string, std::vector<std::string>>& incoming, std::map <std::string, std::vector<std::string>>& outgoing);


#endif
