#ifndef __REFERENCE
#define __REFERENCE

#include <string>
#include <set>
#include <map>
#include <vector>
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
	int len = 0;
	int offset = 0;
	std::string contig;
	//rGFA SR: 0 is the reference, >0 an allele, -1 no tag
	int rank = -1;

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
//GFA link seen from one node: the other node and the L-line orientations (from, to).
//'+' leaves a node through its forward end and enters the next through its forward start.
struct Edge
{
	std::string node;
	char from = '+';
	char to = '+';
};
typedef std::map<std::string, std::vector<Edge>> EdgeMap;

int read_gfa(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, EdgeMap& incoming, EdgeMap& outgoing);


#endif
