#ifndef __GFA
#define __GFA

#include <string>
#include <set>
#include "common.h"

using namespace std;

class gfaNode
{
private:
public:

	string name;
	string sequence;
	int len;
	int offset;
	string contig;	
	
	gfaNode()
	{
	}
	gfaNode(string _name, string _sequence, int _len, string _contig, int _offset)
	{
		name = _name;
		sequence = _sequence;
		len = _len;
		contig = _contig;
		offset = _offset;
	}
	virtual ~gfaNode() {};
};

std::map<std::string, gfaNode*> read_gfa(parameters* params, std::set<std::string>& contigs);

#endif
