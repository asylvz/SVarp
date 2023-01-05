#ifndef __GFA
#define __GFA

#include <string>

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

#endif
