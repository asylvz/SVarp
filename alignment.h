#ifndef __ALIGNMENT
#define __ALIGNMENT

#include <vector>
#include "common.h"
#include "gfa.h" 
#include "sv.h" 

using namespace std;

#define MINMAPQ 20

class alignment
{
private:
public:

	string read_name;
	int start;
	int end;
	char strand;
	string node;
	string path;	
		
	alignment()
	{
	}
	alignment(string _read_name, int _start, int _end, char _strand, string _node, string _path)
	{
		read_name = _read_name;
		start = _start;
		end = _end;
		strand = _strand;
		node = _node;
		path = _path;
	}
	virtual ~alignment() {};
	void display();
};

void alignment_within_gfa(map<string, alignment*>& gaf, map<string, gfaNode*> gfa, vector <std::string> tokens);

int read_alignments(parameters *params, std::map<std::string, gfaNode*> ref, std::vector<variant*>& insertions);

#endif
