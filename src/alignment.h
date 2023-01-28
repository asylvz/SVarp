#ifndef __ALIGNMENT
#define __ALIGNMENT

#include <vector>
#include "common.h"
#include <map>
#include "reference.h" 
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


void find_supporting_reads(std::map<std::string, gfaNode*> ref, std::multimap<std::string, alignment*> aln, std::set<std::string> contigs, std::multimap<std::string, variant*>& insertions);

void alignment_within_gfa(multimap<string, alignment*>& gaf, map<string, gfaNode*> gfa, vector <std::string> tokens);

int read_alignments(parameters *params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*> gfa, std::map<std::string, variant*>& insertions);
#endif
