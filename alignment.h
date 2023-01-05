#ifndef __ALIGNMENT
#define __ALIGNMENT

#include "common.h"

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

int read_alignments(parameters *params);

#endif
