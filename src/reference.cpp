#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <iterator>
#include "gfa.h"
#include "common.h"

std::map<std::string, gfaNode*> read_gfa(parameters* params, std::set<std::string>& contigs)
{	
	// Read the S lines in the GFA file
	std::cout << "Reading the GFA file"<< std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;	
	std::ifstream fp(params->ref_graph);
	map<string, gfaNode*> ref;
	
	while(fp)
	{
		getline(fp, line);
			
		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();
		
		while(getline(s, tmp_str, '\t'))
        	tokens.push_back(tmp_str);
		
		if (tokens[0] != "S")
			continue;
		
		gfaNode *g = new gfaNode();
		g->name = tokens[1];
		g->sequence = tokens[2];
		g->len = stoi(tokens[3].substr(5));
		g->contig = tokens[4].substr(5);
		g->offset = stoi(tokens[5].substr(5));
		
		contigs.insert(g->contig);
		ref.insert(std::pair<std::string, gfaNode*>(g->name, g));
	}
	//cout <<ref["s10"]->contig<<endl;
	/*map<std::string, gfaNode*>::iterator itr;
    cout << "\nThe map is :"<<cnt<< "\n";
    cout << "\tKEY\tELEMENT\n";
    for (itr = ref.begin(); itr != ref.end(); ++itr) {
        cout << '\t' << itr->first << '\t' << itr->second->offset<< '\n';
    }*/
	
	return ref;
}

