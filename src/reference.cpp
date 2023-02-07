#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <iterator>
#include "reference.h"
#include "common.h"


int read_gfa(parameters* params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa)
{	
	// Read the S lines in the GFA file
	std::cout << "\nReading the GFA file"<< std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;	
	
	std::ifstream fp(params->ref_graph);
	if(!fp.good())
	{
        std::cerr << "Error opening '"<<params->ref_graph<< std::endl;
        return RETURN_ERROR;
    }

	//For overall coverage
	Contig *c = new Contig();
	ref.insert(std::pair<std::string, Contig*>("overall", c));
	
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

		std::map<std::string, Contig*>::iterator it = ref.find(g->contig);
		if (it == ref.end())
		{
			Contig *c = new Contig();
			ref.insert(std::pair<std::string, Contig*>(g->contig, c));
		}
		
		//Find the length of each contig	
		ref[g->contig]->contig_length += g->len;
		ref["overall"]->contig_length += g->len;

		gfa.insert(std::pair<std::string, gfaNode*>(g->name, g));
	}
	return RETURN_SUCCESS;
}

