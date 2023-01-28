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

std::map<std::string, gfaNode*> read_gfa(parameters* params, std::map <std::string, Contig*>& ref)
{	
	// Read the S lines in the GFA file
	std::cout << "Reading the GFA file"<< std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;	
	std::ifstream fp(params->ref_graph);
	std::map<std::string, gfaNode*> gfa;
	
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

	//std::cout<< ref.size()<<std::endl;
	/*std::map<std::string, Contig*>::iterator it;
	for (it=ref.begin(); it != ref.end(); ++it)
	{	
		//std::cout<< it->first<<" LEN= "<<it->second->contig_length <<std::endl;
		//it->second->coverage = (double) it->second->
	}*/	
	//cout <<ref["s10"]->contig<<endl;
	/*map<std::string, gfaNode*>::iterator itr;
    cout << "\nThe map is :"<<cnt<< "\n";
    cout << "\tKEY\tELEMENT\n";
    for (itr = ref.begin(); itr != ref.end(); ++itr) {
        cout << '\t' << itr->first << '\t' << itr->second->offset<< '\n';
    }*/
	
	return gfa;
}

