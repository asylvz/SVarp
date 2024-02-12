#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <iterator>
#include "reference.h"
#include "common.h"


int contig_coverage(std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, Gaf& line)
{
	char *path_copy = strdup(line.path.c_str());
	
	char *mytoken = strtok(path_copy,"><");
	
	int node_count = 0, total_so_far = 0;
	int path_start = line.path_start;
	int path_end = line.path_end;
	int total_path_length = path_end - path_start;
	while(mytoken) 
	{
		node_count += 1;
		std::string node = mytoken;

		mytoken = strtok(NULL, "><");
		
		if ((node_count == 1) && (!mytoken)) //the single node
		{

			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += path_end - path_start;	
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += path_end - path_start;
		
			free(path_copy);
			return path_end - path_start;
		}
		else if((node_count == 1) && (mytoken)) //first node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += gfa[node]->len - path_start;	
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += gfa[node]->len - path_start;
			
			total_so_far += gfa[node]->len - path_start;
		}
		else if(!mytoken) //Last node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += total_path_length - total_so_far;
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += total_path_length - total_so_far;
			
			free(path_copy);
			return total_path_length;
			//if(total_path_length - total_so_far < 0)
			//	std::cout<< "middles: "<<middles <<" " <<total_path_length<<" "<<total_so_far<<std::endl;
		}
		else //middle node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += gfa[node]->len;		
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += gfa[node]->len;		

			total_so_far += gfa[node]->len;
		}
	}
	std::cout<<"Error in contig_coverage()\n";
	free(path_copy);
	return -1;
}

int read_gfa(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map <std::string, std::vector<std::string>>& incoming, std::map <std::string, std::vector<std::string>>& outgoing)
{

	// Read the S lines in the GFA file
	std::cout << "\nReading the GFA file"<< std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;	
	
	std::ifstream fp(params.ref_graph);
	if(!fp.good())
	{
        std::cerr << "Error opening '"<<params.ref_graph<< std::endl;
        return RETURN_ERROR;
    }

	std::map<std::string, std::vector<std::string>>::iterator it;
	
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
		{
			if (tokens[0] == "L")
			{
				//std::cout<<tokens[1]<<" " << tokens[3]<<"\n";
				//Add incoming
				
				it = incoming.find(tokens[3]);
				if (it != incoming.end())
				{	
					//if(tokens[3] == "s99992")
					//	std::cout<<tokens[1]<<"\n";
					it->second.push_back(tokens[1]);
				}
				else
				{
					std::vector<std::string> v;
					v.clear();
					v.push_back(tokens[1]);
					incoming.insert(std::pair<std::string, std::vector<std::string>>(tokens[3], v));
				}
				//Add outgoing	
				it = outgoing.find(tokens[1]);
				if (it != outgoing.end())
					it->second.push_back(tokens[3]);
				else
				{
					std::vector<std::string> v;
					v.clear();
					v.push_back(tokens[3]);
					outgoing.insert(std::pair<std::string, std::vector<std::string>>(tokens[1], v));
				}
			}
			continue;
		}
		
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

