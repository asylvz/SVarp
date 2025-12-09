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
	const std::string &path = line.path;
	
	int node_count = 0, total_so_far = 0;
	int path_start = line.path_start;
	int path_end = line.path_end;
	int total_path_length = path_end - path_start;
	
	size_t p = 0;
	while (p < path.size())
	{
		++p;  // skip strand char
		size_t q = p;
		while (q < path.size() && path[q] != '>' && path[q] != '<') ++q;
		std::string node = path.substr(p, q - p);
		p = q;
		node_count += 1;
		
		if ((node_count == 1) && (p == path.size())) //the single node
		{

			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += path_end - path_start;	
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += path_end - path_start;
		
			return path_end - path_start;
		}
		else if((node_count == 1) && (p < path.size())) //first node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += gfa[node]->len - path_start;	
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += gfa[node]->len - path_start;
			
			total_so_far += gfa[node]->len - path_start;
		}
		else if(p == path.size()) //Last node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += total_path_length - total_so_far;
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += total_path_length - total_so_far;
			
			return total_path_length;
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
	return RETURN_ERROR;
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
		
		if (tokens.empty())
			continue;
		
		if (tokens[0] != "S")
		{
			if (tokens[0] == "L" && tokens.size() >= 4)
			{
				// Add incoming: use insert with hint to avoid redundant lookups
				auto [it_in, inserted_in] = incoming.insert({tokens[3], {}});
				it_in->second.push_back(tokens[1]);
				
				// Add outgoing: same pattern
				auto [it_out, inserted_out] = outgoing.insert({tokens[1], {}});
				it_out->second.push_back(tokens[3]);
			}
			continue;
		}
		
		if (tokens.size() < 6)
			continue;  // Skip malformed S lines
		
		gfaNode *g = new gfaNode();
		g->name = tokens[1];
		g->sequence = tokens[2];
		
		try {
			g->len = stoi(tokens[3].substr(5));
			g->offset = stoi(tokens[5].substr(5));
		} catch (const std::invalid_argument&) {
			std::cerr << "[read_gfa] Invalid numeric value in GFA line: " << line << std::endl;
			delete g;
			continue;
		}
		
		g->contig = tokens[4].substr(5);

		// Insert Contig if not present, using insert return value
		auto [it_contig, inserted] = ref.insert({g->contig, nullptr});
		if (inserted)
			it_contig->second = new Contig();
		
		//Find the length of each contig	
		it_contig->second->contig_length += g->len;
		ref["overall"]->contig_length += g->len;

		gfa.insert({g->name, g});
	}
	return RETURN_SUCCESS;
}

