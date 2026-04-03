#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <iterator>
#include <zlib.h>
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

		auto git = gfa.find(node);
		if (git == gfa.end())
			continue;

		gfaNode* gn = git->second;
		Contig* rc = ref[gn->contig];
		Contig* overall = ref["overall"];

		if ((node_count == 1) && (p == path.size())) //the single node
		{
			rc->mapped_bases += path_end - path_start;
			rc->mapped_reads++;

			overall->mapped_reads++;
			overall->mapped_bases += path_end - path_start;

			return path_end - path_start;
		}
		else if((node_count == 1) && (p < path.size())) //first node
		{
			rc->mapped_bases += gn->len - path_start;
			rc->mapped_reads++;

			overall->mapped_reads++;
			overall->mapped_bases += gn->len - path_start;

			total_so_far += gn->len - path_start;
		}
		else if(p == path.size()) //Last node
		{
			rc->mapped_bases += total_path_length - total_so_far;
			rc->mapped_reads++;

			overall->mapped_reads++;
			overall->mapped_bases += total_path_length - total_so_far;

			return total_path_length;
		}
		else //middle node
		{
			rc->mapped_bases += gn->len;
			rc->mapped_reads++;

			overall->mapped_reads++;
			overall->mapped_bases += gn->len;

			total_so_far += gn->len;
		}
	}
	std::cerr<<"Error in contig_coverage()\n";
	return RETURN_ERROR;
}

int read_gfa(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map <std::string, std::vector<std::string>>& incoming, std::map <std::string, std::vector<std::string>>& outgoing)
{

	// Read the S lines in the GFA file
	std::cout << "\nReading the GFA file"<< std::endl;
	
	std::string line;
	std::vector <std::string> tokens;

	gzFile fp = gzopen(params.ref_graph.c_str(), "rb");
	if(!fp)
	{
        std::cerr << "Error opening '"<<params.ref_graph<< std::endl;
        return RETURN_ERROR;
    }

	std::map<std::string, std::vector<std::string>>::iterator it;

	//For overall coverage
	Contig *c = new Contig();
	ref.insert(std::pair<std::string, Contig*>("overall", c));

	const int CHUNK = 1048576;
	char buf[CHUNK];
	line.clear();
	while(gzgets(fp, buf, CHUNK) != nullptr)
	{
		line.append(buf);
		// If line doesn't end with newline, gzgets hit the buffer limit — keep reading
		if (line.empty() || (line.back() != '\n'))
			continue;
		line.pop_back();
		
		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();
		
		while(getline(s, tmp_str, '\t'))
        		tokens.push_back(tmp_str);
		
		if (tokens.empty())
		{
			line.clear();
			continue;
		}

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
			line.clear();
			continue;
		}

		if (tokens.size() < 6)
		{
			line.clear();
			continue;
		}

		gfaNode *g = new gfaNode();
		g->name = tokens[1];
		g->sequence = tokens[2];

		try {
			g->len = stoi(tokens[3].substr(5));
			g->offset = stoi(tokens[5].substr(5));
		} catch (const std::invalid_argument&) {
			std::cerr << "[read_gfa] Invalid numeric value in GFA line: " << line << std::endl;
			delete g;
			line.clear();
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
		line.clear();
	}
	gzclose(fp);
	return RETURN_SUCCESS;
}

