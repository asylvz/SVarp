#include <iostream>
#include <filesystem>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <iterator>
#include "remap.h"
#include "common.h"

int read_remappings(std::map<std::string, gfaNode*> gfa)
{
	int secondary = 0, primary = 0, line_count = 0, sv_count = 0;
	long total_sv_length = 0;

	std::cout<<"Reading remappings"<<std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;
	std::vector<int> cigarLen;
	std::vector<char> cigarOp;
	std::set<std::string> distinct_read;

	std::string cwd = std::filesystem::current_path().string();
	std::string remap_output_path = cwd + "/" + REMAP_OUTPUT;
	std::ifstream fp(remap_output_path);
	
	if(!fp.good())
	{
        std::cerr << "Error opening '"<<remap_output_path<< std::endl;
        return RETURN_ERROR;
    }

	while(fp)
	{
		getline(fp, line);
		line_count++;

		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();
		
		while(getline(s, tmp_str, '\t'))
        	tokens.push_back(tmp_str);
		
		if(stoi(tokens[11]) < MINMAPQ)
			continue;
			
		bool isPrimary = true;	
		std::string cigar;
		for (auto& tok : tokens) 
		{
			if(strstr(tok.c_str(), "tp:A:"))
			{
				if (tok.substr(5, 6) != "P")
				{
					isPrimary = false;
					secondary++;
					break;
				}
				primary++;
			}
			else if(strstr(tok.c_str(), "cg:Z:"))
			{
				cigarLen.clear();
				cigarOp.clear();

				cigar = tok.substr(5);
				int cigar_cnt = decompose_cigars(cigar, cigarLen, cigarOp);
				int ref_pos = 0;
				for (int c = 0; c < cigar_cnt; c++)
				{
					if (cigarOp[c] == INSERTION && cigarLen[c] > MINSVSIZE)
					{
						std::cout<<tokens[0]<<" "<<cigarLen[c]<<"\n";
						total_sv_length += cigarLen[c];
						sv_count++;
						distinct_read.insert(tokens[0]);
					}
					if (cigarOp[c] != 'I')
						ref_pos += cigarLen[c];
				}
			}
   		}

		if(!isPrimary)
			continue;
	}
	std::cout<<"\n--->there are "<<primary<<" primary and "<<line_count-primary<<" secondary mappings\n";
	std::cout<<"--->there are "<<sv_count<< " final SVs with average "<<(total_sv_length/sv_count)<<"\n\n"<<distinct_read.size();

	return RETURN_SUCCESS;
}


int remap_assemblies(parameters* params)
{

	std::cout<<"\nRemapping assemblies using minigraph..."<<std::endl;		
	std::string cwd = std::filesystem::current_path().string();
	std::string fasta_file_path = cwd + "/log/" + FASTA_OUTPUT;
	std::string remap_output_path = cwd + "/log/" + REMAP_OUTPUT;

	std::string minigraph_cmd = "minigraph -cx lr " + params->ref_graph + " " + fasta_file_path + " -t 16 --vc > " + remap_output_path;
	
	system(minigraph_cmd.c_str());
	
	std::cout<<"--->output written to "<<remap_output_path <<std::endl;		

	return RETURN_SUCCESS;
}
