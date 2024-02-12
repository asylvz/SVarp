#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "phasing.h"


int read_phase_file(parameters& params, std::map<std::string, phase*>& phased_reads)
{
	if (params.phase_tags.empty())
		return RETURN_ERROR;
	
	std::cout<<"--->reading .tsv file"<<std::endl;	
	std::ifstream fp(params.phase_tags);
	std::vector <std::string> tokens; 
	std::string line;	
	
	while(getline(fp, line))
	{
		if (line[0] == '#')
			continue;
		
		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();
		
		while(getline(s, tmp_str, '\t'))
        	tokens.push_back(tmp_str);	
		
		phase *temp = new phase;
    	temp->read_name = tokens[0];
    	temp->haplotype = tokens[1];
    	temp->phase_set = tokens[2];
    	temp->contig = tokens[3];

		phased_reads.insert(std::pair<std::string, phase*>(temp->read_name, temp));
	}
	return RETURN_SUCCESS;
}

void phase_svs(std::map<std::string, phase*> phased_reads, std::map<std::string, std::vector<Svtig*>>& vars)
{
	std::map<std::string, std::vector<Svtig*>>::iterator itr;

	int h1_cnt = 0, h2_cnt = 0, none_cnt = 0;

	std::cout<<"--->reading SVs to phase"<<std::endl;	
	for (itr=vars.begin(); itr != vars.end(); ++itr)
	{	
		for (auto &sv: itr->second) 
		{
			for (auto &read: sv->reads_untagged)
			{
				if (phased_reads.find(read) == phased_reads.end())
				{
					//std::cout <<read<< " not found in your .tsv file...\n";
					continue;
				}	
				//std::cout<<phased_reads[read]->haplotype<<" "<< phased_reads[read]->phase_set <<std::endl;
				if(phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
				{
					//none_reads.insert(read);
					none_cnt++;
					continue;
				}

				if (phased_reads[read]->haplotype == "H1")
				{
					sv->reads_h1.insert(read);
					h1_cnt++;
				}
				else if (phased_reads[read]->haplotype == "H2")
				{
					sv->reads_h2.insert(read);
					h2_cnt++;
				}
			}
			//Erase reads added to the second set, from the first set
			for (auto &read: sv->reads_h1)
				sv->reads_untagged.erase(read);
			for (auto &read: sv->reads_h2)
				sv->reads_untagged.erase(read);
					
			//std::cout<<"h1: "<<sv->reads_h1.size()<< "\th2: "<<sv->reads_h2.size()<<std::endl;
		}
	}
	//std::cout<<"--->"<<phased<<" SVtigs phased - "<<not_phased<< " not...\n";	
	std::cout<<"--->"<<h1_cnt<<" reads are haplotype 1, "<<h2_cnt<< " haplotype 2 and "<<none_cnt<< " are untagged...\n";	
}


