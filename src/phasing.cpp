#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "phasing.h"


int read_phase_file(parameters *params, std::map<std::string, phase*>& phased_reads)
{
	if (params->phase_tags.empty())
		return RETURN_ERROR;
	
	std::cout<<"Reading .tsv file"<<std::endl;	
	std::ifstream fp(params->phase_tags);
	std::vector <std::string> tokens; 
	std::string line;	
	
	while(fp)
	{
		getline(fp, line);

		if (line[0] == '#')
			continue;
		
		//std::cout<<line<<std::endl;
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

void phase_svs(std::map<std::string, phase*> phased_reads, std::map<std::string, variant*>& insertions)
{
	std::map<std::string, variant*>::iterator itr;
	std::map<std::string, phase*>::iterator itr2;

	std::cout<<"Phasing"<<std::endl;	
	for (itr=insertions.begin(); itr != insertions.end(); ++itr)
	{
		if (itr->second->reads_h1.size() <= MIN_READ_SUPPORT)
			continue;
		
		//first check if it can be phased
		std::string tmp_phase;
		int tmp_cnt = 0;
		itr->second->phased = true;
		for (auto &read: itr->second->reads_h1)
		{
			//std::cout<<phased_reads[read]->haplotype<<" "<< phased_reads[read]->phase_set <<std::endl;
			if(phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
			   continue;

			if (tmp_cnt == 0)
				tmp_phase = phased_reads[read]->phase_set;
			else
				if(tmp_phase != phased_reads[read]->phase_set)
				{
					itr->second->phased = false;
					break;
				}	
			tmp_cnt++;
		}
		//phasing
		if (itr->second->phased == true)
		{
			//std::cout<<"Second pass\n";
			for (auto &read: itr->second->reads_h1)
			{
				//std::cout<<phased_reads[read]->haplotype<<" "<< phased_reads[read]->phase_set <<std::endl;
				if(phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
			   	continue;
				if (phased_reads[read]->haplotype == "H2")
					itr->second->reads_h2.insert(read);
			}
			//Erase reads added to the second set, from the first set
			for (auto &read: itr->second->reads_h2)
				itr->second->reads_h1.erase(read);
		}
		//std::cout<<"h1: "<<itr->second->reads_h1.size()<< "\th2: "<<itr->second->reads_h2.size()<<std::endl;
	}
}


