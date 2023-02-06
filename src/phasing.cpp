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
	
	std::cout<<"--->reading .tsv file"<<std::endl;	
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

void phase_svs(std::map<std::string, phase*> phased_reads, std::map<std::string, std::vector<svtig*>>& insertions)
{
	std::map<std::string, std::vector<svtig*>>::iterator itr;
	std::set <std::string> none_reads;

	int phased = 0, not_phased = 0, skipped_reads = 0, read_not_found = 0;
	int phased_read_count;

	std::cout<<"--->reading SVs to phase"<<std::endl;	
	for (itr=insertions.begin(); itr != insertions.end(); ++itr)
	{	
		for (auto &sv: itr->second) 
		{
			//first check if it can be phased
			std::string tmp_phase;
			phased_read_count = 0;
			sv->phased = false;
				
			//std::cout<<"First pass "<<sv->reads_h1.size()<<"\n";
			for (auto &read: sv->reads_h1)
			{	
				std::map<std::string, phase*>::iterator it = phased_reads.find(read);
				if (it == phased_reads.end())
				{
					read_not_found++;
					continue;
				}

				//std::cout<<phased_reads[read]->haplotype<<" "<< phased_reads[read]->phase_set <<std::endl;
				if(phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
				{
					skipped_reads++;
					continue;
				}

				if (phased_read_count == 0)
				{
					tmp_phase = phased_reads[read]->phase_set;
					sv->phased = true;
				}
				else
					if(tmp_phase != phased_reads[read]->phase_set)
					{
						sv->phased = false;
						break;
					}	
				phased_read_count++;
			}
		
			//phasing	
			//std::cout<<"Second pass\n";
			if (sv->phased == true)
			{
				for (auto &read: sv->reads_h1)
				{
					//std::cout<<phased_reads[read]->haplotype<<" "<< phased_reads[read]->phase_set <<std::endl;
					if(phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
					{
						none_reads.insert(read);
						continue;
					}
					if (phased_reads[read]->haplotype == "H2")
						sv->reads_h2.insert(read);
				}
				//Erase reads added to the second set, from the first set
				for (auto &read: sv->reads_h2)
					sv->reads_h1.erase(read);
				
				//Erase none reads from the first set
				for (auto &read: none_reads)
					sv->reads_h1.erase(read);
				
				none_reads.clear();

				phased++;
				//std::cout<<"h1: "<<sv->reads_h1.size()<< "\th2: "<<sv->reads_h2.size()<<std::endl;
			}
			else
				not_phased++;
			
		}
	}
	std::cout<<"--->"<<phased<<" SVs phased - "<<not_phased<< " not...\n";
}


