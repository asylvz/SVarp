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
	
	std::cout<<"--->Reading .tsv file"<<std::endl;	
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
	std::map<std::string, phase*>::iterator itr2;

	int phased = 0, not_phased = 0;
	
	std::cout<<"--->Checking SVs"<<std::endl;	
	for (itr=insertions.begin(); itr != insertions.end(); ++itr)
	{	
		for (auto &sv: itr->second) 
		{
			//first check if it can be phased
			std::string tmp_phase;
			int tmp_cnt = 0;
			sv->phased = false;
				
			//std::cout<<"First pass "<<sv->reads_h1.size()<<"\n";
			for (auto &read: sv->reads_h1)
			{	
				std::map<std::string, phase*>::iterator it = phased_reads.find(read);
				if (it == phased_reads.end())
					continue;

				//std::cout<<phased_reads[read]->haplotype<<" "<< phased_reads[read]->phase_set <<std::endl;
				if(phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
				{
			   		sv->phased = false;
					break;
				}

				if (tmp_cnt == 0)
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
				tmp_cnt++;
			}
		
			//phasing
			if (sv->phased == true)
			{
				//std::cout<<"Second pass\n";
				for (auto &read: sv->reads_h1)
				{
					//std::cout<<phased_reads[read]->haplotype<<" "<< phased_reads[read]->phase_set <<std::endl;
					if(phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
			   		continue;
					if (phased_reads[read]->haplotype == "H2")
						sv->reads_h2.insert(read);
				}
				//Erase reads added to the second set, from the first set
				for (auto &read: sv->reads_h2)
					sv->reads_h1.erase(read);

				phased++;
			}
			else
				not_phased++;
			//std::cout<<"h1: "<<sv->reads_h1.size()<< "\th2: "<<sv->reads_h2.size()<<std::endl;
		}
	}
	std::cout<<"--->"<<phased<<" SVs could be phased, but "<<not_phased<< " could not...\n";
}


