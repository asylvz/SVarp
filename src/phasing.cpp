#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "phasing.h"
#include "variant.h"


int read_phase_file(parameters& params, std::map<std::string, phase*>& phased_reads)
{
	if (params.phase_tags.empty())
		return RETURN_ERROR;
	
	std::cout<<"--> reading .tsv file"<<std::endl;
	std::ifstream fp(params.phase_tags);
	if (!fp.is_open())
	{
		std::cerr << "Error: could not open phase file: " << params.phase_tags << std::endl;
		return RETURN_ERROR;
	}
	std::vector <std::string> tokens;
	std::string line;	
	
	while(getline(fp, line))
	{
		if (line.empty() || line[0] == '#')
			continue;

		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();

		while(getline(s, tmp_str, '\t'))
        	tokens.push_back(tmp_str);

		if (tokens.size() < 4)
			continue;

		phase *temp = new phase;
    	temp->read_name = tokens[0];
    	temp->haplotype = tokens[1];
    	temp->phase_set = tokens[2];
    	temp->contig = tokens[3];

		phased_reads.insert(std::pair<std::string, phase*>(temp->read_name, temp));
	}
	return RETURN_SUCCESS;
}

void phase_svs(const std::map<std::string, phase*>& phased_reads, std::map<std::string, std::vector<SVCluster*>>& vars)
{
	std::map<std::string, std::vector<SVCluster*>>::iterator itr;

	// A read supports several clusters, so tally names rather than occurrences.
	std::set<std::string> h1_reads, h2_reads, untagged_reads;

	std::cout<<"--> reading SVs to phase"<<std::endl;
	for (itr=vars.begin(); itr != vars.end(); ++itr)
	{
		for (auto &sv: itr->second) 
		{
			for (auto &read: sv->reads_untagged)
			{
				auto pit = phased_reads.find(read);
				if (pit == phased_reads.end())
				{
					// Absent from the haplotag file, so untagged like any other
					// read the phaser could not place.
					untagged_reads.insert(read);
					continue;
				}

				if(pit->second->haplotype == "none" || pit->second->phase_set == "none")
				{
					untagged_reads.insert(read);
					continue;
				}

				if (pit->second->haplotype == "H1")
				{
					sv->reads_h1.insert(read);
					h1_reads.insert(read);
				}
				else if (pit->second->haplotype == "H2")
				{
					sv->reads_h2.insert(read);
					h2_reads.insert(read);
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

	int svtig_h1 = 0, svtig_h2 = 0, svtig_untagged = 0;
	for (itr=vars.begin(); itr != vars.end(); ++itr)
	{
		for (auto &sv: itr->second) 
		{
			if (sv->reads_h1.size() > 0)
				svtig_h1++;
			if (sv->reads_h2.size() > 0)
				svtig_h2++;
			if (sv->reads_untagged.size() > 0)
				svtig_untagged++;
		}

	}
	std::cout<<"--> " <<svtig_h1+svtig_h2+svtig_untagged<<" phased read clusters ("<<svtig_h1<<" H1, " <<svtig_h2 << " H2, "<< svtig_untagged <<" untagged based on read-support > 0)\n";
	std::cout<<"--> " <<h1_reads.size()<<" reads in haplotype 1, "<<h2_reads.size()<< " in haplotype 2 and "<<untagged_reads.size()<< " are untagged in total\n";
}


