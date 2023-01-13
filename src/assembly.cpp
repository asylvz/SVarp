#include <iostream>
#include <fstream>
#include "assembly.h"


void generate_fastq_file(parameters* params, variant* var)
{
	std::cout<<"Generating fastq file"<<std::endl;

	std::ifstream fp_read(params->ref_graph);
	std::ofstream fp_write("log/tmp.fasta");
}


void run_assembly(parameters* params, std::multimap<std::string, variant*>& insertions)
{	
	std::cout<<"Assembly..."<<std::endl;


	std::multimap<std::string, variant*>::iterator itr;


	for (itr=insertions.begin(); itr != insertions.end(); ++itr)
	{
		//std::cout << itr->first << '\t'<< itr->second->ref_start << '\t'<< itr->second->ref_end << " ("<<itr->second->reads.size()<<" read support)\n";
		//for (auto &a: itr->second->reads)			
			//std::cout << '\t' << a << '\n';
		if (itr->second->reads.size() == 0)
			std::cout<< "ALARM - there is no read supporting the SV ("<<itr->second->ref_start<<" - "<<itr->second->ref_end<<") in "<< itr->second->contig <<std::endl;

		if (itr->second->reads.size() > 5)
		{
			generate_fastq_file(params, itr->second);	
		}
	}	

}

