#include <iostream>
#include <map>
#include "cmdline.h"
#include "reference.h"
#include "alignment.h"
#include "assembly.h"
#include "phasing.h"
#include "remap.h"
//#include "asm.h"


int main(int argc, char** argv)
{
	std::map <std::string, Variant*> tmp_var;
	std::map <std::string, gfaNode*> gfa;
	std::map <std::string, std::vector<std::string>> incoming, outgoing;
	std::set <std::string> unmapped_reads;

	std::map <std::string, phase*> phased_reads;
	std::map <std::string, std::vector<Svtig*>> vars;
	std::map <std::string, Contig*> depth;
	std::map <std::string, FinalSvtig*> final_svtigs;
	
	parameters params = parameters();
	if (parse_command_line(argc, argv, params) != RETURN_SUCCESS)
		return RETURN_ERROR;
	
	init_logs(params);
	if (read_gfa(params, depth, gfa, incoming, outgoing) != RETURN_SUCCESS)
		return RETURN_ERROR;	
		
	/*if(params.asm_mode)
	{
		if (read_alignments_asm(params, depth, gfa, unmapped_reads) != RETURN_SUCCESS)
			return RETURN_ERROR;
	}
	*/
	if (read_alignments(params, depth, gfa, tmp_var, unmapped_reads) != RETURN_SUCCESS)
		return RETURN_ERROR;
		
	merge_svs(params, gfa, tmp_var, vars, incoming, outgoing);
	
	//Read the TSV file and phase the reads
	if (!(params.phase_tags).empty())
	{
		std::cout<<"Phasing"<<std::endl;	
		if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
			phase_svs(phased_reads, vars);
	}
	Assembly s;
	s.run_assembly(params, depth, vars, unmapped_reads, final_svtigs);

	//Filter SVtigs by remapping to the graph
	filter_svtigs(params, gfa, final_svtigs);
	

	std::cout<<"\nThank you for using SVarp... Tschüs, güle güle, adios, bye...\n" <<std::endl;
	
	params.fp_logs.close();

	return RETURN_SUCCESS;
}

