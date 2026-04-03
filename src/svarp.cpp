#include <iostream>
#include <fstream>
#include <map>
#include <chrono>
#include "cmdline.h"
#include "reference.h"
#include "alignment.h"
#include "phasing.h"
#include "remap.h"
#include "assembly.h"


int main(int argc, char** argv)
{
	std::cout << std::unitbuf; // flush stdout after every write

	auto t_start = std::chrono::steady_clock::now();

	std::map <std::string, Variant*> tmp_var;
	std::map <std::string, gfaNode*> gfa;
	std::map <std::string, std::vector<std::string>> incoming, outgoing;
	std::set <std::string> unmapped_reads;

	std::map <std::string, phase*> phased_reads;
	std::map <std::string, std::vector<SVCluster*>> vars;
	std::map <std::string, Contig*> depth;
	std::map <std::string, SVtig*> svtigs;

	parameters params = parameters();
	if (parse_command_line(argc, argv, params) != RETURN_SUCCESS)
		return RETURN_ERROR;

	init_logs(params);

	log_step(params.fp_logs, "Reading reference graph");
	if (read_gfa(params, depth, gfa, incoming, outgoing) != RETURN_SUCCESS)
		return RETURN_ERROR;

	log_step(params.fp_logs, "Reading alignments");
	if (read_alignments(params, depth, gfa, tmp_var, unmapped_reads) != RETURN_SUCCESS)
		return RETURN_ERROR;

	log_step(params.fp_logs, "Merging SV signals");
	merge_svs(params, gfa, tmp_var, vars, incoming, outgoing);

	//Read the TSV file and phase the reads
	if (!(params.phase_tags).empty())
	{
		log_step(params.fp_logs, "Phasing");
		std::cout<<"Phasing"<<std::endl;
		if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
			phase_svs(phased_reads, vars);
	}

	log_step(params.fp_logs, "Assembly");
	Assembly s;
	s.run_assembly(params, depth, vars, unmapped_reads, svtigs);

	log_step(params.fp_logs, "Filtering svtigs (GraphAligner remapping)");
	//Filter SVtigs by remapping to the graph
	filter_svtigs(params, gfa, svtigs);

	auto t_end = std::chrono::steady_clock::now();
	double total_sec = std::chrono::duration<double>(t_end - t_start).count();

	// Write summary to log
	if (params.fp_logs.is_open()) {
		params.fp_logs << "\n==============================\n";
		params.fp_logs << "Run finished: " << current_timestamp() << "\n";
		params.fp_logs << "Total wall time: " << format_duration(total_sec) << "\n";
		params.fp_logs << "==============================\n";
	}

	std::cout<<"\nThank you for using SVarp... Tschüs, güle güle, adios, bye...\n" <<std::endl;

	if (params.fp_logs.is_open()) params.fp_logs.close();
	if (params.fp_asm_log.is_open()) params.fp_asm_log.close();
	if (params.fp_remap_log.is_open()) params.fp_remap_log.close();

	// Cleanup heap-allocated objects
	for (auto& p : tmp_var)
		delete p.second;
	for (auto& p : gfa)
		delete p.second;
	for (auto& p : phased_reads)
		delete p.second;
	for (auto& p : vars)
		for (auto* sv : p.second)
			delete sv;
	for (auto& p : depth)
		delete p.second;
	for (auto& p : svtigs)
		delete p.second;

	return RETURN_SUCCESS;
}

