#ifndef __ASSEMBLY
#define __ASSEMBLY

#include <map>
#include "common.h"
#include "reference.h"
#include "variant.h"
#include <htslib/faidx.h>

int remap_assemblies(parameters& params);


class Assembly
{
private:
		std::string paf_H1_svtig1_path;

public:
		int filter_hicov = 0;
		int filter_lowcov = 0;
		int filter_support = 0;
		int unassembled_cnt = 0;
		std::set <std::string> raw_svtigs;

		void run_assembly(parameters& params, std::map <std::string, Contig*>& depth, std::map<std::string, std::vector<SVCluster*>>& vars, std::set <std::string>& unmapped, std::map <std::string, SVtig*>& final_svtigs);
		void generate_fasta_file(parameters& params, faidx_t*& fasta_index, std::set <std::string>& reads, std::string file_path);
		int write_svtigs(std::string& f_path, const std::string& f_name, int pos, std::string& contig, int coverage, std::ofstream& fp_write);
		int merge_svtigs(parameters& params);
		int final_assembly(parameters& params, faidx_t*& fasta_index, std::set <std::string>& read_set, std::string& svtig_name, double& contig_depth, SVCluster*& sv, std::map <std::string, SVtig*>& final_svtigs);
		int assemble_clusters(parameters& params, faidx_t*& fasta_index, std::vector<SVCluster*>& sv_cluster, std::map <std::string, Contig*>& depth, std::map <std::string, SVtig*>& final_svtigs);

};

#endif
