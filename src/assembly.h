#ifndef __ASSEMBLY
#define __ASSEMBLY

#include <map>
#include <mutex>
#include "common.h"
#include "reference.h"
#include "variant.h"
struct faidx_t;

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
		int no_contig_cnt = 0;
		std::set <std::string> raw_svtigs;
		std::mutex mtx; //final_svtigs, counters and the svtig output file

		static double cluster_depth(const SVCluster* sv, std::map <std::string, Contig*>& depth);
		void run_assembly(parameters& params, std::map <std::string, Contig*>& depth, std::map<std::string, std::vector<SVCluster*>>& vars, std::set <std::string>& unmapped, std::map <std::string, SVtig*>& final_svtigs);
		void generate_fasta_file(parameters& params, faidx_t*& fasta_index, std::set <std::string>& reads, std::string file_path);
		int write_svtigs(std::string& f_path, const std::string& f_name, int pos, std::string& contig, int coverage, std::ostream& fp_write);
		int merge_svtigs(parameters& params, const std::string& dir);
		int final_assembly(parameters& params, faidx_t*& fasta_index, std::set <std::string>& read_set, std::string& svtig_name, double& contig_depth, SVCluster*& sv, std::map <std::string, SVtig*>& final_svtigs, int threads = -1);

};

#endif
