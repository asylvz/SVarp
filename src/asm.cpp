#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <chrono>
#include "phasing.h"
#include "asm.h"
#include "common.h"
#include "alignment.h"


int primary_cnt_asm = 0, secondary_cnt_asm = 0, insertion_cnt_asm = 0, deletion_cnt_asm = 0;

inline bool cmp_asm(Variant* i1, Variant* i2)
{
	if (i1->contig != i2->contig)
		return (i1->contig < i2->contig);
		

	if (i1->pos_in_ref != i2->pos_in_ref)
    	return (i1->pos_in_ref < i2->pos_in_ref);
	else
		return (i1->pos_in_ref + i1->sv_size > i2->pos_in_ref + i2->sv_size);
}


int find_var_asm(std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::vector<Variant*>& vars, Gaf& line, std::set <std::string>& unmapped)
{

	std::vector<int> cigarLen;
	std::vector<char> cigarOp;
	
	std::string rname = line.query_name.substr(0, line.query_name.find(' '));
	std::string cigar;
	
	//Check if the read is unmapped
	if (line.query_start == 0 && line.query_end == 0)
		unmapped.insert(rname);
	
	if(line.mapping_quality < MINMAPQ)
		return RETURN_ERROR;
	
	
	if (line.is_primary == false)
	{
		secondary_cnt_asm++;
		return -1;
	}
	else
		primary_cnt_asm++;
	
	cigarLen.clear();
	cigarOp.clear();

	int cigar_cnt = decompose_cigars(line.cigar, cigarLen, cigarOp);
	int base_pos = 0;
	for (int c = 0; c < cigar_cnt; c++)
	{
		if (cigarOp[c] == INSERTION && cigarLen[c] > MINSVSIZE)
		{
			Variant* var = generate_sv_node(gfa, line, base_pos + 1, cigarLen[c], INSERTION);
					
			if (var)
			{
				var->sv_size = cigarLen[c];
				var->sv_type = INSERTION;
				var->reads_untagged.insert(rname);
				vars.push_back(var);
				insertion_cnt_asm++;
				//std::cout<<"INSERTION " <<var->contig<<"\t"<<var->pos_in_ref<<"\t"<<var->sv_size<< "\t"<<var_name<<"\n";
			}
			else
				std::cout<<"RETURNED NULL\n";
		}
		else if (cigarOp[c] == DELETION && cigarLen[c] > MINSVSIZE)
		{
			Variant* var = generate_sv_node(gfa, line, base_pos + 1, cigarLen[c], DELETION);
			if (var)
			{
				var->sv_type = DELETION;
				var->sv_size = cigarLen[c];	
				var->reads_untagged.insert(rname);
				vars.push_back(var);
				deletion_cnt_asm++;
				//std::cout<<"DELETION " <<var->contig<<"\t"<<var->pos_in_ref<<"\t"<<var->sv_size<< "\t"<<var_name<<"\n";
			}
			else
				std::cout<<"RETURNED NULL\n";
		}

		if (cigarOp[c] != 'I')
			base_pos += cigarLen[c];
	}
	//Check the read depth/coverage here
	contig_coverage(ref, gfa, line);

	return RETURN_SUCCESS;
}

int output_vcf(parameters& params, std::vector<Variant*>& svs, int &dup_cnt, int &ins_cnt, int &del_cnt)
{

	std::map <std::string, phase*> phased_reads;
	std::ofstream fp_write(params.vcf_path);
	
	//Sort the svs by contig, position
	std::sort(svs.begin(), svs.end(), cmp_asm);

	
	//Phase variants	
	for (auto &a: svs)
		a->genotype = "./.";

	/*if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
	{	
		for (auto &a: svs)
		{
			std::string read = *(a->reads_untagged).begin();
			if (phased_reads.find(read) == phased_reads.end() || phased_reads[read]->haplotype == "none" || phased_reads[read]->phase_set == "none")
				continue;				

			if (phased_reads[read]->haplotype == "H1")
				a->genotype = "1/0";
			else if (phased_reads[read]->haplotype == "H2")
				a->genotype = "0/1";	
		}
	}*/

	//Find and remove overlaps
	for(unsigned int i = 0; i < svs.size(); i++)
	{
		Variant *v = svs[i];
		if(v->duplicate)
			continue;

		for(unsigned int j = 0; j < svs.size(); j++)
		{
			if(svs[i]->contig != svs[j]->contig)
				continue;
			if(svs[j]->duplicate || i == j)
				continue;
			
			double overlap = overlap_ratio(svs[i]->pos_in_ref, svs[i]->pos_in_ref + svs[i]->sv_size, svs[j]->pos_in_ref, svs[j]->pos_in_ref + svs[j]->sv_size);
			
			if (overlap > 0.9)
			{
				//if((tmp_svtig[i]->end - tmp_svtig[i]->start) >= (tmp_svtig[j]->end - tmp_svtig[j]->start))
				/*if (svs[i]->genotype != "./." && svs[j]->genotype != "./." && svs[i]->genotype != svs[j]->genotype)
				{
					svs[i]->genotype = "1/1";
					svs[j]->genotype = "1/1";
				}	
				*/

				if(svs[i]->sv_size >= svs[j]->sv_size)
				{
					svs[j]->duplicate = true;
					continue;	
				}
				else
				{
					v->duplicate = true;
					break;
				}
			}
		}
	}

	//Print the non-overlapping SVs
	for (auto &a: svs)
	{
		if (a->duplicate == false)
		{
			if(a->sv_type == DELETION)
			{
				fp_write<<a->contig<<"\t"<<a->pos_in_ref<<"\t"<<a->pos_in_ref + a->sv_size<<"\t"<<a->genotype<<"\tDEL\n";
				del_cnt++;
			}
			else
			{
				fp_write<<a->contig<<"\t"<<a->pos_in_ref<<"\t"<<a->pos_in_ref + a->sv_size<<"\t"<<a->genotype<<"\tINS\n";
				ins_cnt++;
			}
		}
		else
			dup_cnt++;
	}

	fp_write.close();

	return RETURN_SUCCESS;
}


int read_alignments_asm(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::set <std::string>& unmapped)
{
	int line_count = 0;
	std::string line;
	std::vector <Variant*> svs;

	std::cout<<"Reading the GAF file"<<std::endl;
	auto t1 = std::chrono::steady_clock::now();
	
	std::ifstream fp(params.gaf);
	if(!fp.good())
	{
		std::cerr << "Error opening '"<<params.gaf<< std::endl;
        return RETURN_ERROR;
	}

	while(getline(fp, line))
	{
		line_count++;

		Gaf g = parse_gaf_line(line);		
		find_var_asm(ref, gfa, svs, g, unmapped);

		if(line_count > TEST_SAMPLE_SIZE)
			break;
	}
			
	auto t2 = std::chrono::steady_clock::now();	
	int dup_cnt = 0, ins_cnt = 0, del_cnt = 0; 	
	output_vcf(params, svs, dup_cnt, ins_cnt, del_cnt);

	std::cout<<"--->execution time: "<<std::chrono::duration<double> (t2 - t1).count()<<"sec.\n";
	std::cout<<"--->"<<primary_cnt_asm<<" primary mappings and "<<insertion_cnt_asm<<" insertion, "<<deletion_cnt_asm<<" deletion loci in the cigar\n";
	std::cout<<"--->there are "<<unmapped.size()<<" unmapped alignments\n";
	std::cout<<"--->found "<<ins_cnt<<" insertions and "<<del_cnt<<" deletions...(filtered "<<dup_cnt<<" duplicate SVs)\n";

	std::map<std::string, Contig*>::iterator it2;
	for (it2=ref.begin(); it2 != ref.end(); ++it2)
	{
		it2->second->coverage = (double) it2->second->mapped_bases / it2->second->contig_length;	
		params.fp_logs << it2->first<<"---> LEN= "<<it2->second->contig_length<<" - Mapped bases= "<< it2->second->mapped_bases<<" - Mapped reads= "<<it2->second->mapped_reads << " - Cov= "<<it2->second->coverage <<std::endl;		
	}			
	
	return RETURN_SUCCESS;
}
