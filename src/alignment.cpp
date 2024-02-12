#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <zlib.h>
#include <chrono>
#include "alignment.h"

int primary_cnt = 0, secondary_cnt = 0, inter_cnt = 0, intra_cnt = 0, insertion_cnt = 0, deletion_cnt = 0;



Gaf parse_gaf_line(std::string& line)
{
	Gaf gafline;	
	std::vector <std::string> tokens;	
	
	std::string tmp_str;
	std::stringstream s(line);
	while(getline(s, tmp_str, '\t'))
		tokens.push_back(tmp_str);
	
	gafline.query_name = tokens[0].substr(0, tokens[0].find(' '));
	gafline.query_length = stoi(tokens[1]);
	gafline.query_start = stoi(tokens[2]);
	gafline.query_end = stoi(tokens[3]);
	gafline.strand = tokens[4];
	gafline.path = tokens[5];
	gafline.path_length = stoi(tokens[6]);
	gafline.path_start = stoi(tokens[7]);
	gafline.path_end = stoi(tokens[8]);
	gafline.residue_matches = stoi(tokens[9]);
	gafline.alignment_block_length = stoi(tokens[10]);
	gafline.mapping_quality = stoi(tokens[11]);
	gafline.is_primary = true;


	for (auto& tok : tokens) 
	{
		if(strstr(tok.c_str(), "tp:A:"))
		{
			if (tok.substr(5, 6) != "P")
				gafline.is_primary = false;
		}
	
		if(strstr(tok.c_str(), "cg:Z:"))
			gafline.cigar = tok.substr(5);
	}

	return gafline;
}


int find_var(std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, Gaf& line, std::map <std::string, int>& read_freq, std::set <std::string>& unmapped)
{
	std::vector<int> cigarLen;
	std::vector<char> cigarOp;
	
	//Check if the read is unmapped
	if ((line.query_start == 0) && (line.query_end == 0))
		unmapped.insert(line.query_name);
	
	if(line.mapping_quality < MINMAPQ)
		return RETURN_ERROR;
	
	//If the read has multiple mappings, then the ends are putative SV loci
	std::map<std::string, int>::iterator it = read_freq.find(line.query_name);
	if (it != read_freq.end())
		inter_cnt += mapping_start_end(gfa, line, vars);
	
	if (line.is_primary == false)
	{
		secondary_cnt++;
		return -1;
	}
	else
		primary_cnt++;
	
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
				//std::cout<<var->contig<<" - "<< var->pos_in_ref<<"\n";
				var->sv_size = cigarLen[c];	
				intra_cnt++;
				std::string var_name = var->node + ":" + std::to_string(var->pos_in_node);

				std::map<std::string, Variant*>::iterator it = vars.find(var_name);
				if (it != vars.end())
					it->second->reads_untagged.insert(line.query_name);
				else
				{
					var->reads_untagged.insert(line.query_name);
					vars.insert(std::pair<std::string, Variant*>(var_name, var));
					insertion_cnt++;
				}
			}
			else
				std::cout<<"RETURNED NULL\n";
		}
		else if (cigarOp[c] == DELETION && cigarLen[c] > MINSVSIZE)
		{
			Variant* var = generate_sv_node(gfa, line, base_pos + 1, cigarLen[c], DELETION);
						
			if (var)
			{
				//std::cout<<var->contig<<" - "<< var->pos_in_ref<<"\n";
				var->sv_size = cigarLen[c];	
				intra_cnt++;
						
				std::string var_name = var->node + ":" + std::to_string(var->pos_in_node);
				std::map<std::string, Variant*>::iterator it = vars.find(var_name);
						
				if (it != vars.end())
					it->second->reads_untagged.insert(line.query_name);
				else
				{
					var->reads_untagged.insert(line.query_name);
					vars.insert(std::pair<std::string, Variant*>(var_name, var));
					deletion_cnt++;
				}
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

int read_gz(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::set <std::string>& unmapped, std::map <std::string, int>& read_freq)
{

	std::map<std::string, int>::iterator it;	
	std::map <std::string, int> read_freq_tmp;
	std::string line;
	
	gzFile myfile = gzopen((params.gaf).c_str(), "rb");		
	if (!myfile)
		error("Error: Error opening .gz file");

	char buffer[BUFLEN];	
	char *offset = buffer;
	
	while(1)
	{
		int err, len = sizeof(buffer) - (offset - buffer);
		if (len == 0) 
			error("Error: Buffer too small");

		len	= gzread(myfile, offset, len);
		if (len == 0) 
			break;   
    	if (len < 0) 
			error(gzerror(myfile, &err));
		
		char* cur = buffer;
    	char* end = offset + len;	
	
 	   	for (char* eol; (cur < end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1)
    	{
			line = std::string(cur, eol);
			Gaf g = parse_gaf_line(line);

			
			it = read_freq_tmp.find(g.query_name);
			if (it != read_freq_tmp.end())
				it->second++;
			else
				read_freq_tmp.insert(std::pair<std::string, int>(g.query_name, 1));

		}
    	offset = std::copy(cur, end, buffer);
	}

	int cnt_multiple = 0, cnt_single = 0;	
	for (it=read_freq_tmp.begin(); it != read_freq_tmp.end(); ++it)
	{
		if (it->second > 1)
		{
			cnt_multiple++;
			read_freq.insert(std::pair<std::string, int>(it->first, it->second));
		}
		else
			cnt_single++;
	}
	//std::cout<<"Single: "<<cnt_single<<" Multiple: "<<cnt_multiple<<"\n";
	int line_count = 0;
	gzrewind(myfile);
	
	memset(buffer, '\0', sizeof(buffer));
	offset = buffer;
	bool tst = false;
	while(1)
	{
		int err, len = sizeof(buffer) - (offset - buffer);
		if (len == 0) 
		{
			std::cout<<"Line = "<< line_count<<"\n";
			error("Error: Buffer too small");
		}

		len	= gzread(myfile, offset, len);

		if (len == 0) break;   
    	if (len < 0) error(gzerror(myfile, &err));
		
		char* cur = buffer;
    	char* end = offset + len;	
	
 	   	for (char* eol; (cur < end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1)
    	{
			line = std::string(cur, eol);
			line_count++;
			
			Gaf g = parse_gaf_line(line);		
			find_var(ref, gfa, vars, g, read_freq, unmapped);
			
			if(line_count > TEST_SAMPLE_SIZE)
			{
				tst = true;
				break;
			}
		}
		if (tst)
			break;
			
    	offset = std::copy(cur, end, buffer);
	}

   	if (gzclose(myfile) != Z_OK) 
		error("ERROR: GZCLOSE FAILED");
	
	return RETURN_SUCCESS;
}


int read_alignments(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::set <std::string>& unmapped)
{

	std::cout<<"Reading the GAF file"<<std::endl;
	auto t1 = std::chrono::steady_clock::now();
		
	std::map <std::string, int> read_freq_tmp, read_freq;
	std::map<std::string, int>::iterator it;
	std::filesystem::path gaf_path = params.gaf;
	
	if (gaf_path.extension() == ".gz")
		read_gz(params, ref, gfa, vars, unmapped, read_freq);
	else
	{
		std::string line;
		std::ifstream fp(params.gaf);
		
		if(!fp.good())
		{
			std::cerr << "Error opening '"<<params.gaf<< std::endl;
        	return RETURN_ERROR;
		}

		while(getline(fp, line))
		{
			Gaf g = parse_gaf_line(line);
				
			it = read_freq_tmp.find(g.query_name);
			if (it != read_freq_tmp.end())
				it->second++;
			else
				read_freq_tmp.insert(std::pair<std::string, int>(g.query_name, 1));
		}

		int cnt_multiple = 0, cnt_single = 0;	
		for (it=read_freq_tmp.begin(); it != read_freq_tmp.end(); ++it)
		{
			if (it->second > 1)
			{
				cnt_multiple++;
				read_freq.insert(std::pair<std::string, int>(it->first, it->second));
			}
			else
				cnt_single++;
		}
		//std::cout<<"Single: "<<cnt_single<<" Multiple: "<<cnt_multiple<<"\n";
		int line_count = 0;


		fp.clear() ; // clear the failed state of the stream
		fp.seekg(0) ; // seek to the first character in the file
				
		while(getline(fp, line))
		{
			Gaf g = parse_gaf_line(line);
			line_count++;
			find_var(ref, gfa, vars, g, read_freq, unmapped);
			
			if(line_count > TEST_SAMPLE_SIZE)
				break;
		}
	}	
	
	/*for (auto &a: vars)
	{
		//if (a.second->reads_h1.size() > 1 || a.second->reads_h1.size() == 0)
		//	std::cout<<a.first<<"\n"<< a.second->reads_h1.size()<<"\n\n";
		if (a.second->type == INTRA && a.second->sv_type == DELETION)
			std::cout<<a.second->contig<<"\t"<<a.second->pos_in_ref<<"\t"<<a.second->sv_size<<"\t"<<a.second->reads_untagged.size()<<"\n";
	}*/
	auto t2 = std::chrono::steady_clock::now();

	std::cout<<"--->execution time: "<<std::chrono::duration<double> (t2 - t1).count()<<"sec.\n";
	std::cout<<"--->"<<primary_cnt<<" primary mappings and "<<insertion_cnt<<" insertion, "<<deletion_cnt<<" deletion loci in the cigar\n";
	std::cout<<"--->there are "<<inter_cnt + intra_cnt<<" SV signal (" <<inter_cnt<< " inter alignment and "<<intra_cnt<<" intra alignment)\n";
	std::cout<<"--->there are "<<unmapped.size()<<" unmapped alignments\n";
	
	params.fp_logs<<"--->"<<primary_cnt<<" primary mappings and "<<insertion_cnt<<" insertion, "<<deletion_cnt<<" deletion loci in the cigar\n";
	params.fp_logs<<"--->there are "<<inter_cnt + intra_cnt<<" SV signal (" <<inter_cnt<< " inter alignment and "<<intra_cnt<<" intra alignment)\n";
	params.fp_logs<<"--->there are "<<unmapped.size()<<" unmapped alignments\n";
	
	std::map<std::string, Contig*>::iterator it2;
	
	params.fp_logs<<"------->Contig coverage\n\n";
	for (it2=ref.begin(); it2 != ref.end(); ++it2)
	{
		it2->second->coverage = (double) it2->second->mapped_bases / it2->second->contig_length;	
		params.fp_logs<< it2->first<<"---> LEN= "<<it2->second->contig_length<<" - Mapped bases= "<< it2->second->mapped_bases<<" - Mapped reads= "<<it2->second->mapped_reads << " - Cov= "<<it2->second->coverage <<std::endl;		
	}			
	
	return RETURN_SUCCESS;
}

