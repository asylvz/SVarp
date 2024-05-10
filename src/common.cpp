#include <iostream>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <algorithm>
#include <sstream>
#include "common.h"


//Returns a vector of incoming or outgoing nodes based on the input map object
/*const std::vector<std::string>& find_prev_next_nodes(std::map <std::string, std::vector<std::string>> inout_nodes, std::string node)
{
	std::map<std::string, std::vector<std::string>>::iterator it;
	for (it=inout_nodes.begin(); it != inout_nodes.end(); ++it)
	{
		if(it->first == node)
			return it->second;
	}
}
*/

int parse_gaf_line(std::string& line, Gaf& gafline)
{
	//Gaf gafline;	
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
	gafline.mapping_quality = stoi(tokens[11]);
	gafline.residue_matches = stoi(tokens[9]);
	gafline.alignment_block_length = stoi(tokens[10]);
	gafline.is_primary = true;

	for (auto& tok : tokens) 
	{
		if(strstr(tok.c_str(), "tp:A:"))
		{
			if (tok.substr(5, 6) != "P")
				gafline.is_primary = false;
		}
		if(strstr(tok.c_str(), "AS:f:"))
			gafline.aln_score = std::stof(tok.substr(5));
		
	
		if(strstr(tok.c_str(), "cg:Z:"))
			gafline.cigar = tok.substr(5);
	}

	return RETURN_SUCCESS;
}


double overlap_ratio(int x_start, int x_end, int y_start, int y_end)
{
	int overlap = std::max(0, std::min(x_end, y_end) - std::max(x_start, y_start));
	int total_length = x_end - x_start + y_end - y_start;
	int x_length = x_end - x_start;
	int y_length = y_end - y_start;
	
	double a = (double) 2 * (overlap / (double) total_length);
	double b = (double) overlap / (double) x_length;
	double c = (double) (overlap / (double) y_length);
	if (a > b)
	{
		if (c > a)
			return c;
		else
			return a;
	}
	else
	{
		if (c > b)
			return c;
		else
			return b;
	}
	//std::cout<<a<<" "<<b<<" "<<c<<" "<<max_overlap<<"\n";

	return -1;
}


void error(const char* const msg)
{
	std::cerr<<msg<<std::endl;
    exit(EXIT_FAILURE);
}

std::string exec(std::string command, bool return_out) 
{
	FILE* pipe = popen(command.c_str(), "r");
	if (!pipe)
		return "Error";

	if (return_out)
   	{
		char buffer[128];
   	   	std::string result = "";
	   	while (!feof(pipe)) 
			if (fgets(buffer, 128, pipe) != NULL)
				result += buffer;
   		
		pclose(pipe);
   		return result;
	}
	return "Success";
}


/* The codes for the dictionary is taken from The C Programming language 
 * 2nd edition - Brian Kernighan and Dennis Ritchie */
unsigned hash(char *s)
{
    unsigned hashval;
    for (hashval = 0; *s != '\0'; s++)
      hashval = *s + 31 * hashval;
    return hashval % HASHSIZE;
}

void* getMem(size_t size)
{
	void* ret;

	ret = malloc(size);
	if(ret == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(0);
	}

	return ret;
}

/*void init_params(parameters** params)
{
	// initialize parameters
	*params = (parameters*) getMem(sizeof(parameters));
	(*params)->ref_graph = NULL;
	(*params)->gaf = NULL;
}*/

void set_str(char** target, char* source)
{
	if(*target != NULL)
	{
		free((*target));
	}

	if(source != NULL)
	{
		(*target) = (char*) getMem(sizeof(char) * (strlen(source) + 1));
		strncpy( (*target), source, (strlen(source) + 1));
	}
	else
	{
		(*target) = NULL;
	}
}

void print_error( char* msg)
{
	/* print error message and exit */
	printf("\n%s\n", msg);
	printf("Invoke parameter -h for help.\n");
	exit(EXIT_COMMON);
}


int decompose_cigars(std::string cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp)
{
	size_t cigar_offset = 0, str_offset = 0, cigar_cnt = 0;
	char* cigar_copy = (char*) cigar.c_str();	
	//std::cout<<cigar_copy<<"\n";

	while(cigar_offset < cigar.length())
	{
		if (isdigit(*(cigar_copy + cigar_offset)) == 0)
		{	
			std::string s = "";
			for (int z = str_offset; z > 0; z--)
				s += *(cigar_copy + cigar_offset - z);
			
			cigarOp.push_back (*(cigar_copy + cigar_offset));
			cigarLen.push_back (stoi(s));
			
			//std::cout<<*(cigar_copy + cigar_offset) << " "<< stoi(s) <<"\n";
			str_offset = 0;
			cigar_cnt++;			
		}
		else 
			str_offset++;

		cigar_offset++;
		
	}
	return cigar_cnt;
}


std::string& reverse_complement(std::string& seq)
{
	std::reverse(seq.begin(), seq.end());
	for (std::size_t i = 0; i < seq.length(); ++i)
	{
		switch (seq[i])
		{
			case 'A':
				seq[i] = 'T';
				break;
			case 'C':
				seq[i] = 'G';
				break;
			case 'G':
				seq[i] = 'C';
				break;
			case 'T':
				seq[i] = 'A';
				break;
		}
	}
	return seq;
}

/* void calculate_n50_phase(parameters* params, std::map <std::string, phase*> phased_reads)
{
	std::cout<<"\nCalculate N50\n";
	std::map<std::string, unsigned long> fasta_index;
	index_fasta(params, fasta_index);	
	std::map<std::string, phase*>::iterator itr;
	
	std::vector <int> phased_sizes;
	std::vector <int> unphased_sizes;

	std::ifstream fp_read(params->fasta);
	
	long unphased_total = 0, phased_total = 0, all_total = 0, total_line = 0, none_line = 0;
	for (itr=phased_reads.begin(); itr != phased_reads.end(); ++itr)
	{	
			
		long read_size = 0;	
		long char_pos;	
		std::string line;	
		std::string read = itr->first;
		

		if (fasta_index.find(read)!=fasta_index.end())
		{
			char_pos = fasta_index[read];

			fp_read.seekg(char_pos, std::ios::beg);
			
			//std::cout<<read<<" "<<char_pos<<std::endl;	
			getline(fp_read, line);
			//cout<<read<<" - "<<line<<std::endl;
			line.clear();
			while(getline(fp_read, line))
			{
				//std::cout<<line<<std::endl;
				if(line[0] == '>')
					break;
			
				//std::cout<<line.length()<<std::endl;
				read_size += line.length();
				line.clear();
				//std::cout<<read_size<<std::endl;
			}
			//std::cout<<read<<" "<<read_size<<std::endl;
		}
		else
		{
			std::cout<<"Not found "<<read<<std::endl;	
			continue;
		}
	

		if(itr->second->haplotype == "none" || itr->second->phase_set == "none")
		{
			none_line++;
			unphased_total += read_size;
			unphased_sizes.push_back(read_size);
			all_total += read_size;
		}
		else
		{
			total_line++;
			phased_total += read_size;
			phased_sizes.push_back(read_size);
			all_total += read_size;
		}

		//std::cout<<phased_total <<" "<<unphased_total<<" "<<read_size<<" " <<all_total <<std::endl;	
		//cout<<none_line<<" - "<<total_line<<std::endl;
	}

	std::cout<<"Unphased = "<<(double) unphased_total / none_line <<"\nPhased = "<<phased_total/total_line<<std::endl;

	sort(phased_sizes.begin(), phased_sizes.end(), std::greater<>());
	sort(unphased_sizes.begin(), unphased_sizes.end(), std::greater<>());
	
	long tmp_total = 0;
	for(size_t i = 0; i< phased_sizes.size(); i++)
	{
		//std::cout<<"PHASED "<<phased_sizes[i]<<std::endl;
		tmp_total += phased_sizes[i];

		if (tmp_total >= (phased_total / 2))
		{
			std::cout<<"(phased total) N50 of phased is "<<phased_sizes[i]<<std::endl;
			break;
		}
	}

	tmp_total = 0;
	for(size_t i = 0; i< unphased_sizes.size(); i++)
	{
		//std::cout<<"UNPHASED "<<unphased_sizes[i]<<std::endl;
		tmp_total += unphased_sizes[i];

		if (tmp_total >= (unphased_total / 2))
		{
			std::cout<<"(unphased total) N50 of unphased is "<<unphased_sizes[i]<<std::endl;
			break;
		}
	}			
}*/
