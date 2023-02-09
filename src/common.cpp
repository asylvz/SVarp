#include <iostream>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "common.h"


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
	/*get the Cigar*/
	size_t cigar_offset = 0, str_offset = 0, cigar_cnt = 0;
	char* cigar_copy = (char*) cigar.c_str();	
	char *tmp_str = new char[6];
	while(cigar_offset < cigar.length())
	{
		if (isdigit(*(cigar_copy + cigar_offset)) == 0)
		{
			cigarOp.push_back (*(cigar_copy + cigar_offset));
			cigarLen.push_back (atoi(tmp_str));
			
			delete[] tmp_str; 	
			tmp_str = new char[6];
			str_offset = 0;
			cigar_cnt++;			
		}
		else 
		{
			*(tmp_str + str_offset) = *(cigar_copy + cigar_offset);
			str_offset++;
		}
		cigar_offset++;
	}
	delete[] tmp_str;
	
	return cigar_cnt;
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
