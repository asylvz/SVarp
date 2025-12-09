#include <iostream>
#include <stdlib.h>
#include <errno.h>
#include <thread>
#include <chrono>
#include <sys/wait.h>
#include <stdarg.h>
#include <string.h>
#include <algorithm>
#include <sstream>
#include "common.h"
#include <filesystem>


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
    gafline.aln_score = 0.0f;

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

    //std::cout<<gafline.aln_score<<" - " << gafline.mapping_quality<<"\n";
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

std::string exec(const std::string& command, bool return_out) 
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


int run_and_log(const std::string& cmd, parameters& params,
                const std::string& label, int retries,
                int backoff_seconds, bool fatal)
{
	if (params.debug && params.fp_logs.is_open()) {
        params.fp_logs << "[run_and_log] " << label << " CMD: " << cmd << "\n";
    }
    int attempt = 0;

    while (true) {
        int rc = system(cmd.c_str());
        if (rc == 0) return 0;

        // Hatalar sadece log dosyasına yazılır, ekrana asla basılmaz
        if (params.fp_logs.is_open()) {
            if (rc == -1) {
                params.fp_logs << "Failed to run '" << label << "' (" << cmd
                               << ") system() error: " << strerror(errno) << "\n";
            }
#ifdef __unix__
            else if (WIFEXITED(rc)) {
                params.fp_logs << "Command '" << label << "' (" << cmd
                               << ") exited with status " << WEXITSTATUS(rc) << "\n";
            } else if (WIFSIGNALED(rc)) {
                params.fp_logs << "Command '" << label << "' (" << cmd
                               << ") terminated by signal " << WTERMSIG(rc) << "\n";
            } else {
                params.fp_logs << "Command '" << label << "' (" << cmd
                               << ") returned " << rc << "\n";
            }
#else
            else {
                params.fp_logs << "Command '" << label << "' returned " << rc << "\n";
            }
#endif
        }

        if (attempt < retries) {
            int sleep_seconds = backoff_seconds * (1 << attempt);
            if (params.fp_logs.is_open())
                params.fp_logs << "Retrying in " << sleep_seconds
                               << " seconds... (attempt "
                               << (attempt + 1) << ")\n";
            std::this_thread::sleep_for(std::chrono::seconds(sleep_seconds));
            attempt++;
            continue;
        }

        if (fatal) {
            if (params.fp_logs.is_open())
                params.fp_logs << "Fatal: command '" << cmd
                               << "' failed after " << (attempt + 1)
                               << " attempts\n";
            exit(EXIT_COMMON);
        }

        return rc;
    }
}

// Helper: find an executable in extra dirs, current dir and PATH
std::string find_executable(const std::string &progname,
                            const std::vector<std::string> &extra_dirs)
{
    namespace fs = std::filesystem;

    // 1) Absolute path ise ve çalıştırılabilir ise direkt kullan
    if (!progname.empty() && progname[0] == '/') {
        if (fs::exists(progname) && access(progname.c_str(), X_OK) == 0)
            return progname;
    }

    // 2) Ek klasörler
    for (const auto& d : extra_dirs) {
        fs::path p = fs::path(d) / progname;
        if (fs::exists(p) && access(p.c_str(), X_OK) == 0)
            return p.string();
    }

    // 3) Çalışılan dizin ./progname
    {
        fs::path p = fs::current_path() / progname;
        if (fs::exists(p) && access(p.c_str(), X_OK) == 0)
            return p.string();
    }

    // 4) PATH içinde ara
    const char* path_env = std::getenv("PATH");
    if (path_env) {
        std::string path(path_env);
        std::string::size_type start = 0;
        while (true) {
            auto pos = path.find(':', start);
            std::string dir = (pos == std::string::npos)
                                  ? path.substr(start)
                                  : path.substr(start, pos - start);
            if (!dir.empty()) {
                fs::path p = fs::path(dir) / progname;
                if (fs::exists(p) && access(p.c_str(), X_OK) == 0)
                    return p.string();
            }
            if (pos == std::string::npos) break;
            start = pos + 1;
        }
    }

    return ""; // bulunamadı
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


int decompose_cigars(const std::string& cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp)
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
