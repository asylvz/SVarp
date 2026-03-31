#include <iostream>
#include <stdlib.h>
#include <errno.h>
#include <thread>
#include <chrono>
#include <sys/wait.h>
#include <string.h>
#include <algorithm>
#include <sstream>
#include "common.h"
#include <filesystem>
#include <unistd.h>


int parse_gaf_line(std::string& line, Gaf& gafline)
{
	//Gaf gafline;	
	std::vector <std::string> tokens;	
	
	std::string tmp_str;
	std::stringstream s(line);
	while(getline(s, tmp_str, '\t'))
		tokens.push_back(tmp_str);

	if (tokens.size() < 12)
		return RETURN_ERROR;

	try {
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
	} catch (const std::exception&) {
		return RETURN_ERROR;
	}
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
		{
			try { gafline.aln_score = std::stof(tok.substr(5)); }
			catch (const std::exception&) { return RETURN_ERROR; }
		}

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
			if (fgets(buffer, 128, pipe) != nullptr)
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

std::string find_executable(const std::string &progname,
                            const std::vector<std::string> &extra_dirs)
{
    namespace fs = std::filesystem;

    if (!progname.empty() && progname[0] == '/') {
        if (fs::exists(progname) && access(progname.c_str(), X_OK) == 0)
            return progname;
    }

    for (const auto& d : extra_dirs) {
        fs::path p = fs::path(d) / progname;
        if (fs::exists(p) && access(p.c_str(), X_OK) == 0)
            return p.string();
    }

    {
        fs::path p = fs::current_path() / progname;
        if (fs::exists(p) && access(p.c_str(), X_OK) == 0)
            return p.string();
    }

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

    return "";
}

int decompose_cigars(const std::string& cigar, std::vector<int>& cigarLen, std::vector<char>& cigarOp)
{
	size_t cigar_offset = 0, str_offset = 0, cigar_cnt = 0;
	const char* cigar_ptr = cigar.c_str();

	while(cigar_offset < cigar.length())
	{
		if (isdigit(*(cigar_ptr + cigar_offset)) == 0)
		{
			std::string s = "";
			for (int z = str_offset; z > 0; z--)
				s += *(cigar_ptr + cigar_offset - z);

			cigarOp.push_back (*(cigar_ptr + cigar_offset));
			try { cigarLen.push_back (stoi(s)); }
			catch (const std::exception&) { return -1; }
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
