#include <iostream>
#include <stdlib.h>
#include <errno.h>
#include <thread>
#include <chrono>
#include <sys/wait.h>
#include <signal.h>
#include <spawn.h>
#include <string.h>
#include <algorithm>
#include <sstream>
#include "common.h"
#include <filesystem>
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif


extern char **environ;

std::mutex g_log_mtx;

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
		if (tok.rfind("tp:A:", 0) == 0)
		{
			if (tok.substr(5, 6) != "P")
				gafline.is_primary = false;
		}
		else if (tok.rfind("AS:f:", 0) == 0)
		{
			try { gafline.aln_score = std::stof(tok.substr(5)); }
			catch (const std::exception&) { return RETURN_ERROR; }
		}
		else if (tok.rfind("cg:Z:", 0) == 0)
			gafline.cigar = tok.substr(5);
		else if (tok.rfind("id:f:", 0) == 0)
		{
			try { gafline.identity = std::stof(tok.substr(5)); }
			catch (const std::exception&) { return RETURN_ERROR; }
		}
	}

	return RETURN_SUCCESS;
}


// Largest of three ratios, so an interval contained in another scores 1. A
// reciprocal criterion (the smaller of the two one-sided ratios) is the usual
// choice when comparing SV calls and is the alternative worth revisiting.
double overlap_ratio(int x_start, int x_end, int y_start, int y_end)
{
	int overlap = std::max(0, std::min(x_end, y_end) - std::max(x_start, y_start));
	int total_length = x_end - x_start + y_end - y_start;
	int x_length = x_end - x_start;
	int y_length = y_end - y_start;

	if (x_length <= 0 || y_length <= 0 || total_length <= 0)
		return 0;

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
}


//First line of `bin flag`; tools that answer with their usage (old samtools) are run bare and
//their "Version" line is taken. "unknown" when nothing useful is printed.
std::string tool_version(const std::string& bin, const std::string& flag)
{
	if (bin.empty())
		return "unknown";
	auto first_line = [](const std::string& out) {
		std::string::size_type b = out.find_first_not_of(" \t\r\n");
		if (b == std::string::npos)
			return std::string();
		std::string::size_type e = out.find('\n', b);
		return out.substr(b, e == std::string::npos ? std::string::npos : e - b);
	};
	std::string line = first_line(exec(bin + " " + flag + " 2>&1", true));
	bool usage = line.rfind("Usage", 0) == 0 || line.rfind("Program", 0) == 0 || line.rfind("[main]", 0) == 0 || line.find("unrecognized") != std::string::npos || line.find("invalid option") != std::string::npos;
	if (usage)
	{
		std::istringstream in(exec(bin + " 2>&1", true));
		std::string l; line.clear();
		while (std::getline(in, l))
			if (l.find("ersion") != std::string::npos) { line = first_line(l); break; }
	}
	if (line.empty() || line.find("not found") != std::string::npos || line.find("No such file") != std::string::npos)
		return "unknown";
	if (line.size() > 100)
		line.resize(100);
	return line;
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
	pclose(pipe);
	return "Success";
}


//Run cmd under /bin/sh in its own process group; after timeout_seconds (0 = none) the whole group is
//killed and the status reads as exit 124, like coreutils timeout. Returns a wait status, -1 on failure.
//posix_spawn, as in glibc's system(), avoids duplicating the page tables of a large process.
static int run_shell(const std::string& cmd, int timeout_seconds)
{
    posix_spawnattr_t attr;
    posix_spawnattr_init(&attr);
    posix_spawnattr_setflags(&attr, POSIX_SPAWN_SETPGROUP);
    posix_spawnattr_setpgroup(&attr, 0);
    const char* argv[] = {"sh", "-c", cmd.c_str(), nullptr};
    pid_t pid = -1;
    int err = posix_spawn(&pid, "/bin/sh", nullptr, &attr, const_cast<char* const*>(argv), environ);
    posix_spawnattr_destroy(&attr);
    if (err != 0) {
        errno = err;
        return -1;
    }

    int status = 0;
    if (timeout_seconds <= 0)
        return waitpid(pid, &status, 0) == pid ? status : -1;

    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(timeout_seconds);
    while (true) {
        pid_t r = waitpid(pid, &status, WNOHANG);
        if (r == pid) return status;
        if (r < 0) return -1;
        if (std::chrono::steady_clock::now() >= deadline) {
            kill(-pid, SIGKILL);
            waitpid(pid, &status, 0);
            return 124 << 8;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
}

int run_and_log(const std::string& cmd, parameters& params,
                const std::string& label, int retries,
                int backoff_seconds, bool fatal, int timeout_seconds)
{
	if (std::lock_guard<std::mutex> lk(g_log_mtx); params.debug && params.fp_logs.is_open()) {
        params.fp_logs << "[run_and_log] " << label << " CMD: " << cmd << "\n";
    }
    int attempt = 0;

    while (true) {
        int rc = run_shell(cmd, timeout_seconds);
        if (rc == 0) return 0;

        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open()) {
            if (rc == -1) {
                params.fp_logs << "Failed to run '" << label << "' (" << cmd
                               << ") fork/wait error: " << strerror(errno) << "\n";
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
            if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
                params.fp_logs << "Retrying in " << sleep_seconds
                               << " seconds... (attempt "
                               << (attempt + 1) << ")\n";
            std::this_thread::sleep_for(std::chrono::seconds(sleep_seconds));
            attempt++;
            continue;
        }

        if (fatal) {
            if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
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

    // Resolve extra_dirs relative to the svarp binary's directory
    fs::path exe_dir;
    {
        char buf[4096];
#ifdef __linux__
        ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
#elif defined(__APPLE__)
        uint32_t sz = sizeof(buf);
        int ret = _NSGetExecutablePath(buf, &sz);
        ssize_t len = (ret == 0) ? strlen(buf) : -1;
#else
        ssize_t len = -1;
#endif
        if (len > 0) {
            buf[len] = '\0';
            exe_dir = fs::path(buf).parent_path().parent_path(); // build/../ = SVarp root
        }
    }

    for (const auto& d : extra_dirs) {
        fs::path p = fs::path(d) / progname;
        if (fs::exists(p) && access(p.c_str(), X_OK) == 0)
            return p.string();
        // Try relative to svarp binary's root directory
        if (!exe_dir.empty()) {
            fs::path p2 = exe_dir / d / progname;
            if (fs::exists(p2) && access(p2.c_str(), X_OK) == 0)
                return p2.string();
        }
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
			case 'a':
				seq[i] = 't';
				break;
			case 'c':
				seq[i] = 'g';
				break;
			case 'g':
				seq[i] = 'c';
				break;
			case 't':
				seq[i] = 'a';
				break;
		}
	}
	return seq;
}

std::string current_timestamp()
{
	auto now = std::chrono::system_clock::now();
	std::time_t t = std::chrono::system_clock::to_time_t(now);
	char buf[64];
	std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
	return std::string(buf);
}

std::string format_duration(double seconds)
{
	int total = static_cast<int>(seconds);
	int h = total / 3600;
	int m = (total % 3600) / 60;
	int s = total % 60;
	char buf[32];
	if (h > 0)
		std::snprintf(buf, sizeof(buf), "%dh %dm %ds", h, m, s);
	else if (m > 0)
		std::snprintf(buf, sizeof(buf), "%dm %ds", m, s);
	else
		std::snprintf(buf, sizeof(buf), "%.1fs", seconds);
	return std::string(buf);
}

void log_step(LogFile& fp_logs, const std::string& step)
{
	if (fp_logs.is_open())
		fp_logs << "\n[" << current_timestamp() << "] " << step << "\n";
}
