#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <sys/wait.h>
#include <htslib/faidx.h>
// memory included indirectly
#include <zlib.h>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include "assembly.h"

// Removes a cluster's work directory on every exit path; in debug mode the
// read FASTA is kept under in/ first.
struct JobDir
{
	std::string dir, fasta, keep;
	~JobDir()
	{
		std::error_code ec;
		if (!keep.empty() && std::filesystem::exists(fasta, ec))
			std::filesystem::copy_file(fasta, keep, std::filesystem::copy_options::overwrite_existing, ec);
		std::filesystem::remove_all(dir, ec);
	}
};



//Bundled binaries first, then PATH
static std::string tool_bin(const std::string& name)
{
    static const std::vector<std::string> extra_dirs = {"third_party/wtdbg2", "dep/wtdbg2", "dep/minimap2", "dep/samtools"};
    std::string bin = find_executable(name, extra_dirs);
    return bin.empty() ? find_executable(name) : bin;
}

void Assembly::generate_fasta_file(parameters &params, faidx_t *&fasta_index, std::set<std::string> &reads, std::string file_path)
{
	std::ofstream fp_write(file_path);

	hts_pos_t loc_length;
	std::string line;
	const size_t line_len = 60;
	std::set<std::string> read_seqs;

	for (auto &read : reads)
	{
		hts_pos_t n = faidx_seq_len64(fasta_index, read.c_str());
		if (n <= 0)
			continue;
		char *tmp = faidx_fetch_seq64(fasta_index, read.c_str(), 0, n - 1, &loc_length);
		if (tmp == nullptr)
			continue;

		std::string seq(tmp);
		free(tmp);

		fp_write << ">" << read << std::endl;

		for (unsigned int i = 0; i < seq.size(); i += line_len)
		{
			std::string sub = (seq.substr(i, line_len));
			fp_write << sub << std::endl;
		}
		read_seqs.insert(seq);
	}

	fp_write.close();
}

int Assembly::write_svtigs(std::string &f_path, const std::string &f_name, int pos, std::string &contig, int coverage, std::ostream &fp_write)
{
	std::ifstream fp_read(f_path);
	std::string line;
	int contig_cnt = 0;
	std::string svtig_name = f_name;

	if (!fp_read)
	{
		std::cerr << "Error opening " << f_path << std::endl;
		error("Error opening file in write_svtigs()");
	}
	while (getline(fp_read, line))
	{
		if (line[0] == '>')
		{
			if (contig_cnt > 0)
			{
				// continue;
				svtig_name = f_name + "_" + std::to_string(contig_cnt + 1);
			}

			contig_cnt++;
			if (contig == "" && pos == 0 && coverage == 0)
				fp_write << ">" << svtig_name << std::endl;
			else
				fp_write << ">" << svtig_name << " contig=" << contig << " pos=" << pos << " support=" << coverage << std::endl;
		}
		else
			fp_write << line << std::endl;

		line.clear();
	}
	return contig_cnt;
}

int Assembly::merge_svtigs(parameters &params, const std::string &dir)
{
	int cnt = 0;
	std::lock_guard<std::mutex> lk(mtx);
	for (const auto &entry : std::filesystem::directory_iterator(dir))
	{
		if (entry.path().extension() == ".fa" && entry.path().stem().extension() == ".cns")
		{
			if (entry.file_size() == 0)
			{
				unassembled_cnt++;
				continue;
			}
			std::string file_path = entry.path();
			std::string file_name = entry.path().stem().stem();
			std::string tmp = "";
			cnt += write_svtigs(file_path, file_name, 0, tmp, 0, params.fp_svtigs.out());

			if (params.debug)
			{
				std::error_code ec;
				std::filesystem::copy_file(file_path, params.log_path + "out/" + entry.path().filename().string(), std::filesystem::copy_options::overwrite_existing, ec);
			}
		}
	}
	return cnt;
}

// wtdbg2 can write a non-empty but unreadable .ctg.lay.gz, which wtpoa-cns
// segfaults on instead of failing cleanly.
static bool layout_has_contigs(const std::string& layout_gz)
{
    gzFile fp = gzopen(layout_gz.c_str(), "rb");
    if (!fp)
        return false;

    char buf[256];
    const char* line = gzgets(fp, buf, sizeof(buf));
    bool ok = (line != nullptr && buf[0] == '>');
    gzclose(fp);
    return ok;
}

// Assemble SV clusters
int Assembly::final_assembly(parameters& params, faidx_t*& fasta_index,
                             std::set<std::string>& read_set,
                             std::string& svtig_name,
                             double& contig_depth,
                             SVCluster*& sv,
                             std::map<std::string, SVtig*>& final_svtigs, int threads)
{
    unsigned int support_threshold = static_cast<unsigned int>(params.support) / 2;
    if (support_threshold < 3) support_threshold = 3;
    int svtig_tmp_cnt = 0;

    // Whether a cluster is plausible for the region's coverage is a property of
    // the cluster, so the two coverage checks weigh all of its reads. Splitting
    // them by haplotype must not change that verdict; only the support check
    // below is per haplotype. Without a phase file the untagged set holds every
    // read, so this is the read set itself.
    size_t cluster_reads = (sv != nullptr)
        ? sv->reads_h1.size() + sv->reads_h2.size() + sv->reads_untagged.size()
        : read_set.size();

    double zscore = (contig_depth > 0) ? (static_cast<double>(read_set.size()) - contig_depth) / std::sqrt(contig_depth) : 0.0;

    if (contig_depth * 2 < cluster_reads)
    {
        { std::lock_guard<std::mutex> lk(mtx); this->filter_hicov++; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open())
            params.fp_asm_log << svtig_name << "\tFILTERED\treason=high_coverage\treads=" << read_set.size() << "\tcluster_reads=" << cluster_reads << "\tcontig_depth=" << contig_depth << "\tzscore=" << zscore << "\n";
        return 0;
    }
    else if (contig_depth > 5 * cluster_reads)
    {
        { std::lock_guard<std::mutex> lk(mtx); this->filter_lowcov++; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open())
            params.fp_asm_log << svtig_name << "\tFILTERED\treason=low_coverage\treads=" << read_set.size() << "\tcluster_reads=" << cluster_reads << "\tcontig_depth=" << contig_depth << "\tzscore=" << zscore << "\n";
        return 0;
    }
    else if (read_set.size() < support_threshold)
    {
        { std::lock_guard<std::mutex> lk(mtx); this->filter_support++; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open())
            params.fp_asm_log << svtig_name << "\tFILTERED\treason=low_support\treads=" << read_set.size() << "\tthreshold=" << support_threshold << "\tzscore=" << zscore << "\n";
        return 0;
    }
    else if (contig_depth > MAX_CONTIG_DEPTH)
    {
        { std::lock_guard<std::mutex> lk(mtx); this->filter_hicov++; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open())
            params.fp_asm_log << svtig_name << "\tFILTERED\treason=max_contig_depth\treads=" << read_set.size() << "\tcontig_depth=" << contig_depth << "\tzscore=" << zscore << "\n";
        return 0;
    }

    if (sv != nullptr)
    {
        SVtig* tmp = new SVtig;
        tmp->name = svtig_name;
        tmp->pos = sv->ref_pos;
        (tmp->reads).insert(read_set.begin(), read_set.end());
        tmp->contig = sv->contig;
        std::lock_guard<std::mutex> lk(mtx);
        final_svtigs.insert(std::pair<std::string, SVtig*>(tmp->name, tmp));
    }

    std::string job_dir     = params.log_path + "tmp/" + svtig_name + "/";
    std::filesystem::create_directories(job_dir);
    std::string file_path   = job_dir + svtig_name + ".fasta";
    std::string output_path = job_dir + svtig_name;
    JobDir guard{job_dir, file_path, params.debug ? params.log_path + "in/" + svtig_name + ".fasta" : ""};

    generate_fasta_file(params, fasta_index, read_set, file_path);

    int var_size = 4;

    // wtdbg2 -> wtpoa-cns (raw) -> minimap2 | samtools sort -> wtpoa-cns (polish)
    std::string wtdbg2_bin = tool_bin("wtdbg2");
    std::string wtpoa_bin = tool_bin("wtpoa-cns");
    std::string minimap2_bin = tool_bin("minimap2");
    std::string samtools_bin = tool_bin("samtools");

    if (wtdbg2_bin.empty() || wtpoa_bin.empty() ||
        minimap2_bin.empty() || samtools_bin.empty())
    {
        std::string msg = "[final_assembly] Failed to find required tools: ";
        if (wtdbg2_bin.empty())   msg += "wtdbg2 ";
        if (wtpoa_bin.empty())    msg += "wtpoa-cns ";
        if (minimap2_bin.empty()) msg += "minimap2 ";
        if (samtools_bin.empty()) msg += "samtools ";
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open()) params.fp_logs << msg << std::endl;
        error(msg.c_str());
    }

    if (threads < 0)
        threads = params.threads;
#if defined(__APPLE__) && defined(__aarch64__)
    // wtdbg2 + sse2neon crashes with multiple threads on macOS ARM
    int wtdbg2_threads = 1;
#else
    int wtdbg2_threads = threads;
#endif
    int smt_threads = (threads < 1 ? 1 : threads);
    std::string genome_opt = std::to_string(var_size) + "m";

    // Set presets based on read type
    std::string wtdbg2_preset, minimap2_preset;
    if (params.read_type == "hifi") {
        wtdbg2_preset = "ccs";
        minimap2_preset = "map-hifi";
    } else if (params.read_type == "clr") {
        wtdbg2_preset = "rs";
        minimap2_preset = "map-pb";
    } else {
        wtdbg2_preset = "ont";
        minimap2_preset = "map-ont";
    }

    std::string layout_gz = output_path + ".ctg.lay.gz";
    std::string raw_fa    = output_path + ".raw.fa";
    std::string bam_path  = output_path + ".bam";
    std::string cns_fa    = output_path + ".cns.fa";

    // In debug mode, capture stderr; otherwise discard.
    //
    // wtpoa-cns aborts inside glibc on some clusters. The shell that runs the step
    // announces that abort on its OWN stderr, and glibc writes its heap report
    // straight to the terminal, so neither obeys a redirection written inside the
    // command and both surface as if SVarp had crashed. Redirecting the shell
    // itself with exec covers the announcement, and LIBC_FATAL_STDERR_ sends the
    // glibc report along the same path.
    std::string stderr_file = output_path + ".stderr";
    std::string err_target  = params.debug ? stderr_file : std::string("/dev/null");
    std::string quiet       = "exec 2>" + err_target + "; LIBC_FATAL_STDERR_=1 ";
    std::string out_redir   = params.debug ? std::string("") : std::string(" >/dev/null");

    auto asm_t1 = std::chrono::steady_clock::now();

    // Per-step wall-clock caps (seconds): wtpoa-cns can hang for a day on one bad
    // cluster. run_and_log kills the step and reports exit 124, caught below.
    const int TIMEOUT_WTDBG2    = 600;
    const int TIMEOUT_WTPOA_RAW = 120;
    const int TIMEOUT_MM2_SORT  = 300;
    const int TIMEOUT_WTPOA_CNS = 120;

    // 1) wtdbg2 assembler
    std::string asm_cmd = quiet + wtdbg2_bin +
        std::string(" -t ") + std::to_string(wtdbg2_threads) +
        " -x " + wtdbg2_preset +
        " -g " + genome_opt +
        " -fo " + output_path +
        " -i " + file_path +
        out_redir;

    int rc = run_and_log(asm_cmd, params, "wtdbg2_asm", 0, 1, false, TIMEOUT_WTDBG2);

    // A non-zero rc is a real failure and is reported per cluster. A clean exit
    // with no usable layout only means too few or too short reads, so it is
    // counted and reported once in the summary.
    if (rc != 0)
    {
        const char* reason = (WEXITSTATUS(rc) == 124) ? " (timeout)" : "";
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
            params.fp_logs << "[warning] wtdbg2 assembly failed for "
                           << svtig_name << " (rc=" << rc << ")"
                           << reason << std::endl;
        { std::lock_guard<std::mutex> lk(g_log_mtx); std::cout << "[warning] wtdbg2 assembly failed for "
                  << svtig_name << reason << std::endl; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open()) {
            params.fp_asm_log << svtig_name << "\tFAILED\tstep=wtdbg2\trc=" << rc
                              << reason
                              << "\treads=" << read_set.size() << "\n";
            if (std::filesystem::exists(stderr_file)) {
                std::ifstream ef(stderr_file); std::string el;
                while (std::getline(ef, el))
                    params.fp_asm_log << "  stderr: " << el << "\n";
            }
        }
        return 0;
    }

    if (!std::filesystem::exists(layout_gz) ||
        std::filesystem::file_size(layout_gz) == 0 ||
        !layout_has_contigs(layout_gz))
    {
        { std::lock_guard<std::mutex> lk(mtx); this->no_contig_cnt++; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open())
            params.fp_asm_log << svtig_name << "\tFAILED\tstep=wtdbg2\treason=no_contig"
                              << "\treads=" << read_set.size() << "\n";
        return 0;
    }

    // 2) raw consensus
    std::string cns_raw_cmd = quiet + wtpoa_bin +
        std::string(" -t ") + std::to_string(threads) +
        " -i " + layout_gz +
        " -fo " + raw_fa +
        out_redir;

    rc = run_and_log(cns_raw_cmd, params, "wtpoa_raw", 0, 1, false, TIMEOUT_WTPOA_RAW);
    if (rc != 0 || !std::filesystem::exists(raw_fa) || std::filesystem::file_size(raw_fa) == 0)
    {
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
            params.fp_logs << "[warning] wtpoa-cns raw consensus failed for "
                           << svtig_name << " (rc=" << rc << ")" << std::endl;
        { std::lock_guard<std::mutex> lk(g_log_mtx); std::cout << "[warning] wtpoa-cns raw consensus failed for "
                  << svtig_name << std::endl; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open()) {
            params.fp_asm_log << svtig_name << "\tFAILED\tstep=wtpoa-raw\trc=" << rc
                              << "\treads=" << read_set.size() << "\n";
            if (std::filesystem::exists(stderr_file)) {
                std::ifstream ef(stderr_file); std::string el;
                while (std::getline(ef, el))
                    params.fp_asm_log << "  stderr: " << el << "\n";
            }
        }
        return 0;
    }

    // 3) minimap2 + samtools sort
    // Wrapped in sh -c so LIBC_FATAL_STDERR_ reaches both commands of the pipeline.
    std::string inner_mm =
        minimap2_bin +
        " -ax " + minimap2_preset +
        " -t" + std::to_string(threads) +
        " -r2k " + raw_fa + " " + file_path +
        " | " +
        samtools_bin +
        " sort -m 512m -@" + std::to_string(smt_threads) +
        " -o " + bam_path;
    std::string map_sort_cmd =
        quiet + "sh -c \"" + inner_mm + "\"" + out_redir;

    rc = run_and_log(map_sort_cmd, params, "mm2_samtools", 0, 1, false, TIMEOUT_MM2_SORT);
    if (rc != 0 || !std::filesystem::exists(bam_path))
    {
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
            params.fp_logs << "[warning] minimap2+samtools sort failed for "
                           << svtig_name << " (rc=" << rc << ")" << std::endl;
        { std::lock_guard<std::mutex> lk(g_log_mtx); std::cout << "[warning] minimap2+samtools sort failed for "
                  << svtig_name << std::endl; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open())
            params.fp_asm_log << svtig_name << "\tFAILED\tstep=mm2+samtools\trc=" << rc
                              << "\treads=" << read_set.size() << "\n";
        return 0;
    }

    // 4) polishing consensus
    std::string inner_pol =
        samtools_bin +
        " view -F0x900 " + bam_path +
        " | " +
        wtpoa_bin +
        " -t " + std::to_string(threads) +
        " -d " + raw_fa +
        " -i - -fo " + cns_fa;
    std::string polish_cmd =
        quiet + "sh -c \"" + inner_pol + "\"" + out_redir;

    rc = run_and_log(polish_cmd, params, "wtpoa_cns_polish", 0, 1, false, TIMEOUT_WTPOA_CNS);

    if (rc != 0 || !std::filesystem::exists(cns_fa) || std::filesystem::file_size(cns_fa) == 0)
    {
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
            params.fp_logs << "[warning] wtdbg2 pipeline failed for "
                           << svtig_name << " (rc=" << rc << ", empty="
                           << (std::filesystem::exists(cns_fa) && std::filesystem::file_size(cns_fa) == 0)
                           << ")" << std::endl;
        { std::lock_guard<std::mutex> lk(g_log_mtx); std::cout << "[warning] wtdbg2 pipeline failed for "
                  << svtig_name << std::endl; }
        if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open())
            params.fp_asm_log << svtig_name << "\tFAILED\tstep=polish\trc=" << rc
                              << "\treads=" << read_set.size() << "\n";
        return 0;
    }

    auto asm_t2 = std::chrono::steady_clock::now();
    if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_asm_log.is_open()) {
        std::string contig_name = (sv != nullptr) ? sv->contig : "unknown";
        int pos = (sv != nullptr) ? sv->ref_pos : 0;
        params.fp_asm_log << svtig_name << "\tOK"
                          << "\treads=" << read_set.size()
                          << "\tcontig=" << contig_name
                          << "\tpos=" << pos
                          << "\tzscore=" << zscore
                          << "\ttime=" << format_duration(std::chrono::duration<double>(asm_t2 - asm_t1).count())
                          << "\n";
    }
	
    svtig_tmp_cnt = merge_svtigs(params, job_dir);

    return svtig_tmp_cnt;
}


// Expected depth for a cluster. Reference nodes (SR 0) use their contig's depth;
// alt nodes use the genome-wide depth, since a donor contig's own depth only
// counts reads on its nodes. Without SR tags a contig at genome depth is taken
// as reference.
double Assembly::cluster_depth(const SVCluster* sv, std::map <std::string, Contig*>& depth)
{
    double overall = depth.count("overall") ? depth["overall"]->coverage : 0.0;
    double contig = depth.count(sv->contig) ? depth[sv->contig]->coverage : 0.0;
    double lambda = overall;
    if (sv->rank == 0 || (sv->rank < 0 && contig >= 0.5 * overall))
        lambda = contig;
    if (lambda < 5)
        lambda = 5;
    return lambda;
}

void Assembly::run_assembly(parameters &params, std::map<std::string, Contig *> &depth, std::map<std::string, std::vector<SVCluster *>> &vars, std::map<std::string, SVtig *> &final_svtigs)
{
	int initial_svtigs_cnt = 0;
	std::map<std::string, std::vector<SVCluster *>>::iterator itr;

	auto t1 = std::chrono::steady_clock::now();
	std::cout << "\nAssembly..." << std::endl;
	std::cout << "--> assembling reads using " << params.assembler << std::endl;
	if (params.fp_logs.is_open())
	{
		params.fp_logs << "--> wtdbg2 " << tool_version(tool_bin("wtdbg2"), "-V") << "\n";
		params.fp_logs << "--> wtpoa-cns " << tool_version(tool_bin("wtpoa-cns"), "-V") << "\n";
		params.fp_logs << "--> minimap2 " << tool_version(tool_bin("minimap2")) << "\n";
		params.fp_logs << "--> samtools " << tool_version(tool_bin("samtools")) << "\n";
	}

	std::string svtigs_tmp_path = params.log_path + params.sample_name + "_svtigs_tmp.fa";
	params.fp_svtigs.open(svtigs_tmp_path);

	faidx_t *fasta_index = fai_load((params.fasta).c_str());
	if (!fasta_index)
		error("Error loading FASTA index: file not found or corrupted");

	// One job per cluster and haplotype; clusters are assembled in parallel,
	// each job running the external tools with threads/asm_jobs threads.
	struct Job { SVCluster* sv; std::set<std::string>* reads; std::string name; double depth; };
	std::vector<Job> jobs;
	for (itr = vars.begin(); itr != vars.end(); ++itr)
		for (auto &sv : itr->second)
		{
			double d = cluster_depth(sv, depth);
			std::string base = sv->node + "_" + std::to_string(sv->start_pos);
			if ((params.phase_tags).empty())
				jobs.push_back({sv, &sv->reads_untagged, base, d});
			else
			{
				jobs.push_back({sv, &sv->reads_h1, "H1-" + base, d});
				jobs.push_back({sv, &sv->reads_h2, "H2-" + base, d});
				if (!params.skip_untagged)
					jobs.push_back({sv, &sv->reads_untagged, "None-" + base, d});
			}
		}
	int n_jobs = params.asm_jobs > 0 ? params.asm_jobs : 1;
	if (n_jobs > (int) jobs.size()) n_jobs = jobs.size() > 0 ? jobs.size() : 1;
	int job_threads = params.threads / n_jobs;
	if (job_threads < 1) job_threads = 1;
	{
		std::lock_guard<std::mutex> lk(g_log_mtx);
		std::cout << "--> " << jobs.size() << " assembly jobs, " << n_jobs << " in parallel, " << job_threads << " thread(s) each" << std::endl;
	}

	std::atomic<size_t> next{0};
	std::atomic<int> done{0};
	std::atomic<int> produced{0};
	auto worker = [&]() {
		faidx_t *idx = fai_load((params.fasta).c_str());
		if (!idx)
			error("Error loading FASTA index: file not found or corrupted");
		while (true)
		{
			size_t i = next.fetch_add(1);
			if (i >= jobs.size()) break;
			Job &j = jobs[i];
			produced += final_assembly(params, idx, *j.reads, j.name, j.depth, j.sv, final_svtigs, job_threads);
			int k = ++done;
			if (k % 100 == 0)
			{
				std::lock_guard<std::mutex> lk(g_log_mtx);
				std::cout << "\r--> assembled " << k << "/" << jobs.size() << " jobs" << std::flush;
			}
		}
		fai_destroy(idx);
	};
	std::vector<std::thread> pool;
	for (int t = 0; t < n_jobs; t++) pool.emplace_back(worker);
	for (auto &t : pool) t.join();
	initial_svtigs_cnt += produced;
	if (jobs.size() >= 100)
		std::cout << "\r--> assembled " << jobs.size() << "/" << jobs.size() << " jobs\n";

	if (std::filesystem::exists(params.log_path + "tmp/"))
		std::filesystem::remove_all(params.log_path + "tmp/");

	if (params.fp_svtigs.is_open())
		params.fp_svtigs.close();
	fai_destroy(fasta_index);

	std::cout << "--> " << (filter_hicov) + (this->filter_lowcov) + (this->filter_support) << " filtered (" << this->filter_hicov << " high, " << this->filter_lowcov << " low coverage read clusters and " << this->filter_support << " low read support)\n";

	if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
		params.fp_logs << "--> " << (this->filter_hicov) + (this->filter_lowcov) + (this->filter_support) << " filtered (" << this->filter_hicov << " high, " << this->filter_lowcov << " low coverage read clusters and " << this->filter_support << " low read support)\n";

	std::cout << "--> " << this->no_contig_cnt << " clusters produced no contig (too few or too short reads)\n";
	if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
		params.fp_logs << "--> " << this->no_contig_cnt << " clusters produced no contig (too few or too short reads)\n";

	std::cout << "--> " << unassembled_cnt << " clusters cannot be assembled\n";
	if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
		params.fp_logs << "--> " << unassembled_cnt << " clusters cannot be assembled\n";

	std::cout << "--> " << initial_svtigs_cnt << " svtigs before final filtering\n";
	if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
		params.fp_logs << "--> " << initial_svtigs_cnt << " svtigs before final filtering\n";

	auto t2 = std::chrono::steady_clock::now();
	std::string asm_dur = format_duration(std::chrono::duration<double>(t2 - t1).count());

	std::cout << "--> assembly execution time: " << asm_dur << "\n";
	if (std::lock_guard<std::mutex> lk(g_log_mtx); params.fp_logs.is_open())
		params.fp_logs << "--> assembly execution time: " << asm_dur << "\n";
}
