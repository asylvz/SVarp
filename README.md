# SVarp
**Pangenome-based structural variant discovery**

SVarp discovers **haplotype-resolved structural variants (SVs)** on pangenome graphs using long-read sequencing data. It outputs **local assemblies** of SV alleles (*svtigs*).

*For questions, please open an issue.*

---

## Quick Start

### **Bioconda (Linux only)**
```bash
conda create -n svarp -c conda-forge -c bioconda svarp
conda activate svarp
svarp --help
```

### **Building from source (Linux and macOS)**
#### 1. Clone the repository

```bash
git clone https://github.com/asylvz/SVarp.git
cd SVarp
```
#### 2. Install third-party dependencies
SVarp depends on:
* **zlib** *(sudo apt-get install zlib1g-dev)*
* **GraphAligner** *(https://github.com/maickrau/GraphAligner)*
* **Samtools** *(https://github.com/samtools/samtools)*
* **HTSlib** *(https://github.com/samtools/htslib/)*
* **WFA2-lib** *(https://github.com/smarco/WFA2-lib)*
* **Minimap2** *(https://github.com/lh3/minimap2)*

> **Note:** wtdbg2, Minimap2, and Samtools are bundled under `third_party/` and built automatically. Their paths are resolved relative to the `svarp` executable, so the binary is portable.

You can download and build **HTSlib, WFA2-lib, wtdbg2, Minimap2, Samtools** using:
```bash
make libs
```

#### macOS (Apple Silicon) notes
SVarp is developed and tested on Linux. It compiles and runs on macOS ARM,
but has not been fully tested due to GraphAligner lacking native ARM support.
`gcc` is required to build the bundled wtdbg2 (`brew install gcc`).

#### 3. Build SVarp
```bash
make
```
>Binary will be at **build/svarp**.

#### 4. Run

**ONT (default):**
```bash
build/svarp \
    --gaf sample.gaf \
    --graph pangenome.gfa \
    --fasta reads.fasta.gz \
    --sample SAMPLE1 \
    --out output_dir
```

**PacBio HiFi:**
```bash
build/svarp \
    --gaf sample.gaf \
    --graph pangenome.gfa \
    --fasta reads.fasta.gz \
    --reads hifi \
    --sample SAMPLE1 \
    --out output_dir
```

## Input Requirements

* Reads FASTA must be bgzipped (i.e., `bgzip NA12878.fasta`)
* GAF must follow the unstable rGFA coordinate system (https://github.com/lh3/gfatools/blob/master/doc/rGFA.md)
* Graph must be provided in GFA format. We tested SVarp with **HPRC Minigraph pangenome** (https://zenodo.org/record/6983934)
	- GAF alignments must be generated using the same GFA file that you use as input
* Phasing tags (optional) should be generated using WhatsHap: (https://whatshap.readthedocs.io/)
```bash
whatshap haplotag input.vcf.gz input.bam \
    --reference ref.fasta \
    --output-haplotag-list tags.tsv \
    -o phased.bam
```

#### Example Usage
```bash
build/svarp \
    --gaf sample.gaf \
    --graph pangenome.gfa \
    --fasta reads.fasta.gz \
    --reads hifi \
    --phase read_tags.tsv \
    --sample SAMPLE1 \
    --out output_dir
```

## All parameters

	Required arguments:
	--gaf (-a)                  : Alignment file in GAF format
	--graph (-g)                : Pangenome file in GFA format
	--fasta (-f)                : Fasta sequence file


	Optional arguments:
	--sample (-i)               : Sample name.
	--out (-o)                  : Output folder.
	--phase (-p)                : WhatsHap haplotag file in .tsv (https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
	--support (-s)              : A cluster is assembled only if it has MORE than this many reads (default=5 for diploid samples)
	--dist-threshold (-d)       : Distance threshold to merge SV breakpoints (default=100)
	--reads (-w)                : Read type: ont (default), hifi, or clr. Sets wtdbg2 preset and assembly parameters.
	--threads (-t)              : Number of threads for assembly and realignment (default=16)
	--asm-jobs                  : Clusters assembled in parallel (default=min(threads, 8))
	--keep-untagged             : Also assemble and write clusters of untagged reads (off by default)
	--min-clip                  : Unaligned read end (bp) that counts as a breakpoint on a single alignment (default=500)
	--min-graph-cov             : Fraction of an svtig covered by graph alignments (identity >= --min-identity) needed to keep it (default=0.90)
	--min-svtig-len             : Shortest svtig written, in bp (default=5000)
	--min-identity              : Remap records below this identity do not count as graph coverage (default=0.90)
	--as                        : GraphAligner minimum alignment score for remapping (default=1000)
	--pc                        : GraphAligner --precise-clipping for remapping (default: GraphAligner default)
	--no-remap (-r)             : Skip remapping (not suggested)
	--keep-remap                : Keep <sample>_svtigs_tmp.fa and <sample>_remap.gaf (the inputs of the final filter)
	--write-unmapped            : Write the names of reads without a GAF record to <out>/<sample>_unmapped_reads.txt
	--debug (-u)                : Keep all intermediate files and log every command
	--version (-v)              : Print version
	--help (-h)                 : Print this help menu

## Output

Svtigs are written per haplotype (`<sample>_svtigs_H1.fa`, `_H2.fa`; `_untagged.fa` with `--keep-untagged`).
Each header carries the cluster locus and read support plus the result of remapping the svtig onto the graph:

	>H1-s1065813_569 contig=CHM13#0#chr16 pos=4106 support=6 path=>s1065813>s1065814 graph_cov=0.998 max_gap=0 max_indel=12 graph_explained=yes

`graph_cov` is the fraction of the svtig covered by graph alignments, `max_gap` the largest uncovered stretch and
`max_indel` the largest indel inside an alignment (bp). `graph_explained=yes` means the graph already holds this
sequence without an SV-sized (>=50 bp) difference; `no` marks an allele the graph lacks. Only svtigs anchored in the
graph (`graph_cov` >= `--min-graph-cov`) and at least `--min-svtig-len` long are written; every dropped svtig is listed
with its reason in `<sample>_remap.log`, every cluster in `<sample>_assembly.log`.

## Citation

*Soylev, A., Ebler, J., Pani, S., Rausch, T., Korbel, J., & Marschall, T. (2024). **SVarp: pangenome-based structural variant discovery**. bioRxiv, 2024-02.*
