# SVarp 
**Pangenome-based structural variant discovery**

SVarp discovers **haplotype-resolved structural variants (SVs)** on pangenome graphs using long-read sequencing data. It outputs **local assemblies** of SV alleles (*svtigs*).

*For questions, please open an issue or feel free to send me an e-mail (asoylev@gmail.com)*

---

## Quick Start

### **Bioconda (recommended)**  
```bash
conda create -n svarp -c conda-forge -c bioconda svarp
conda activate svarp
svarp --help
```
>Note: SVarp currently supports 64-bit Linux only,


### **Building from source** 
#### 1. Clone the repository


```bash
git clone https://github.com/asylvz/SVarp.git
cd svarp
```
#### 2. Install third-party dependencies
SVarp depends on:
* **zlib** *(sudo apt-get install zlib1g-d)*
* **GraphAligner** *(https://github.com/maickrau/GraphAligner)*
* **Samtools** *(https://github.com/samtools/samtools)*
* **HTSlib** *(https://github.com/samtools/htslib/)*
* **WFA2-lib** *(https://github.com/smarco/WFA2-lib)*
* **wtdbg2** *(https://github.com/ruanjue/wtdbg2)*
* **Minimap2** *(https://github.com/lh3/minimap2)*

You can download **HTSlib, WFA2-lib, wtdbg2, Minimap2** using:
```bash
make libs
```
#### 3.Build SVarp 
```bash
make
```
>Binary will be at **build/svarp**. 

#### 4. Run
```bash
build/svarp \
    --gaf sample.gaf \
    --graph pangenome.gfa \
    --fasta reads.fasta.gz \
    --sample SAMPLE1 \
    --out output_dir
```

## Input Requirements

SVarp is developed and tested using Linux Ubuntu operating system

* Reads FASTA must be bgzipped (i.e., bgzip NA12878.fasta)
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
	--debug                     : Output multiple log files for debugging purpose.
	--skip-untagged             : Output only phased variants (~30% faster).
	--dist_threshold (-d)       : Distance threshold to merge SV breakpoints (default=100)
	--out (-o)                  : Output folder path
	--phase (-p)                : WhatsHap haplotag file in .tsv (https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
	--reads(-r)                 : Bgzipped FASTA file of reads for extensive mode (needed for WFA realignment)
	--sample (-i)               : Sample (Individual) name
	--support (-s)              : Minimum support for a cluster to be assembled (default=5 for diploid samples)
	--threads(-t)               : Number of threads for assembly and realignment (default:32)
	--help                      : Print this help menu

## Citation

*Soylev, A., Ebler, J., Pani, S., Rausch, T., Korbel, J., & Marschall, T. (2024). **SVarp: pangenome-based structural variant discovery**. bioRxiv, 2024-02.*
